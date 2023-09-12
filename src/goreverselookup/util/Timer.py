import time
from datetime import datetime, timedelta

import logging
#from logging import config
#config.fileConfig("../logging_config.py")
logger = logging.getLogger(__name__)

class Timer:
    def __init__(self):
        self.start_time = time.time()
    
    def set_start_time(self):
        """
        Sets a new reference start time.
        """
        self.start_time = time.time()
    
    def get_elapsed_seconds(self) -> int:
        """
        Returns the amount of seconds unformatted (contains decimal places)
        """
        return time.time() - self.start_time
    
    def get_elapsed_time(self) -> str:
        """
        Gets elapsed time in hh mm ss format.
        """
        sec = int(self.get_elapsed_seconds())
        td = timedelta(seconds=sec)
        return str(td)
    
    def print_elapsed_time(self, useLogger: bool = True, prefix: str = "Elapsed: "):
        """
        Prints the elapsed time in hh mm ss format. 
        
        Args:
          - useLogger: if True, then logger.info is used. If false, then print is used.
          - prefix: the string you want to use as a prefix
        """
        if useLogger:
            logger.info(f"{prefix}{self.get_elapsed_time()}")
        else:
            print(f"{prefix}{self.get_elapsed_time()}")
    
    @classmethod
    def get_current_time(cls):
        """
        Gets the current time and returns it in the format "%Y-%m-%d %H:%M:%S"
        """
        current_time = datetime.now()
        formatted_time = current_time.strftime("%Y-%m-%d %H:%M:%S")
        return formatted_time
    
    @classmethod
    def compare_time(cls, timestamp_one:str, timestamp_two:str) -> bool:
        """
        Compares timestamp_two against timestamp_one. If timestamp_two is greater than timestamp_one
        (aka timestamp_two was recorded at a time later than timestamp_one), the function returns True.
        
        The input timestamps must be supplied in the format "%Y-%m-%d %H:%M:%S"
        """
        format_str = "%Y-%m-%d %H:%M:%S"
        time_one = datetime.strptime(timestamp_one, format_str)
        time_two = datetime.strptime(timestamp_two, format_str)
    
        return time_two > time_one