import time
from datetime import datetime, timedelta

import logging
logger = logging.getLogger(__name__)
#from goreverselookup import logger


class Timer:
    def __init__(self, millisecond_init=False):
        """
        If you will be measuring function execution times, use millisecond_init.
        """
        if millisecond_init is False:
            self.start_time = time.time()
        else:
            self.start_time = time.time() * 1000
        self.elapsed_time = 0.0

        if millisecond_init is False:
            self.init_type = "seconds"
        else:
            self.init_type = "milliseconds"

    def set_start_time(self):
        """
        Sets a new reference start time.
        """
        self.start_time = time.time()

    def get_elapsed(self, format: str = "seconds") -> int:
        """
        Returns the amount of seconds unformatted (contains decimal places)

        'format' can be "seconds" or "milliseconds"
        """
        if format == "seconds":
            if self.init_type == "milliseconds":
                self.start_time = self.start_time / 1000  # prevent wrong init values
            return time.time() - self.start_time
        elif format == "milliseconds":
            if self.init_type == "seconds":
                self.start_time = self.start_time * 1000  # prevent wrong init values
            return time.time() * 1000 - self.start_time

    def get_elapsed_formatted(
        self, format: str = "seconds", reset_start_time: bool = False
    ) -> str:
        """
        If format is 'seconds':
            Gets elapsed time in hh mm ss format
        If format is 'milliseconds':
            Gets elapsed time in ss:ms format

        'format' can be "seconds" or "milliseconds"
        If 'reset_start_time', self.start_time will be reset to the time of calling this function.
        """
        return_value = ""
        if format == "seconds":
            sec = int(self.get_elapsed("seconds"))
            td = timedelta(seconds=sec)
            return_value = str(td)
        elif format == "milliseconds":
            elapsed_time_ms = int(self.get_elapsed("milliseconds"))
            seconds, milliseconds = divmod(elapsed_time_ms, 1000)
            formatted_time = f"{seconds}s:{milliseconds:03d}ms"
            return_value = formatted_time

        if reset_start_time:
            if format == "seconds":
                self.start_time = time.time()
            else:
                self.start_time = time.time() * 1000

        return return_value

    def print_elapsed_time(self, useLogger: bool = True, prefix: str = "Elapsed: "):
        """
        Prints the elapsed time in hh mm ss format.

        Args:
          - useLogger: if True, then logger.info is used. If false, then print is used.
          - prefix: the string you want to use as a prefix
        """
        if useLogger:
            logger.info(f"{prefix}{self.get_elapsed_formatted()}")
        else:
            print(f"{prefix}{self.get_elapsed_formatted()}")

    @classmethod
    def get_current_time(cls):
        """
        Gets the current time and returns it in the format "%Y-%m-%d %H:%M:%S"
        """
        current_time = datetime.now()
        formatted_time = current_time.strftime("%Y-%m-%d %H:%M:%S")
        return formatted_time

    @classmethod
    def compare_time(cls, timestamp_one: str, timestamp_two: str) -> bool:
        """
        Compares timestamp_two against timestamp_one. If timestamp_two is greater than timestamp_one
        (aka timestamp_two was recorded at a time later than timestamp_one), the function returns True.

        The input timestamps must be supplied in the format "%Y-%m-%d %H:%M:%S"
        """
        format_str = "%Y-%m-%d %H:%M:%S"
        time_one = datetime.strptime(timestamp_one, format_str)
        time_two = datetime.strptime(timestamp_two, format_str)

        return time_two > time_one
