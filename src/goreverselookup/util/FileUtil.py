import os
import requests
import shutil
import urllib
import gzip

import logging

# from logging import config
# config.fileConfig("../logging_config.py")
logger = logging.getLogger(__name__)


class FileUtil:
    project_root_path = ""

    def __init__(self, root=""):
        """
        A file utility class. The construction initializes project_root_path to the project's root directory, or sets it to the value
        of 'root'.
        """
        if root == "":
            self.project_root_path = os.path.dirname(os.path.abspath(__file__))
        else:
            self.project_root_path = root
        logger.info(f"Project root path set to: {self.project_root_path}")

    def find_file(self, filepath: str, backtrace: int = 6):
        """
        Attempts a search for a relative filepath, by conducting a backwards-folder search. This function, unlike find_win_abs_filepath,
        cannot be statically called, meaning find_file(...) must be called on an instantiated FileUtil instance.

        Parameters:
          - filepath: a relative filepath to find. Can be a single file (eg. data.json) or a file in a folder(s) eg. program_data/data/data.json
          - backtrace: the amount of levels to scan backwards

        Calling:
            from .FileUtil import FileUtil
            fu = new FileUtil()
            fu.find_file(FILEPATH)

        Returns the correct path to filepath, if it is found during the amount of 'backtraces' performed on the root filepath computed during init of this class
        by os.path.dirname(os.path.abspath(__file__)).

        A better alternative for finding a relative path to a file on windows is by using (FileUtil).find_win_abs_path(...)
        """

        def check_contains(folder_path, element_name):
            """
            Checks if 'filename' is in 'folder_path'.

            Parameters:
              - folder_path: the path to the folder where to search
              - element_name: either a filename or a folder name

            Return True if found or False if not found
            """
            folder_elements = os.listdir(folder_path)
            for folder_element in folder_elements:
                if folder_element == element_name:
                    return True
            return False

        # bugfix: sometimes when using os.path.join("folder1/folder2", "file1"), the result can be: folder1/folder2\file1
        if "\\" in filepath and "/" in filepath:
            filepath = filepath.replace("\\", "/")

        # first, see if is a single file or a file inside folder(s)
        folders = []
        file = ""
        if os.sep in filepath:
            folders = filepath.split(os.sep)[
                :-1
            ]  # folders are all but the last element
            file = filepath.split(os.sep)[-1]  # file is the last in filepath
        elif "/" in filepath:
            folders = filepath.split("/")[:-1]  # folders are all but the last element
            file = filepath.split("/")[-1]  # file is the last in filepath
        else:
            file = filepath

        current_path = self.project_root_path
        if len(folders) == 0:  # filepath is a single file
            for i in range(
                backtrace
            ):  # ascend 'backtrace' levels to the parent directory
                # we ascend a directory using os.path.dirname on a filepath
                if i != 0:  # don't ascend up 1 directory on the first iteration
                    parent_path = os.path.dirname(current_path)
                    current_path = parent_path
                elif i == 0:  # stay on the root for the first iteration
                    parent_path = current_path
                if check_contains(parent_path, file):
                    # file was found in parent_path
                    base_path = self.project_root_path
                    for j in range(i):  # append i number of backtraces
                        base_path = os.path.join(base_path, "..")
                    return os.path.join(base_path, file)  # return the correct path

        current_path = self.project_root_path  # reset
        if len(folders) != 0:  # filepath contains folders, primary search is for folder
            for i in range(backtrace):
                if i != 0:
                    parent_path = os.path.dirname(current_path)
                    current_path = parent_path
                elif i == 0:
                    parent_path = current_path
                if check_contains(parent_path, folders[0]):
                    # folder was found in parent path
                    base_path = self.project_root_path
                    for j in range(i):  # append i number of backtraces
                        base_path = os.path.join(base_path, "..")
                    # append the folders
                    for folder in folders:
                        if check_contains(base_path, folder):
                            base_path = os.path.join(base_path, folder)
                        else:
                            logger.info(
                                "ERROR during file search. Base path"
                                f" {base_path} doesn't contain folder {folder}"
                            )
                    # finally, append the file
                    return os.path.join(base_path, file)

    @classmethod
    def find_win_abs_filepath(cls, relative_filepath: str):
        """
        Finds the absolute filepath to specified 'relative_filepath' using os.getcwd(). Note, that the function does not include
        the "self" parameter, therefore it can be statically called using FileUtil.find_win_abs_filepath(...) without instantiating
        a FilePath instance.

        Parameters:
          - (str) relative_filepath: a relative filepath from the root of the project's directory to the file

        Calling:
            from .FileUtil import FileUtil
            FileUtil.find_win_abs_filepath(FILEPATH)

        Returns:
          - an absolute filepath to the specified relative filepath: os.getcwd + relative filepath
        """
        if "/" in relative_filepath:
            relative_filepath = relative_filepath.replace("/", os.sep)

        return os.path.join(os.getcwd(), relative_filepath)

    @classmethod
    def get_workspace_dir(cls):
        """
        Returns the workspace-specific directory using os.getcwd().
        """
        return os.getcwd()

    @classmethod
    def write_to_txt_file(cls, list: list, filepath: str):
        """
        Writes the input list, line by line, to the specified filepath.
        """
        filepath = cls.find_win_abs_filepath(filepath)
        with open(filepath, "w") as f:
            for element in list:
                line = element
                if "\n" not in line:
                    line = f"{line}\n"
                f.write(line)

    @classmethod
    def is_file_empty(cls, filepath: str):
        """
        Checks if 'filepath' is an empty file. Returns True if the file is empty.
        """
        return True if os.path.getsize(filepath) == 0 else False

    @classmethod
    def check_paths(cls, paths: list[str], are_files=True, auto_create=True):
        """
        Checks if the directories of the given paths exists. If they don't, this function automatically creates the
        directories (all the folders leading to the path).

        Params:
          - (list[str]) paths: The paths to the files or folders
          - (bool) are_files: If True, paths point to a files
                              If False, paths points to a folders
          - (bool) auto_create: If True, will automatically create the path if it doesn't exist.

        Warning: 'path' is by default presumably pointing to files (specified by 'is_file' set to True).
        If you are passing paths pointing to folders, make sure to set 'is_file' to False.
        """
        for path in paths:
            cls.check_path(path=path, is_file=are_files, auto_create=auto_create)

    @classmethod
    def check_path(cls, path: str, is_file=True, auto_create=True):
        """
        Checks if the directory of the given path exists and create the path.

        Params:
          - (str) path: The path to the file or folder
          - (bool) is_file: If True, path points to a file
                            If False, path points to a folder
          - (bool) auto_create: If True, will automatically create the path if it doesn't exist.

        Warning: 'path' is by default presumably pointing to a file (specified by 'is_file' set to True).
        If you are passing a path pointing to a directory, make sure to set 'is_file' to False.
        """
        if is_file is True:
            dir_path = os.path.dirname(path)
        else:
            dir_path = path

        # check if directory exist - if it doesn't create it
        if not os.path.exists(dir_path) and auto_create is True:
            os.makedirs(dir_path)

        # check if the file exists - if it doesn't, create it
        if is_file is True:
            if not os.path.exists(path) and auto_create is True:
                with open(path, "w"):  # creates an empty file
                    pass

        if os.path.exists(path):
            return True
        else:
            return False

    @classmethod
    def clear_file(cls, filepath: str, replacement_text: str = ""):
        """
        Clears the file contents while preserving the file.

        Params:
          - (str) filepath: the file to clear
          - (str) replacement_text: any optional text to write to the file after the clear
        """
        try:
            with open(filepath, "w") as f:
                f.truncate(0)
                f.write(replacement_text)
                return True
        except FileNotFoundError:
            logger.warning(f"File {filepath} wasn't found!")
    
    @classmethod
    def download_file(cls, filepath:str, download_url:str):
        """
        Downloads the file specified by 'download_url' to 'filepath'.
        
        Warning: For .txt web files, rather use the download_text_file function.
        """
        if os.path.exists(filepath) and not cls.is_file_empty(filepath):
            logger.info(f"File {filepath} exists and isn't empty.")
            return
        # file doesn't exist or is empty -> create filepath tree and download
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        response = requests.get(download_url, stream=True)
        if response.status_code == 200:
            with open(filepath, 'wb') as file:
                shutil.copyfileobj(response.raw, file)
            if not cls.is_file_empty(filepath):
                logger.info(f"Successfully downloaded {download_url} to {filepath}!")
        else:
            raise Exception(f"Failed to download {download_url}! Response code: {response.status_code}, reason: {response.reason}")

    @classmethod
    def download_txt_file(cls, filepath:str, download_url:str):
        """
        Use this function to download .txt files from the web.

        Downloads 'download_url' and saves it into 'filepath'
        """
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        if not os.path.exists(filepath):
            url = download_url
            response = requests.get(url)
            with open(filepath, "wb") as f:
                f.write(response.content)
        if not cls.is_file_empty(filepath):
            logger.info(f"Successfully downloaded {download_url} to {filepath}")
    
    @classmethod
    def download_zip_file(cls, filepath:str, download_url:str, zip_specifier:str="rt"):
        """
        Downloads a zip file from 'download_url' into 'filepath'.
        Use the appropriate 'zip_specifier'.

        Example usage:
        download_url = "http://geneontology.org/gene-associations/goa_human.gaf.gz"
        download_path = "data_files/goa_human_test.gaf"
        FileUtil.download_zip_file(download_path, download_url, "rt")
        """
        temp_file, _ = urllib.request.urlretrieve(download_url)

        # read the contents of the gzip file and save it to the txt file
        with gzip.open(temp_file, "rt") as f_in, open(filepath, "w") as f_out:
            for line in f_in:
                f_out.write(line)

        # delete the temporary file
        os.remove(temp_file)