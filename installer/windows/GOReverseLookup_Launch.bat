@echo off
REM Get the path of the current directory (where the .bat file is)
SET "INSTALL_DIR=%~dp0"
SET "LOG_FILE=%INSTALL_DIR%launch_debug.log"

REM --------------------------------------------------------------------------
REM --- START LOGGING ---
REM Clear the log file and record the start time
> "%LOG_FILE%" ECHO [DEBUG START] %DATE% %TIME%
ECHO -------------------------------------------------- >> "%LOG_FILE%"

REM CRITICAL: Disable Python output buffering for real-time console output
SET PYTHONUNBUFFERED=1
ECHO [CONFIG] PYTHONUNBUFFERED set to: 1 >> "%LOG_FILE%"

REM Log the directory path
ECHO [STEP 1] INSTALL_DIR set to: "%INSTALL_DIR%" >> "%LOG_FILE%"

REM Change directory to the installation folder
cd /d "%INSTALL_DIR%"
IF ERRORLEVEL 1 (
    ECHO [ERROR] Failed to change directory to "%INSTALL_DIR%" >> "%LOG_FILE%"
    GOTO :END_SCRIPT
)
ECHO [STEP 2] Successfully changed current directory to: %CD% >> "%LOG_FILE%"

REM Clear the console window for a clean user interface
CLS

REM Display user messages
ECHO.
ECHO Welcome to the GOReverseLookup Command Line Interface!
ECHO.
ECHO ------------------------------------------------------------------
ECHO Python Output Buffering is DISABLED. Program printouts will appear immediately.
ECHO ------------------------------------------------------------------
ECHO.
ECHO You can now run the application by typing:
ECHO GOReverseLookup.exe ^<arguments^>
ECHO.
ECHO For example:
ECHO GOReverseLookup.exe "input.txt" --full_directory_op "C:/path/to/project/"
ECHO.
ECHO Debug log created at: "%LOG_FILE%"
ECHO.

REM --------------------------------------------------------------------------
REM --- LAUNCH INTERACTIVE SHELL ---
ECHO [STEP 3] Attempting to launch persistent CMD shell... >> "%LOG_FILE%"

REM Execute a new instance of CMD, keeping the window open (/k)
CMD /k

REM --------------------------------------------------------------------------
REM The following lines will ONLY execute if CMD /k somehow failed to launch or if the .bat file continued.
ECHO [ERROR] CMD /k command failed or was skipped. >> "%LOG_FILE%"

:END_SCRIPT
REM Log the end time
ECHO [DEBUG END] %DATE% %TIME% >> "%LOG_FILE%"

REM Use PAUSE here to stop the window from closing if CMD /k failed.
PAUSE