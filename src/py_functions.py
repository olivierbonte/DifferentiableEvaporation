import logging
import os
import subprocess
from datetime import datetime


def run_and_log(command, sys_info=False):
    """Executes a shell command and logs its output. Optionally logs system information.

    Parameters
    ----------
    command : str
        The shell command to execute.
    sys_info : bool, optional
        If True, logs the current time, git commit hash, and script path. Defaults to False.
    """
    if sys_info:
        # Get the current time
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        # Get the git commit hash
        git_commit = (
            subprocess.check_output(["git", "rev-parse", "HEAD"])
            .strip()
            .decode("utf-8")
        )
        # Get the path of the script
        script_path = os.path.abspath(__file__)
        # Log the current time, git commit, and script path
        logging.info(f"Script run at: {current_time}")
        logging.info(f"Git commit: {git_commit}")
        logging.info(f"Script path: {script_path}")
    process = subprocess.run(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    if process.stdout:
        logging.info(process.stdout)
    if process.returncode != 0:
        logging.error(
            f"Command '{command}' failed with return code {process.returncode}"
        )
        raise subprocess.CalledProcessError(process.returncode, command)
