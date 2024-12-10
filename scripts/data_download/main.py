import logging
import os
import subprocess
from datetime import datetime
from pathlib import Path

from conf import datarawdir

## Logging set up
log_file = Path(datarawdir / "data_download_log.txt")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
    handlers=[logging.FileHandler(log_file, mode="w"), logging.StreamHandler()],
)


## Function to run a subprocess and log its output
def run_and_log(command, sys_info=False):
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


## Logging of the subprocesses
logging.info("\n Starting data download script \n")
logging.info("\n 1: Downloading eddy covaraince data \n")
run_and_log("python eddy_covariance.py", sys_info=True)
logging.info("\n 2: Downloading data from SoilGrids \n")
run_and_log("python soil_grids.py")
logging.info("\n 3: Download data from HiHydroSoil \n")
run_and_log("python hihydrosoil.py")
logging.info("\n 4: Download data on vegetation cover fractions \n")
run_and_log("python vegetation_VCF.py")
