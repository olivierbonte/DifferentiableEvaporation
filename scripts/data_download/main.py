import logging
import os
import subprocess
from datetime import datetime
from pathlib import Path

from conf import datarawdir

## Info for logging
# Get the current time
current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
# Get the git commit hash
git_commit = (
    subprocess.check_output(["git", "rev-parse", "HEAD"]).strip().decode("utf-8")
)
# Get the path of the script
script_path = os.path.abspath(__file__)

## Logging set up
log_file = Path(datarawdir / "data_download_log.txt")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
    handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
)

# Log the current time, git commit, and script path
logging.info(f"Script run at: {current_time}")
logging.info(f"Git commit: {git_commit}")
logging.info(f"Script path: {script_path}")
# Logging of the subprocesses
logging.info("\n Starting data download script \n")
logging.info("\n 1: Downloading eddy covaraince data \n")
subprocess.run("python eddy_covariance.py", shell=True)
logging.info("\n 2: Downloading data from SoilGrids \n")
subprocess.run("python soil_grids.py", shell=True)
logging.info("\n 3: Download data from HiHydroSoil \n")
subprocess.run("python hihydrosoil.py", shell=True)
logging.info("\n 4: Download data on vegetation cover fractions \n")
subprocess.run("python vegetation_VCF.py", shell=True)
