# %% Imports
import logging
from datetime import datetime
from importlib.machinery import SourceFileLoader
from pathlib import Path

from conf import datarawdir
from icoscp_core.icos import auth

project_folder = Path(__file__).parent.parent.parent.resolve()
py_functions_module = SourceFileLoader(
    "py_functions", (project_folder / "src" / "py_functions.py").as_posix()
).load_module()
# %% Set up log in icos
auth.init_config_file()

# %% Logging set up
datarawdir.mkdir(exist_ok=True, parents=True)
current_datetime = datetime.now().strftime("%Y-%m-%d_%H-%M")
log_file = Path(datarawdir / f"data_download_log_{current_datetime}.txt")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
    handlers=[logging.FileHandler(log_file, mode="w"), logging.StreamHandler()],
)

# %%  Logging of the subprocesses
logging.info("\n Starting data download script \n")
logging.info("\n 1: Downloading eddy covaraince data \n")
py_functions_module.run_and_log("python eddy_covariance.py", sys_info=True)
logging.info("\n 2: Downloading data from SoilGrids \n")
py_functions_module.run_and_log("python soil_grids.py")
logging.info("\n 3: Download data from HiHydroSoil \n")
py_functions_module.run_and_log("python hihydrosoil.py")
logging.info("\n 4: Download data on vegetation cover fractions \n")
py_functions_module.run_and_log("python vegetation_VCF.py")
logging.info("\n 5: Downloading auxiliary data from flux towers")
py_functions_module.run_and_log("python soil_moisture_fluxnet.py")
