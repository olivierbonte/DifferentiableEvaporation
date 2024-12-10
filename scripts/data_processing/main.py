# %% Imports
import logging
from datetime import datetime
from importlib.machinery import SourceFileLoader
from pathlib import Path

from conf import conf_module

project_folder = Path(__file__).parent.parent.parent.resolve()
py_functions_module = SourceFileLoader(
    "py_functions", (project_folder / "src" / "py_functions.py").as_posix()
).load_module()

# %% Logging set up
current_datetime = datetime.now().strftime("%Y-%m-%d_%H-%M")
log_file = Path(conf_module.dataprodir / f"data_process_log_{current_datetime}.txt")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
    handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
)

# %% Logging of the subprocesses
logging.info("\n Starting data processing scripts \n")
logging.info("\n 1: Processing eddy covaraince data \n")
py_functions_module.run_and_log("python eddy_covariance.py", sys_info=True)
logging.info("\n 2: Processing soil data \n")
py_functions_module.run_and_log("python soil.py")
logging.info("\n 3: Processing vegetation data \n")
py_functions_module.run_and_log("python vegetation.py")
logging.info("\n 4: Processing soil moisture data from flux towers \n")
py_functions_module.run_and_log("python soil_moisture_fluxnet.py")
