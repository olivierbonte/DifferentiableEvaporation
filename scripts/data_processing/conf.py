# Objective: allows import from the conf.py as defined in the data_download folder
import importlib.util
from pathlib import Path

_conf_path = Path(__file__).resolve().parent.parent / "data_download" / "conf.py"
_spec = importlib.util.spec_from_file_location("data_download_conf", _conf_path)
conf_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(conf_module)

# Re-export the logger setup so scripts only need to import from `conf`
logger = conf_module.logger
LOG_FORMAT = conf_module.LOG_FORMAT
