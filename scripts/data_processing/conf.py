# Objective: allows import from the conf.py as defined in the data_download folder
import os
from importlib.machinery import SourceFileLoader
from pathlib import Path

scripts_folder = Path(__file__).resolve().parent.parent
conf_module = SourceFileLoader(
    "conf", os.path.join(scripts_folder, "data_download", "conf.py")
).load_module()
