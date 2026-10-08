# %% Imports
import zipfile
from pathlib import Path

import numpy as np
from conf import fluxnet_dir, sites
from icoscp_core.icos import auth, data, meta

auth.init_config_file()  # Only do this once per computer
fluxnet_dir.mkdir(exist_ok=True)

# %% Download FLUXNET data for each site

# Get all possible stations
all_known_stations = np.array(meta.list_stations())
id_list = np.array([station_temp.id for station_temp in all_known_stations])

# Get data for selected stations
for site in sites:
    print(f"site in progress: {site}")
    station_tmp_info = all_known_stations[id_list == site][0]
    # All possible datatypes via meta.list_datatypes()
    data_types = [
        # "http://meta.icos-cp.eu/resources/cpmeta/miscFluxnetArchiveProduct",
        "http://meta.icos-cp.eu/resources/cpmeta/etcArchiveProduct",  # MetaData + more recent data
        "http://meta.icos-cp.eu/resources/cpmeta/miscFluxnetProduct",  # Fluxnet Half-hourly product
    ]
    data_list_site = meta.list_data_objects(
        datatype=data_types, station=station_tmp_info.uri
    )
    if len(data_list_site) != len(data_types):
        raise Warning(
            f"""Not the same amount datasets ({len(data_list_site)}) as 
            dataset types ({len(data_types)}), which is unwanted"""
        )
    for data_item in data_list_site:
        data.save_to_folder(data_item, str(fluxnet_dir))

# %% Unzip the FLUXNET data
for zip_file in fluxnet_dir.iterdir():
    # Check if a .zip file
    if zip_file.as_posix()[-4:] == ".zip":
        with zipfile.ZipFile(zip_file, "r") as zip_ref:
            extractdir = Path(zip_file.as_posix()[:-4])
            extractdir.mkdir(exist_ok=True)
            zip_ref.extractall(extractdir)
        # Delete zip file
        zip_file.unlink()
# %%
