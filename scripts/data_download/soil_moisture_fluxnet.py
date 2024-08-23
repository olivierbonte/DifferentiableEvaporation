# %% Imports
import zipfile
from pathlib import Path

from conf import fluxnet_dir, sites
from icoscp.station import station
from icoscp_core.icos import auth, data, meta

auth.init_config_file()
fluxnet_dir.mkdir(exist_ok=True)

# %% Download FLUXNET data for each site
for site in sites:
    print(f"site in progress: {site}")
    station_tmp = station.get(site)
    data_type = "http://meta.icos-cp.eu/resources/cpmeta/miscFluxnetArchiveProduct"
    data_list_site = meta.list_data_objects(
        datatype=data_type, station=station_tmp.info()["uri"]
    )
    if len(data_list_site) > 1:
        raise Warning(
            "More than 1 dataset for the combination of dataset and site, which is unwanted"
        )
    data.save_to_folder(data_list_site[0], str(fluxnet_dir))

# %% Unzip the FLUXNET data
for zip_file in fluxnet_dir.iterdir():
    # Check if a .zip file
    if zip_file.as_posix()[-4:] == ".zip":
        with zipfile.ZipFile(zip_file, "r") as zip_ref:
            extractdir = Path(zip_file.as_posix()[:-4])
            extractdir.mkdir(exist_ok=True)
            zip_ref.extractall(extractdir)
