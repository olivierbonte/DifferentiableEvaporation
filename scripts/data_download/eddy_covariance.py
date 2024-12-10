# %% imports
import os
import shutil
import subprocess
import tarfile

import xarray as xr
from conf import ec_dir, sites

ec_dir.mkdir(exist_ok=True)

# %% PLUMBER 2 data
# download [here](http://doi.org/10.25914/5fdb0902607e1) via NCI THREDDS Data Server
# Related paper is [Ukkola et al. (2022)](https://doi.org/10.5194/essd-14-449-2022)
# Site: BE-Bra
flux_url = r"https://thredds.nci.org.au/thredds/dodsC/ks32/CLEX_Data/PLUMBER2/v1-0/Flux/BE-Bra_2004-2014_FLUXNET2015_Flux.nc"
met_url = r"https://thredds.nci.org.au/thredds/dodsC/ks32/CLEX_Data/PLUMBER2/v1-0/Met/BE-Bra_2004-2014_FLUXNET2015_Met.nc"
urls = [flux_url, met_url]


# write data
for url in urls:
    (xr.open_dataset(url)).to_netcdf(ec_dir / url.split("/")[-1])


# %% FluxDataKit data
# Available on Zenodo at https://doi.org/10.5281/zenodo.12818273
def check_download_move(file_name, zenodo_doi):
    if not (ec_dir / file_name).exists():
        subprocess.run(["zenodo_get", zenodo_doi, "-g", file_name])
        shutil.move(file_name, ec_dir / file_name)
    else:
        print(f"No new download, a version of {file_name} was present on disk")


zenodo_repo = "10.5281/zenodo.12818273"
check_download_move("FLUXDATAKIT_LSM.tar.gz", zenodo_repo)  # Data
check_download_move("fdk_site_info.csv", zenodo_repo)  # Tabular metadata
check_download_move(
    "fdk_site_fullyearsequence.csv", zenodo_repo
)  # Good quality sequences

# Extract to folder
with tarfile.open(ec_dir / "FLUXDATAKIT_LSM.tar.gz", "r") as tar:
    for member in tar.getmembers():
        if any(pattern in member.name for pattern in sites):
            tar.extract(member, path=ec_dir)


# %%
