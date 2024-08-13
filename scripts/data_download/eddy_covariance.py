# %% imports
import xarray as xr
import os
import tarfile
import shutil
from conf import ec_dir, sites

# %% PLUMBER 2 data
# download [here](http://doi.org/10.25914/5fdb0902607e1) via NCI THREDDS Data Server
# Related paper is [Ukkola et al. (2022)](https://doi.org/10.5194/essd-14-449-2022)
# Site: BE-Bra
flux_url = r"https://thredds.nci.org.au/thredds/dodsC/ks32/CLEX_Data/PLUMBER2/v1-0/Flux/BE-Bra_2004-2014_FLUXNET2015_Flux.nc"
met_url = r"https://thredds.nci.org.au/thredds/dodsC/ks32/CLEX_Data/PLUMBER2/v1-0/Met/BE-Bra_2004-2014_FLUXNET2015_Met.nc"
urls = [flux_url, met_url]

if not os.path.exists(ec_dir):
    os.makedirs(ec_dir)
# write data
for url in urls:
    (xr.open_dataset(url)).to_netcdf(os.path.join(ec_dir, url.split("/")[-1]))

# %% FluxDataKit data
# Available on Zenodo at https://doi.org/10.5281/zenodo.12818273
# Download FLUXDATAKIT_LSM.tar.gz
os.system("zenodo_get 10.5281/zenodo.12818273 -g FLUXDATAKIT_LSM.tar.gz")
# Download metadata in tabular form
os.system("zenodo_get 10.5281/zenodo.12818273 -g fdk_site_info.csv")
# Download data containing sequences of good-quality data
os.system("zenodo_get 10.5281/zenodo.12818273 -g fdk_site_fullyearsequence.csv")

# Extract to folder
with tarfile.open("FLUXDATAKIT_LSM.tar.gz", "r") as tar:
    for member in tar.getmembers():
        if any(pattern in member.name for pattern in sites):
            tar.extract(member, path=ec_dir)

# Remove the .tar.gz and download check file
os.remove("FLUXDATAKIT_LSM.tar.gz")
os.remove("md5sums.txt")

# Move metadata to correct folder
shutil.move("fdk_site_info.csv", os.path.join(ec_dir, "fdk_site_info.csv"))
shutil.move(
    "fdk_site_fullyearsequence.csv",
    os.path.join(ec_dir, "fdk_site_fullyearsequence.csv"),
)
