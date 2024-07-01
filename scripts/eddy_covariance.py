#%% imports
import xarray as xr
import os 
import tarfile

#%% PLUMBER 2 data
# download [here](http://doi.org/10.25914/5fdb0902607e1) via NCI THREDDS Data Server
# Related paper is [Ukkola et al. (2022)](https://doi.org/10.5194/essd-14-449-2022)
#Site: BE-Bra
flux_url = r"https://thredds.nci.org.au/thredds/dodsC/ks32/CLEX_Data/PLUMBER2/v1-0/Flux/BE-Bra_2004-2014_FLUXNET2015_Flux.nc"
met_url = r"https://thredds.nci.org.au/thredds/dodsC/ks32/CLEX_Data/PLUMBER2/v1-0/Met/BE-Bra_2004-2014_FLUXNET2015_Met.nc"
urls = [flux_url, met_url]

datadir = os.path.join("..","data")
data_ec_dir = os.path.join(datadir, "exp_raw","eddy_covariance")
if not os.path.exists(data_ec_dir):
    os.makedirs(data_ec_dir)
#write data
for url in urls:
    (xr.open_dataset(url)).to_netcdf(
        os.path.join(data_ec_dir,url.split("/")[-1])
    )

#%% FluxDataKit data
# Available on Zenodo at https://doi.org/10.5281/zenodo.11370417
# Only download FLUXDATAKIT_LSM.tar.gz
os.system("zenodo_get 10.5281/zenodo.11370417 -g FLUXDATAKIT_LSM.tar.gz")

# Extract to folder
sites = ["ES-LM1","ES-LM2","ES-LMa"]
with tarfile.open("FLUXDATAKIT_LSM.tar.gz","r") as tar:
    for member in tar.getmembers():
        if any(pattern in member.name for pattern in sites):
            tar.extract(member, path=data_ec_dir)

# Remove the .tar.gz and download check file 
os.remove("FLUXDATAKIT_LSM.tar.gz")
os.remove("md5sums.txt")

