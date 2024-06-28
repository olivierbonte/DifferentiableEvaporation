#%% Imports
import subprocess
subprocess.run("earthengine authenticate --quiet", shell=True, capture_output=True, text=True)
import xarray as xr
import rioxarray
import ee  
import os
ee.Initialize(opt_url = "https://earthengine-highvolume.googleapis.com")

#%% Extract to match with data from soilgrids
sites = ["BE-Bra","ES-LM1"]
datadir = os.path.join("..","data")
datarawdir = os.path.join(datadir,"exp_raw")

for site in sites:
    soilgrid_dir = os.path.join(datarawdir,"soilgrids",site)
    tif_files = [k for k in os.listdir(soilgrid_dir) if ".tif" in k]
    ds_soilgrids_site = rioxarray.open_rasterio(
        os.path.join(soilgrid_dir, tif_files[0]), mask_and_scale=True
    )
    ds_hydrosoil = xr.open_dataset(
        "ee://projects/sat-io/open-datasets/HiHydroSoilv2_0/wcres",
        engine = "ee",
        geometry = tuple(ds_soilgrids_site.rio.bounds()),
        projection = ee.Projection(
            ds_soilgrids_site.rio.crs.to_wkt(),
            #transform = ds_soilgrids_site.rio.transform()[:6]
        )
    )
##IDEA FOR WIP:
#1.) Get bbox from soilgrids in EPSG:4326
#2.) get hihydorsoil data in EPSG:4326 at this resolution
#3.) Convert back to resolution of 
#LOOK FOR CONSISTENCY WITH LSM PRODUCT!

#%%

