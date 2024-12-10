# %% Imports
import glob

import cubo
import ee
import xarray as xr
from conf import ec_dir, sites, veg_dir

veg_dir.mkdir(exist_ok=True)
ee.Initialize(opt_url="https://earthengine-highvolume.googleapis.com")

# %% Define variables of interest
var_list_vcf = [
    "Percent_Tree_Cover",
    "Percent_NonTree_Vegetation",
    "Percent_NonVegetated",
]
collection_vcf = "MODIS/006/MOD44B"


# %% Extract data for 1250m x 1250m grid
for site in sites:
    print(f"site in progress: {site}")
    ec_file = glob.glob(str(ec_dir / ("*" + site + "*Flux.nc")))
    ds_ec = xr.open_dataset(ec_file[0], decode_coords="all")
    # Get data in UTM grid (see https://cubo.readthedocs.io/en/stable/index.html#how-does-it-work)
    da_modis_vcf = cubo.create(
        lat=ds_ec["latitude"].item(),
        lon=ds_ec["longitude"].item(),
        collection=collection_vcf,
        bands=var_list_vcf,
        start_date=ds_ec.time[0].values.astype("datetime64[D]").astype(str),
        end_date=ds_ec.time[-1].values.astype("datetime64[D]").astype(str),
        edge_size=5,  # nr of pixels, analogous to the 1250m x 1250m of the soilgrids data
        resolution=250,  # resolution in m (matching original)
        gee=True,
    )
    # rename array to avoid writing issues
    da_modis_vcf.name = "MODIS_VCF"
    # take the spatial mean over the
    da_modis_vcf_mean = da_modis_vcf.mean(dim=["x", "y"], keep_attrs=True)
    # drop incorrect attributes after mean operation
    [da_modis_vcf_mean.attrs.pop(attr) for attr in ["resolution", "edge_size"]]
    # write to disk
    da_modis_vcf.to_netcdf(veg_dir / (site + "_cube.nc"))
    da_modis_vcf_mean.to_netcdf(veg_dir / (site + "_horizontal_agg.nc"))

# Note: CUBO performs the nearest neighbour resampling be default, see
# https://developers.google.com/earth-engine/guides/resample  and
# https://developers.google.com/earth-engine/guides/resample

# %%
