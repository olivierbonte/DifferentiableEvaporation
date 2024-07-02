# %% Imports
import subprocess

subprocess.run(
    "earthengine authenticate --quiet", shell=True, capture_output=True, text=True
)
import xarray as xr
import rioxarray
import ee
import os
import glob
from rasterio.enums import Resampling
from conf import soilgrids_dir, hihydrosoil_dir, sites

ee.Initialize(opt_url="https://earthengine-highvolume.googleapis.com")

# %% Define variables of interest from HiHydroSoilv2
# https://gee-community-catalog.org/projects/hihydro_soil/
var_list_gee = ["N", "wcsat", "wcres", "wcpf2", "wcpf3", "wcpf4-2"]
var_dict = {
    "N": "n",
    "wcsat": "w_sat",
    "wcres": "w_res",
    "wcpf2": "w_fc",
    "wcpf3": "w_crit",
    "wcpf4-2": "w_wp",
}  # key = GEE, value = own_name
full_names_dict = {
    "n": "N parameter for Mualem Van Genuchten Equation",
    "w_sat": "Saturated soil moisture",
    "w_res": "Residual soil moisture",
    "w_fc": "Soil moisture at field capacity (pF2)",
    "w_crit": "Soil moisture at critical point (pF3)",
    "w_wp": "Soil moisture at permanent wilting point (pF4.2)",
}
units_dict = {
    "n": "[-]",
    "w_sat": "[m^3/m^3]",
    "w_res": "[m^3/m^3]",
    "w_fc": "[m^3/m^3]",
    "w_crit": "[m^3/m^3]",
    "w_wp": "[m^3/m^3]",
}
depth_values = ["0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm"]
base_url = "ee://projects/sat-io/open-datasets/HiHydroSoilv2_0/"

# %% Extract to match with data from soilgrids
if not os.path.exists(hihydrosoil_dir):
    os.makedirs(hihydrosoil_dir)

for site in sites:
    # Get bounding box from soilgrids data
    soilgrids_file = glob.glob(os.path.join(soilgrids_dir, site, "*" + site + "*.nc"))
    ds_soilgrids = xr.open_dataset(soilgrids_file[0], decode_coords="all")
    ds_soilgrids_wgs84 = ds_soilgrids.rio.reproject(
        dst_crs="EPSG:4326", resampling=Resampling.bilinear
    )
    bbox_wgs84 = ds_soilgrids_wgs84.rio.bounds()
    scale = ds_soilgrids_wgs84.rio.resolution()  # spatial spacing in °

    da_list = []
    for gee_var, var in var_dict.items():
        # See https://github.com/google/Xee?tab=readme-ov-file#how-to-use
        ds_temp = xr.open_dataset(
            base_url + gee_var,
            engine="ee",
            geometry=tuple(bbox_wgs84),
            projection=ee.Projection(
                crs=str(ds_soilgrids_wgs84.rio.crs),
                transform=ds_soilgrids_wgs84.rio.transform()[:6],
            ),
            scale=scale,
        )
        # Convert to correct units by multiplying by 0.0001, see:
        # https://gee-community-catalog.org/projects/hihydro_soil/
        ds_temp = ds_temp * 0.0001
        ds_temp = ds_temp.rename({list(ds_temp.data_vars.keys())[0]: var})
        ds_temp[var].attrs["GEE_var_name"] = gee_var
        ds_temp[var].attrs["units"] = units_dict[var]
        ds_temp[var].attrs["full_name"] = full_names_dict[var]
        ds_temp[var].attrs[
            "url"
        ] = "https://gee-community-catalog.org/projects/hihydro_soil/"
        da_list.append(ds_temp[var])
    # Data cube creation
    ds_hihydrosoil = xr.merge(da_list)
    # Rename depth dimension
    ds_hihydrosoil = ds_hihydrosoil.rename({"time": "depth"})
    ds_hihydrosoil.depth.attrs["Description"] = (
        "Intervals of soil depth below the surface"
    )
    ds_hihydrosoil = ds_hihydrosoil.assign_coords(depth=depth_values)
    # Clean up metadata
    [ds_hihydrosoil.attrs.pop(attr) for attr in ["GEE_var_name", "units", "full_name"]]
    ds_hihydrosoil.attrs["url"] = (
        "https://gee-community-catalog.org/projects/hihydro_soil/"
    )
    ds_hihydrosoil.attrs["url_report"] = (
        "https://www.futurewater.nl/wp-content/uploads/2020/10/HiHydroSoil-v2.0-"
        "High-Resolution-Soil-Maps-of-Global-Hydraulic-Properties_v2.pdf"
    )
    # Write to disk as NetCDF
    ds_hihydrosoil.to_netcdf(os.path.join(hihydrosoil_dir, site + ".nc"))
