# %%
import rioxarray
import os
import glob
import xarray as xr
from rasterio.enums import Resampling

datadir = os.path.join("..", "..", "data")
soilgrids_dir = os.path.join(datadir, "exp_raw", "soilgrids")
hihydrosoil_dir = os.path.join(datadir, "exp_raw", "hihydrosoil")
soil_pro_dir = os.path.join(datadir, "exp_pro", "soil")
if not os.path.exists(soil_pro_dir):
    os.makedirs(soil_pro_dir)
sites = os.listdir(soilgrids_dir)

# %% Datacube creation
for site in sites:
    print(f"site in progress: {site}")
    ds_soilgrids = xr.open_dataset(
        glob.glob(os.path.join(soilgrids_dir, site, "*.nc"))[0], decode_coords="all"
    )
    ds_hihydrosoil = xr.open_dataset(
        os.path.join(hihydrosoil_dir, site + ".nc"), decode_coords="all"
    )

    ## Convert to 1 cube: (depth, lon, lat) + crs:ESRI_54052
    # Order needed for reprojection
    ds_hihydrosoil = ds_hihydrosoil.transpose("depth", "lat", "lon")
    # To allow reprojection with rasterio: x,y coordinate names
    ds_hihydrosoil = ds_hihydrosoil.rename({"lon": "x", "lat": "y"})
    ds_hihydrosoil = ds_hihydrosoil.rio.reproject_match(
        ds_soilgrids, resampling=Resampling.bilinear
    )
    # Add depth dimension to soilgrids data + keep attributes to reassign
    da_list = []
    for var in ds_soilgrids.data_vars:
        print(var.split("_"))
        attr_var_dict = ds_soilgrids[var].attrs
        da_tmp = ds_soilgrids[var] / 10  # g/kg -> %, also needed for uncertainty layer
        da_tmp.attrs = attr_var_dict
        if "units" in da_tmp.attrs:
            da_tmp.attrs["units"] = "%"  # g/kg -> %
        else:  # uncertainty layer
            da_tmp.attrs["description"] = attr_var_dict["description"][:-3]
        da_tmp.name = var.split("_")[2] + "_" + var.split("_")[0] + "_" + "percentage"
        da_tmp = da_tmp.expand_dims({"depth": [var.split("_")[1]]})
        da_list.append(da_tmp)
    ds_soilgrids_depth = xr.merge(da_list)
    # Merge soilgrids and hihydrosoil
    ds_soil = xr.merge([ds_soilgrids_depth, ds_hihydrosoil])
    # Add useful metadata
    [ds_soilgrids.attrs.pop(attr) for attr in ["url_doi", "url_OGC", "units"]]
    ds_soil.attrs = ds_soilgrids.attrs
    ds_soil.to_netcdf(os.path.join(soil_pro_dir, site + "_soil_cube.nc"))

    # Spatial data aggregation
    ds_soil_agg = ds_soil.mean(dim=["x", "y"], keep_attrs=True)
    ds_soil_agg.to_netcdf(os.path.join(soil_pro_dir, site + "_soil_horizontal_agg.nc"))

    # Aggregate with a weighted mean (respective to its depth)
    # -> Do this in function of the depth of layer 2 we use in the model itself!
