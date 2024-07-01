# %%
import rioxarray
import os
import glob
import xarray as xr

datadir = os.path.join("..", "..", "data")
soilgrids_dir = os.path.join(datadir, "exp_raw", "soilgrids")
hihydrosoil_dir = os.path.join(datadir, "exp_raw", "hihydrosoil")
sites = os.listdir(soilgrids_dir)

# %%
for site in sites:
    ds_soilgrids = xr.open_dataset(
        glob.glob(os.path.join(soilgrids_dir, site, "*.nc"))[0], decode_coords="all"
    )
    ds_hihydrosoil = xr.open_dataset(
        os.path.join(hihydrosoil_dir, site + ".nc"), decode_coords="all"
    )

    # Convert to 1 cube: (depth, lon, lat) + crs:ESRI_54052
    ds_hihydrosoil = ds_hihydrosoil.rio.set_spatial_dims("lon", "lat")
    ds_hihydrosoil = ds_hihydrosoil.transpose("depth", "lat", "lon")
    ds_hihydrosoil = ds_hihydrosoil.rio.write_crs(ds_hihydrosoil.attrs["crs"])
    ds_hihydrosoil = ds_hihydrosoil.rio.reproject(
        ds_soilgrids.spatial_ref.attrs["crs_wkt"]
    )


# %%
