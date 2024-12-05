# %%
import glob

import rioxarray
import xarray as xr
from conf import conf_module
from rasterio.enums import Resampling

conf_module.soil_pro_dir.mkdir(exist_ok=True)

# %% Datacube creation
for site in conf_module.sites:
    print(f"site in progress: {site}")
    ds_soilgrids = xr.open_mfdataset(
        glob.glob(str(conf_module.soilgrids_dir / site / "*.nc")),
        decode_coords="all",
    )
    ds_hihydrosoil = xr.open_dataset(
        conf_module.hihydrosoil_dir / (site + ".nc"), decode_coords="all"
    )
    da_root_depth = xr.open_dataarray(
        conf_module.essd_dir / "NetCDF" / (site + "_root_depth.nc"), decode_coords="all"
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
    # Assert correct order in depth (as in hihydrosoil)
    ds_soilgrids_depth = ds_soilgrids_depth.reindex(depth=ds_hihydrosoil.depth)
    # Merge soilgrids and hihydrosoil
    ds_soil = xr.merge([ds_soilgrids_depth, ds_hihydrosoil])
    # Add root_depth
    ds_soil["root_depth"] = float(da_root_depth.values)
    ds_soil["root_depth"].attrs = da_root_depth.attrs
    # Add useful metadata
    [ds_soilgrids.attrs.pop(attr) for attr in ["url_doi", "url_OGC"]]
    ds_soil.attrs = ds_soilgrids.attrs
    ds_soil.to_netcdf(conf_module.soil_pro_dir / (site + "_soil_cube.nc"))

    # Spatially horizontal data aggregation
    ds_soil_agg = ds_soil.mean(dim=["x", "y"], keep_attrs=True)
    ds_soil_agg.to_netcdf(conf_module.soil_pro_dir / (site + "_soil_horizontal_agg.nc"))

    # Add vertical aggregation: weighted (on layer depth) mean per layer,
    # max depth = root_depth
    # ds_soil_agg_vert = ds_soil_agg.mean(dim=["depth"], keep_attrs=True)
    # Select layers which are (partly at least) in the available root depth
    bool_vert = (
        ds_soil_agg.layer_depth.cumsum() - ds_soil_agg.layer_depth
    ) < ds_soil_agg.root_depth.values / 100
    ds_soil_agg_sel = ds_soil_agg.isel(depth=bool_vert)
    weights = ds_soil_agg_sel.layer_depth
    # For the deepest layer, the weight is root depth - bottom depth of previous layer!
    weights[-1] = (
        ds_soil_agg.root_depth / 100 - ds_soil_agg_sel.layer_depth.cumsum()[-2]
    )
    # Usage of weights explained here: https://tutorial.xarray.dev/fundamentals/03.4_weighted.html
    ds_soil_agg_sel_weighted = ds_soil_agg_sel.weighted(weights)
    ds_soil_agg_vert = ds_soil_agg_sel_weighted.mean(dim="depth", keep_attrs=True)
    ds_soil_agg_vert.to_netcdf(conf_module.soil_pro_dir / (site + "_total_agg.nc"))

# %%
