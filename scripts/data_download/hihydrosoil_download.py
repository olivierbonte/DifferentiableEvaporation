# %% Imports
import glob

import ee
import rioxarray  # noqa: F401  (registers the .rio accessor)
import shapely
import xarray as xr
from conf import gee_project_id, hihydrosoil_dir, logger, sites, soilgrids_dir
from xee import helpers


def main():
    try:
        ee.Initialize(
            project=gee_project_id,
            opt_url="https://earthengine-highvolume.googleapis.com",
        )
    except ee.EEException:
        ee.Authenticate(auth_mode="notebook")
        ee.Initialize(
            project=gee_project_id,
            opt_url="https://earthengine-highvolume.googleapis.com",
        )

    hihydrosoil_dir.mkdir(exist_ok=True)

    # %% Define variables of interest from HiHydroSoilv2
    # https://gee-community-catalog.org/projects/hihydro_soil/
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
    layer_depths = [0.05, 0.1, 0.15, 0.3, 0.4, 1.0]  # m
    base_url = "projects/sat-io/open-datasets/HiHydroSoilv2_0/"

    # %% Extract at the native GEE HiHydroSoil resolution (EPSG:4326) covering the
    # soilgrids data; reprojecting to the soilgrids grid happens in data processing
    for site in sites:
        logger.info(f"site in progress: {site}")
        soilgrids_file = glob.glob(str(soilgrids_dir / site / ("*" + site + "*.nc")))
        ds_soilgrids = xr.open_dataset(soilgrids_file[0], decode_coords="all")
        soilgrids_bbox = shapely.box(*ds_soilgrids.rio.bounds())

        da_list = []
        for gee_var, var in var_dict.items():
            collection = ee.ImageCollection(base_url + gee_var)
            native_grid_params = helpers.extract_grid_params(collection)
            native_scale = native_grid_params["crs_transform"][0]
            # Buffer of 1 pixel for the bilinear interpolation at the edges in procesing
            grid_params = helpers.fit_geometry(
                soilgrids_bbox,
                geometry_crs=ds_soilgrids.rio.crs.to_wkt(),  # CRS of soilgrids_bbox
                buffer=native_scale,  # Buffer of 1 pixel (in units of geometry_crs)
                grid_crs=native_grid_params["crs"],  # CRS of target (GEE grid)
                grid_scale=(native_scale, -native_scale),  # grid_crs units
            )
            ds_temp = xr.open_dataset(collection, engine="ee", **grid_params)
            # Convert to correct units by multiplying by 0.0001, see:
            # https://gee-community-catalog.org/projects/hihydro_soil/
            ds_temp = ds_temp * 0.0001
            ds_temp = ds_temp.rename({list(ds_temp.data_vars.keys())[0]: var})
            # write crs with rioxarray
            ds_temp = ds_temp.rio.write_crs(grid_params["crs"])
            ds_temp[var].attrs["GEE_var_name"] = gee_var
            ds_temp[var].attrs["units"] = units_dict[var]
            ds_temp[var].attrs["full_name"] = full_names_dict[var]
            ds_temp[var].attrs["url"] = (
                "https://gee-community-catalog.org/projects/hihydro_soil/"
            )
            da_list.append(ds_temp[var])
        # Data cube creation
        ds_hihydrosoil = xr.merge(da_list, compat="no_conflicts")
        # Rename depth dimension
        ds_hihydrosoil = ds_hihydrosoil.rename({"time": "depth"})
        ds_hihydrosoil.depth.attrs["Description"] = (
            "Intervals of soil depth below the surface"
        )
        ds_hihydrosoil = ds_hihydrosoil.assign_coords(depth=depth_values)
        # Add layer depth variable
        da_depth = xr.DataArray(
            layer_depths, dims=("depth"), coords={"depth": depth_values}
        )
        da_depth.attrs = {"units": "m"}
        ds_hihydrosoil["layer_depth"] = da_depth
        # Clean up metadata
        [
            ds_hihydrosoil.attrs.pop(attr)
            for attr in ["GEE_var_name", "units", "full_name"]
        ]
        ds_hihydrosoil.attrs["url"] = (
            "https://gee-community-catalog.org/projects/hihydro_soil/"
        )
        ds_hihydrosoil.attrs["url_report"] = (
            "https://www.futurewater.nl/wp-content/uploads/2020/10/HiHydroSoil-v2.0-"
            "High-Resolution-Soil-Maps-of-Global-Hydraulic-Properties_v2.pdf"
        )
        # Write to disk as NetCDF
        ds_hihydrosoil.to_netcdf(hihydrosoil_dir / (site + ".nc"))


if __name__ == "__main__":
    main()
