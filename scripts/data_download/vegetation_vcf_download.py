# %% Imports
import glob

import ee
import shapely
import xarray as xr
from conf import ec_dir, gee_project_id, logger, sites, veg_dir
from xee import helpers


def main():
    veg_dir.mkdir(exist_ok=True)
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

    # %% Define variables of interest
    var_list_vcf = [
        "Percent_Tree_Cover",
        "Percent_NonTree_Vegetation",
        "Percent_NonVegetated",
    ]
    collection_vcf = "MODIS/061/MOD44B"
    nr_pixels = 2  # extra pixels to each side of the center pixel
    # MODIS sinusoidal projection (SR-ORG:6974 in GEE, unknown to pyproj)
    crs_modis_sinusoidal = (
        "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
    )
    # source: https://gis.stackexchange.com/questions/272639/using-modis-sinusoidal-projection

    # %% Extract data for 5 x 5 pixels (~1160m x 1160m) on the native MODIS grid
    for site in sites:
        logger.info(f"site in progress: {site}")
        ec_file = glob.glob(str(ec_dir / ("*" + site + "*Flux.nc")))
        ds_ec = xr.open_dataset(ec_file[0], decode_coords="all")
        lat, lon = ds_ec["latitude"].item(), ds_ec["longitude"].item()
        start_date = ds_ec.time[0].values.astype("datetime64[D]").astype(str)
        end_date = ds_ec.time[-1].values.astype("datetime64[D]").astype(str)
        collection = (
            ee.ImageCollection(collection_vcf)
            .filterDate(start_date, end_date)
            .select(var_list_vcf)
        )
        native_grid_params = helpers.extract_grid_params(collection)
        resolution = native_grid_params["crs_transform"][0]
        # Pixel containing the site + nr_pixels to each side (by buffering the
        # site with nr_pixels * resolution), analogous to the soilgrids data.
        # The MODIS grid origin is a whole number of pixels from 0, so the grid of
        # fit_geometry coincides with the native grid (no resampling)
        grid_params = helpers.fit_geometry(
            shapely.Point(lon, lat),
            geometry_crs="EPSG:4326",  # CRS of lon/lat point
            buffer=nr_pixels * resolution,
            grid_crs=crs_modis_sinusoidal,  # CRS of target that pyproj can handle (GEE grid)
            grid_scale=(resolution, -resolution),
        )
        grid_params["crs"] = native_grid_params["crs"]
        da_modis_vcf = (
            xr.open_dataset(collection, engine="ee", **grid_params)
            .to_array("band")
            .transpose("time", "band", "y", "x")
        )
        # rename array to avoid writing issues
        da_modis_vcf.name = "MODIS_VCF"
        da_modis_vcf.attrs = dict(
            collection=collection_vcf,
            crs=grid_params["crs"],
            crs_transform=grid_params["crs_transform"],
            resolution=resolution,
            edge_size=grid_params["shape_2d"][0],
            central_lat=lat,
            central_lon=lon,
            time_coverage_start=start_date,
            time_coverage_end=end_date,
        )
        # take the spatial mean over the cube
        da_modis_vcf_mean = da_modis_vcf.mean(dim=["x", "y"], keep_attrs=True)
        # drop incorrect attributes after mean operation
        [
            da_modis_vcf_mean.attrs.pop(attr)
            for attr in ["crs_transform", "resolution", "edge_size"]
        ]
        # write to disk
        da_modis_vcf.to_netcdf(veg_dir / (site + "_cube.nc"))
        da_modis_vcf_mean.to_netcdf(veg_dir / (site + "_horizontal_agg.nc"))


if __name__ == "__main__":
    main()
