# %%
import rioxarray
import os
import xarray as xr
import numpy as np
from pyproj import CRS, Transformer
from soilgrids import SoilGrids
from owslib.wcs import WebCoverageService
from conf import datadir, ec_dir, sites

# %% Clay information
var_of_interest = "clay"
soil_grids = SoilGrids()
url = soil_grids.MAP_SERVICES[var_of_interest]["link"]
units = soil_grids.MAP_SERVICES[var_of_interest]["units"]
wcs = WebCoverageService(url, version="1.0.0")
mean_names = [k for k in wcs.contents.keys() if k.find("mean") != -1]
uncertainty_names = [k for k in wcs.contents.keys() if k.find("uncertainty") != -1]

# %% Store and process data for each site
for site in sites:
    files = [k for k in os.listdir(ec_dir) if site in k]
    ds_site = xr.open_dataset(os.path.join(ec_dir, files[0]))
    output_folder = os.path.join(datadir, "exp_raw", "soilgrids", site)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # EPSG:4326 coordinates
    lon = float(ds_site["longitude"].values[0][0])
    lat = float(ds_site["latitude"].values[0][0])

    # Homolsine projection is used (https://epsg.io/54052 for proj4), see
    # https://www.isric.org/explore/soilgrids/faq-soilgrids#How_can_I_use_the_Homolosine_projection
    crs_homolsine = CRS.from_proj4(
        "+proj=igh +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
    )
    transformer = Transformer.from_crs(4326, crs_homolsine, always_xy=True)
    coords_homolsine = transformer.transform(lon, lat)

    # Take bounding box of 4 km around site
    dist = 2000  # m
    ds_list = []
    for id in mean_names + uncertainty_names:
        # Download raw data
        out_path_tif = os.path.join(output_folder, id + ".tif")
        out_path_nc = os.path.join(output_folder, id + ".nc")
        data = soil_grids.get_coverage_data(
            service_id="clay",
            coverage_id=id,
            west=coords_homolsine[0] - dist,
            east=coords_homolsine[0] + dist,
            south=coords_homolsine[1] - dist,
            north=coords_homolsine[1] + dist,
            crs="urn:ogc:def:crs:EPSG::152160",  # Homolsine;
            output=out_path_tif,
        )

        da_temp = rioxarray.open_rasterio(out_path_tif, mask_and_scale=True)
        # Extra masking for uncertainty at int16 max of 32767
        if da_temp.max().values == 32767:
            da_temp = da_temp.where(da_temp != da_temp.max())
        da_temp.name = id
        da_temp.attrs = {
            "url_doi": "https://doi.org/10.5194/soil-7-217-2021",
            "url_OGC": "https://maps.isric.org/",
        }
        if "mean" in id:
            da_temp.attrs["units"] = units
        else:
            da_temp.attrs["description"] = (
                "(95th percentile - 5th percentile)/median*10"
            )
        ds_list.append(da_temp.isel(band=0))

    # Select center pixel (around point) + 2 extra pixels to each side
    # --> 1250 x 1250 m cube, close to 1500 x 1500 m MODIS LAI or
    # 1000 x 1000 Copernicus LAI from Ukkola et al. (2019),
    # https://doi.org/10.5194/essd-14-449-2022
    # Make 1 datacube per site
    nr_pixels = 2
    ds_cube = xr.merge(ds_list)
    ds_point = ds_cube.sel(
        x=coords_homolsine[0], y=coords_homolsine[1], method="nearest"
    )
    x_index = np.where((ds_point.x == ds_cube.x).values)[0][0]
    y_index = np.where((ds_point.y == ds_cube.y).values)[0][0]
    ds_cube = ds_cube.isel(
        x=range(x_index - nr_pixels, x_index + nr_pixels + 1),
        y=range(y_index - nr_pixels, y_index + nr_pixels + 1),
    )
    # Add metadata
    ds_cube.attrs["site_code"] = ds_site.attrs["site_code"]
    ds_cube.attrs["site_name"] = ds_site.attrs["site_name"]
    ds_cube.attrs["lon_EPSG_4326"] = lon
    ds_cube.attrs["lat_EPSG_4326"] = lat
    ds_cube.attrs["x_ESRI_54052"] = coords_homolsine[0]
    ds_cube.attrs["y_ESRI_54052"] = coords_homolsine[1]
    ds_cube.attrs.pop("units")
    # Write to NetCDF
    ds_cube.to_netcdf(os.path.join(output_folder, site + "_" + var_of_interest + ".nc"))
