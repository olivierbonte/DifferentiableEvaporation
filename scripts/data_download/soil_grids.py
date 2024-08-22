# %%
import os
import shutil
from glob import glob

import numpy as np
import rioxarray
import xarray as xr
from conf import datarawdir, ec_dir, essd_dir, sites, soilgrids_dir
from osgeo import gdal
from owslib.wcs import WebCoverageService
from pyproj import CRS, Transformer
from soilgrids import SoilGrids

# %% Clay information
var_of_interest = "clay"
soil_grids = SoilGrids()
url = soil_grids.MAP_SERVICES[var_of_interest]["link"]
units = soil_grids.MAP_SERVICES[var_of_interest]["units"]
wcs = WebCoverageService(url, version="1.0.0")
mean_names = [k for k in wcs.contents.keys() if k.find("mean") != -1]
uncertainty_names = [k for k in wcs.contents.keys() if k.find("uncertainty") != -1]

# %% European Soil Database Derived data: root depth
netcdf_dir = essd_dir / "NetCDF"
netcdf_dir.mkdir(exist_ok=True)
netcdf_file = netcdf_dir / "soil_depth_test.nc"
dataset = gdal.Open(essd_dir / "STU_EU_DEPTH_ROOTS.rst")
x_size = dataset.RasterXSize
y_size = dataset.RasterYSize
bands = dataset.RasterCount
geo_transform = dataset.GetGeoTransform()
projection = dataset.GetProjection()
# Create the NetCDF file
netcdf_driver = gdal.GetDriverByName("NetCDF")
netcdf_ds = netcdf_driver.Create(
    str(netcdf_file), x_size, y_size, bands, gdal.GDT_Float32
)
# Allows for multiple bands if present
for i in range(1, bands + 1):
    rst_band = dataset.GetRasterBand(i)
    data = rst_band.ReadAsArray()
    netcdf_band = netcdf_ds.GetRasterBand(i)
    netcdf_band.WriteArray(data)
# Set georeference and projection
# Manually define the CRS, which is ETRS89-extended / LAEA Europe (https://epsg.io/3035),
# which is manually read from metadata (.RDC file)
crs_epsg3035 = CRS.from_epsg(3035)
# proj4_epsg3035 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"
netcdf_ds.SetGeoTransform(geo_transform)
netcdf_ds.SetProjection(crs_epsg3035.to_wkt())
# Close the datasets
dataset = None
netcdf_ds = None

# read netcdf with xarray, append features and write out
da_root_depth = xr.open_dataarray(netcdf_file, decode_coords="all").load()
da_root_depth.close()
da_root_depth = da_root_depth.rename("root_depth")
da_root_depth.attrs = {
    "units": "cm",
    "long_name": "Depth available to roots",
    "url": "https://esdac.jrc.ec.europa.eu/content/european-soil-database-derived-data#tabs-0-description=0",
}
da_root_depth.to_netcdf(datarawdir / "temp.nc")
shutil.move(datarawdir / "temp.nc", netcdf_file)
# %% Store and process data for each site
for site in sites:
    print(f"site in progress: {site}")
    files = glob(str(ec_dir / ("*" + site + "*FLUXDATAKIT_Flux.nc")))
    ds_site = xr.open_dataset(files[0])
    output_folder = soilgrids_dir / site
    output_folder.mkdir(exist_ok=True)

    # EPSG:4326 coordinates
    lon = float(ds_site["longitude"].values[0][0])
    lat = float(ds_site["latitude"].values[0][0])

    # Homolsine projection is used (https://epsg.io/54052 for proj4), see
    # https://www.isric.org/explore/soilgrids/faq-soilgrids#How_can_I_use_the_Homolosine_projection
    crs_homolsine = CRS.from_proj4(
        "+proj=igh +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
    )
    # Transofmer lat/lon to both CRS
    transformer = Transformer.from_crs(4326, crs_homolsine, always_xy=True)
    transformer_epsg3035 = Transformer.from_crs(4326, crs_epsg3035, always_xy=True)
    coords_homolsine = transformer.transform(lon, lat)
    coords_epsg3035 = transformer_epsg3035.transform(lon, lat)

    ## SoilGrids
    # ----------
    # Take bounding box of 4 km around site
    dist = 2000  # m
    ds_list = []
    for id in mean_names + uncertainty_names:
        # Download raw data
        out_path_tif = output_folder / (id + ".tif")
        out_path_nc = output_folder / (id + ".nc")
        data = soil_grids.get_coverage_data(
            service_id="clay",
            coverage_id=id,
            west=coords_homolsine[0] - dist,
            east=coords_homolsine[0] + dist,
            south=coords_homolsine[1] - dist,
            north=coords_homolsine[1] + dist,
            crs="urn:ogc:def:crs:EPSG::152160",  # Homolsine;
            output=str(out_path_tif),
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
    ds_cube.to_netcdf(output_folder / (site + "_" + var_of_interest + ".nc"))

    ## European Soil Database
    # -----------------------
    # Select 1 x 1 km pixel closest to site
    da_root_depth_site = da_root_depth.sel(
        x=coords_epsg3035[0], y=coords_epsg3035[1], method="nearest"
    )
    da_root_depth_site.to_netcdf(netcdf_dir / (site + "_root_depth.nc"))
