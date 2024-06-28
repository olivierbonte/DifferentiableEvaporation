import subprocess
subprocess.run("earthengine authenticate --quiet", shell=True, capture_output=True, text=True)
import xarray as xr
import ee  
import os
ee.Initialize(opt_url = "https://earthengine-highvolume.googleapis.com")

datadir = os.path.join("..","data")
datarawdir = os.path.join(datadir,"exp_raw")
ds_meteo_BE_Bra = xr.open_dataset(
    os.path.join(datarawdir,"BE-Bra_2004-2014_FLUXNET2015_Flux.nc")
)

lon_Bra = float(ds_meteo_BE_Bra["longitude"].values[0][0])
lat_Bra = float(ds_meteo_BE_Bra["latitude"].values[0][0])
# point_BE_Bra = ee.Geometry.Point(
#     [float(ds_meteo_BE_Bra["longitude"].values[0][0]),
#     float(ds_meteo_BE_Bra["latitude"].values[0][0])],
#     "EPSG:4326"
# )
test = ee.Geometry.Rectangle(lon_Bra - 1, lat_Bra - 1, lon_Bra + 1, lat_Bra + 1)
ic_clay = ee.ImageCollection(
    ee.Image("projects/soilgrids-isric/clay_mean")
    )
# ds = xr.open_dataset(
#     ic_clay, engine = "ee",
#     geometry = test,
#     crs = "EPSG:4326",
#     resolution = 250
# )

## Different approach: use WCS 
# (see https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/03-WCS-2.0.ipynb?ref_type=heads)
from owslib.wcs import WebCoverageService
wcs = WebCoverageService("https://maps.isric.org/mapserv?map=/map/clay.map",
                         version = '2.0.1')
print(wcs.contents)

#Take mean (better metrics, see p. 224 of https://doi.org/10.5194/soil-7-217-2021)
mean_names = [k for k in wcs.contents.keys() if k.find("mean") != -1]
uncertainty_names = [k for k in wcs.contents.keys() if k.find("uncertainty") != -1]

from pyproj import Transformer, CRS
#Homolsine projection is used (see https://www.isric.org/explore/soilgrids/faq-soilgrids#How_can_I_use_the_Homolosine_projection)
#crs_152160 = CRS.from_proj4("+proj=igh +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
crs_opengis = "http://www.opengis.net/def/crs/EPSG/0/152160" 
crs_homolsine = CRS.from_wkt(
    """
PROJCS["Homolosine", 
    GEOGCS["WGS 84", 
        DATUM["WGS_1984", 
            SPHEROID["WGS 84",6378137,298.257223563, 
                AUTHORITY["EPSG","7030"]], 
   AUTHORITY["EPSG","6326"]], 
        PRIMEM["Greenwich",0, 
            AUTHORITY["EPSG","8901"]], 
        UNIT["degree",0.0174532925199433, 
            AUTHORITY["EPSG","9122"]], 
        AUTHORITY["EPSG","4326"]], 
    PROJECTION["Interrupted_Goode_Homolosine"], 
    UNIT["Meter",1]]
"""
)
transformer = Transformer.from_crs(4326, crs_homolsine, always_xy=True) #lon/lat
coords_Bra_homolsine = transformer.transform(lon_Bra, lat_Bra)
#Take bounding box of 4 km around the site 
dist = 2000 #m
subsets = [('X', int(coords_Bra_homolsine[0] - dist), int(coords_Bra_homolsine[0] - dist)),
           ('Y', int(coords_Bra_homolsine[1] - dist), int(coords_Bra_homolsine[1] + dist))]
#Test data retrieval
response = wcs.getCoverage(
    identifier = [mean_names[0]],
    crs = crs_opengis,
    subsets = subsets,
    resx = 250, resy = 250,
    format = wcs.contents[mean_names[0]].supportedFormats[0]
)
## v1.0.0
wcs = WebCoverageService("https://maps.isric.org/mapserv?map=/map/clay.map",
                         version = '1.0.0')
bbox = (coords_Bra_homolsine[0] - dist, coords_Bra_homolsine[1] - dist,
        coords_Bra_homolsine[0] + dist, coords_Bra_homolsine[1] + dist)
response = wcs.getCoverage(
    identifier=mean_names[0],
    crs='urn:ogc:def:crs:EPSG::152160',
    bbox=bbox, 
    resx=250, resy=250, 
    format='GEOTIFF_INT16')
with open('test_clay_0-5_mean.tif', 'wb') as file:
    file.write(response.read())
ds_Bra = rioxarray.open_rasterio("test_clay_0-5_mean.tif", mask_and_scale = True)

##TEST
wcs = WebCoverageService('http://maps.isric.org/mapserv?map=/map/phh2o.map',
                         version='2.0.1')
cov_id = 'phh2o_0-5cm_mean'
ph_0_5 = wcs.contents[cov_id]
ph_0_5.supportedFormats 
subsets = [('X', -1784000, -1140000), ('Y', 1356000, 1863000)]
subsets = [('X', coords_Bra_homolsine[0] - dist, coords_Bra_homolsine[0] - dist),
           ('Y', coords_Bra_homolsine[1] - dist, coords_Bra_homolsine[1] + dist)]
crs = "http://www.opengis.net/def/crs/EPSG/0/152160"
response = wcs.getCoverage(
    identifier=[cov_id], 
    crs=crs,
    subsets=subsets, 
    resx=250, resy=250, 
    format=ph_0_5.supportedFormats[0])
with open('Bra-5_mean.tif', 'wb') as file:
    file.write(response.read())
import rioxarray
# ds_senegal = xr.open_dataset("Senegal_pH_0-5_mean.tif", engine = "rasterio")
ds_senegal = rioxarray.open_rasterio("Senegal_pH_0-5_mean.tif", mask_and_scale=True)

##testing soilgrids package!
from soilgrids import SoilGrids