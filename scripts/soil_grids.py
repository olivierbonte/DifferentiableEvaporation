#%%
import rioxarray
import os
import subprocess
import xarray as xr
from pyproj import CRS, Transformer
from soilgrids import SoilGrids
from owslib.wcs import WebCoverageService

datadir = os.path.join("..","data")
data_ec_dir = os.path.join(datadir, "exp_raw","eddy_covariance")

#%% Clay information 
var_of_interest = 'clay'
soil_grids = SoilGrids()
url = soil_grids.MAP_SERVICES[var_of_interest]['link']
wcs = WebCoverageService(url, version = '1.0.0')
mean_names = [k for k in wcs.contents.keys() if k.find("mean") != -1]
uncertainty_names = [k for k in wcs.contents.keys() if k.find("uncertainty") != -1]

# %%Define sites of interest
sites = ["BE-Bra", "ES-LM1"]

# %% Store data for each 
for site in sites:
    files = [k for k in os.listdir(data_ec_dir) if site in k]
    ds_site = xr.open_dataset(os.path.join(data_ec_dir , files[0]))
    output_folder  = os.path.join(datadir,"exp_raw","soilgrids",site)
    output_pro_folder = os.path.join(datadir, "exp_pro","soilgrids",site)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_pro_folder):
        os.makedirs(output_pro_folder)
    
    # EPSG:4326 coordinates 
    lon = float(ds_site["longitude"].values[0][0])
    lat = float(ds_site["latitude"].values[0][0])
    
    # Homolsine projection is used (https://epsg.io/54052 for proj4), see
    # https://www.isric.org/explore/soilgrids/faq-soilgrids#How_can_I_use_the_Homolosine_projection
    crs_homolsine = CRS.from_proj4("+proj=igh +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
    transformer = Transformer.from_crs(4326, crs_homolsine, always_xy=True)
    coords_homolsine = transformer.transform(lon, lat)
    
    # Take bounding box of 4 km around site
    dist = 2000 # m
    ds_list = []
    for id in (mean_names + uncertainty_names):
        out_path_tif = os.path.join(output_folder, id + '.tif')
        out_path_nc = os.path.join(output_folder, id + '.nc')
        data = soil_grids.get_coverage_data(
            service_id="clay",
            coverage_id=id,
            west=coords_homolsine[0] - dist,
            east=coords_homolsine[0] + dist,
            south=coords_homolsine[1] - dist,
            north=coords_homolsine[1] + dist,
            crs="urn:ogc:def:crs:EPSG::152160", # Homolsine;
            output = out_path_tif
        )
        # #Convert to NetCDF with GDAL 
        # subprocess.run(["gdal_translate", "-of" ,"NetCDF",
        #                 os.path.abspath(out_path_tif),
        #                 os.path.abspath(out_path_nc)],
        #                 capture_output=True)
        ds_temp = rioxarray.open_rasterio(out_path_tif,mask_and_scale=True)
        ds_temp.name = id
        ds_list.append(
            ds_temp.isel(band=0)
        )
    #Make 1 datacube per site 
    ds_cube = xr.merge(ds_list)
    ds_cube.to_netcdf(os.path.join(output_pro_folder,site + '_' + var_of_interest + ".nc"))

    

