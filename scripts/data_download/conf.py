# %% Configuration file for the data downloads
import pathlib

# %% Google Earth Engine project ID
gee_project_id = "ee-bonteolivier15"  # Replace by your own project ID!

# %% sites of interest
sites = ["BE-Bra", "DE-Hai", "DE-Tha"]
# DE-Tha, NL-Loo also potentially a good site

# %% List of input variables needed for the model (for quality check)
input_variables = ["Tair", "SWdown", "VPD", "Psurf", "Rnet", "Qle", "Qle_cor", "Qh"]

# %% List of variables of interest from SoilGrids
soil_vars_of_interest = ["clay", "sand"]

# %% folder paths
conf_file_path = pathlib.Path(__file__)
root_project_path = conf_file_path.parent.parent.parent
datadir = root_project_path / "data"
datarawdir = datadir / "exp_raw"
dataprodir = datadir / "exp_pro"

soilgrids_dir = datarawdir / "soilgrids"
hihydrosoil_dir = datarawdir / "hihydrosoil"
ec_dir = datarawdir / "eddy_covariance"
veg_dir = datarawdir / "vegetation"
essd_dir = datarawdir / "STU_EU_Layers"  # European Soil Database Derived data
fluxnet_dir = datarawdir / "fluxnet_for_soil_moisture"
stocker_dir = datarawdir / "root_depth_stocker"

ec_pro_dir = dataprodir / "eddy_covariance"
soil_pro_dir = dataprodir / "soil"
