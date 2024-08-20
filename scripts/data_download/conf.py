# %% Configuration file for the data downloads
import pathlib

# %% folder paths
conf_file_path = pathlib.Path(__file__)
root_project_path = conf_file_path.parent.parent.parent
datadir = root_project_path / "data"
datarawdir = datadir / "exp_raw"

soilgrids_dir = datarawdir / "soilgrids"
hihydrosoil_dir = datarawdir / "hihydrosoil"
ec_dir = datarawdir / "eddy_covariance"
veg_dir = datarawdir / "vegetation"

soil_pro_dir = datadir / "exp_pro" / "soil"
veg_pro_dir = datadir / "exp_pro" / "vegetation"

# %% sites of interest
sites = ["BE-Bra", "ES-LM1", "DE-Hai", "DE-Tha"]
# DE-Tha, NL-Loo also potentially a good site
