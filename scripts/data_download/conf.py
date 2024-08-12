# %% Configuration file for the data downloads
import os

# %% folder paths
datadir = os.path.join("..", "..", "data")
datarawdir = os.path.join(datadir, "exp_raw")

soilgrids_dir = os.path.join(datarawdir, "soilgrids")
hihydrosoil_dir = os.path.join(datarawdir, "hihydrosoil")
ec_dir = os.path.join(datarawdir, "eddy_covariance")
veg_dir = os.path.join(datarawdir, "vegetation")

soil_pro_dir = os.path.join(datadir, "exp_pro", "soil")
veg_pro_dir = os.path.join(datadir, "exp_pro", "vegetation")

# %% sites of interest
sites = ["BE-Bra", "ES-LM1"]
