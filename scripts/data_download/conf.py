#%% Configuration file for the data downloads
import os

#%% folder paths
datadir = os.path.join("..","..","data")
datarawdir = os.path.join(datadir,"exp_raw")

soilgrids_dir = os.path.join(datarawdir, "soilgrids")
hihydrosoil_dir = os.path.join(datarawdir, "hihydrosoil")
ec_dir = os.path.join(datarawdir,"eddy_covariance")

#%% sites of interest
sites = ["BE-Bra","ES-LM1"]