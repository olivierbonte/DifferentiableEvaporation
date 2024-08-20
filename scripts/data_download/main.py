import logging
import subprocess

logging.basicConfig(level=logging.INFO)

logging.info("\n Starting data download script \n")
logging.info("\n 1: Downloading eddy covaraince data \n")
subprocess.run("python eddy_covariance.py", shell=True)
logging.info("\n 2: Downloading data from SoilGrids \n")
subprocess.run("python soil_grids.py", shell=True)
logging.info("\n 3: Download data from HiHydroSoil \n")
subprocess.run("python hihydrosoil.py", shell=True)
logging.info("\n 4: Download data on vegetation cover fractions \n")
subprocess.run("python vegetation_VCF.py", shell=True)
