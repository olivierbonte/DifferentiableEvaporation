# %% Imports
from conf import LOG_FORMAT, datarawdir
from eddy_covariance_download import main as download_eddy_covariance
from hihydrosoil_download import main as download_hihydrosoil
from land_cover_translation_download import main as download_land_cover_translation
from loguru import logger
from root_depth_stocker_download import main as download_root_depth_stocker
from soil_grids_download import main as download_soil_grids
from soil_moisture_fluxnet_download import main as download_soil_moisture_fluxnet
from vegetation_vcf_download import main as download_vegetation_vcf

# %% Run all downloads, logging to file
if __name__ == "__main__":
    datarawdir.mkdir(exist_ok=True, parents=True)
    logger.add(
        str(datarawdir / "data_download_log_{time:YYYY-MM-DD_HH-mm}.txt"),
        format=LOG_FORMAT,
        enqueue=True,
        mode="w",
    )

    logger.info("Starting data download script")
    logger.info("1: Downloading eddy covariance data")
    download_eddy_covariance()
    logger.info("2: Downloading data from SoilGrids")
    download_soil_grids()
    logger.info("3: Download data from HiHydroSoil")
    download_hihydrosoil()
    logger.info("4: Download data on vegetation cover fractions")
    download_vegetation_vcf()
    logger.info("5: Downloading auxiliary data from flux towers")
    download_soil_moisture_fluxnet()
    logger.info("6: Downloading root depth data from Stocker et al. (2023)")
    download_root_depth_stocker()
    logger.info("7: Downloading land cover translation data from GLCC")
    download_land_cover_translation()
