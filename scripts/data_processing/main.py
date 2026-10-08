# %% Imports
from conf import LOG_FORMAT, conf_module, logger
from eddy_covariance_process import main as process_eddy_covariance
from land_cover_translation_process import main as process_land_cover_translation
from soil_moisture_fluxnet_process import main as process_soil_moisture_fluxnet
from soil_process import main as process_soil
from vegetation_process import main as process_vegetation

# %% Run all processing steps, logging to file
if __name__ == "__main__":
    conf_module.dataprodir.mkdir(exist_ok=True, parents=True)
    logger.add(
        str(conf_module.dataprodir / "data_process_log_{time:YYYY-MM-DD_HH-mm}.txt"),
        format=LOG_FORMAT,
        enqueue=True,
        mode="w",
    )

    logger.info("Starting data processing scripts")
    logger.info("1: Processing eddy covariance data")
    process_eddy_covariance()
    logger.info("2: Processing soil data")
    process_soil()
    logger.info("3: Processing vegetation data")
    process_vegetation()
    logger.info("4: Processing soil moisture data from flux towers")
    process_soil_moisture_fluxnet()
    logger.info("5: Processing land cover data for IBGP with BATS mapping")
    process_land_cover_translation()
