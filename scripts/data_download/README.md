# Data downloading

To download the necessary data, first activate the conda python environment as defined by the `environment.yml` (see main `README.md`).  To modify the sites of interest (currently `BE-Bra` and `ES-LM1`), change the list `sites` in `conf.py`. Then, execute the downloads by running scripts from the `data_download folder` (important for correct paths) in the following order:

1. `eddy_covariance.py`: Download flux, meteorological and LAI data from [FluxDataKit](https://doi.org/10.5281/zenodo.11370417) via Zenodo. To download data from [PLUMBER 2](http://doi.org/10.25914/5fdb0902607e1) via OpenDAP, manually add the url of interest to `urls` in this python file. 
2. `soil_grids.py`: Download the percentage of clay in the soil from [SoilGrids](https://doi.org/10.5194/soil-7-217-2021
). Data is accessed via the WebCoverageService (WCS) of OGC, example notebooks for use with SoilGrids are given [here](https://git.wur.nl/isric/soilgrids/soilgrids.notebooks).
3. `hihydrosoil.py`: Download hydraulic soil properties (for full list, see the python file) from the Google Earth Engine commutni catalog (see [here](https://gee-community-catalog.org/projects/hihydro_soil/)) via [Xee](https://github.com/google/Xee/tree/v0.0.14).