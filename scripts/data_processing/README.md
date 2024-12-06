# Data processing

The same computational environment and way of execution files as explained for the data downloading can be applied (see explanation [here](../README.md)). The `conf.py`provided here only serves as a link to the true configuration file under `data_download`. For the full processing chain, `main.py` can be run, which executes following files:

1. `eddy_covariance.py`: Remove periods of low quality/missing data, so that only full years of gap-free forcing data (variables defined in `../data_download/conf.py`) remain. Additionally, it further trims the data to only maintain periods where the quality of both non-corrected and corrected latent heat flux is sufficient (as defined by point 1 to 3 on the [FluxDataKit website](https://geco-bern.github.io/FluxDataKit/articles/04_data_use.html#time-series))
2. `soil.py`: Combine the information from SoilGrids, HiHydroSoil and the European Soil Database Derived data in one datacube.
3. `vegetation.py`: Derive the fraction of vegetation cover, currently providing 2 options (based on [Zhong et al. (2022)](https://doi.org/10.5194/hess-26-5647-2022) and [van Dijk & Bruijnzeel (2001)](https://doi.org/10.1016/S0022-1694(01)00392-4)).
4. `soil_moisture_fluxnet.py`: The soil moisture data in `.csv` format (sources described in [download README](../data_download/README.md)) are converted to `.nc` and combined with metadata describing depth of soil moisture sensor. Data is merged with the data cube created in `eddy_covariance.py`. 

Output can be found in `data/exp_pro`.