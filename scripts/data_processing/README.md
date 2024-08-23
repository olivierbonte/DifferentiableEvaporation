# Data processing

The same computational environment and way of execution files as explained for the data downloading can be applied (see explanation [here](../README.md)). The `conf.py`provided here only serves as a link to the true configuration file under `data_download`. For the full processing chain, `main.py` can be run, which executes following files:
1. `eddy_covariance.py`: Remove periods of low quality/missing data, so that only full years of grap-free forcing data (variables defined in `../data_download/conf.py`) remain.
2. `soil.py`: Combine the information from SoilGrids, HiHydroSoil and the European Soil Database Derived data in one datacube.
3. `vegetation.py`: Derive the fraction of vegetation cover, currently providing 2 options (based on [Zhong et al. (2022)](https://doi.org/10.5194/hess-26-5647-2022) and [van Dijk & Bruijnzeel (2001)](https://doi.org/10.1016/S0022-1694(01)00392-4)).

Output can be found in `data/exp_pro`.
