# %% Imports
import logging
from glob import glob

import hvplot.xarray
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from conf import conf_module

conf_module.ec_pro_dir.mkdir(exist_ok=True)

# %% Remove periods of bad quality data, only keep continuous long sequences
df_quality = pd.read_csv(
    conf_module.ec_dir / "fdk_site_fullyearsequence.csv", index_col=0
)
for site in conf_module.sites:
    print(f"site in progress: {site}")
    # Quality of LE products check
    df_quality_site = df_quality.loc[[site]]
    start_year = int(
        max(
            df_quality_site["year_start_le"].values[0],
            df_quality_site["year_start_lecorr"].values[0],
        )
    )
    end_year = int(
        min(
            df_quality_site["year_end_le"].values[0],
            df_quality_site["year_end_lecorr"].values[0],
        )
    )
    files = glob(str(conf_module.ec_dir / ("*" + site + "*FLUXDATAKIT*.nc")))
    ds_ec = xr.open_mfdataset(files)
    ds_ec_sel = ds_ec.sel(time=slice(str(start_year), str(end_year)))

    # Quality of intput variable to model check
    for var in conf_module.input_variables:
        if ds_ec_sel[var].isnull().any():
            print(f"{var} misses data")
            # Split in groups seperated by Nans
            slice_list = np.ma.clump_unmasked(
                np.ma.masked_invalid(ds_ec_sel[var].values)
            )
            # Determine length of each group
            length_list = [tmp_slice.stop - tmp_slice.start for tmp_slice in slice_list]
            # Select group of longest length
            selected_slice = slice_list[np.argmax(length_list)]
            print(slice_list)
            print(length_list)
            ds_ec_sel = ds_ec_sel.isel(time=selected_slice)
            # fig, ax = plt.subplots()
            # ds_ec_sel.Rnet.plot(ax=ax)
            print(
                f"""Selected data: Start = {ds_ec_sel.time[0].values},
                    End = {ds_ec_sel.time[-1].values}"""
            )
    # Further trim data to only keep full years
    if ds_ec_sel.time[0].dt.dayofyear > 1:
        start_year_trim = int(ds_ec_sel.time[0].dt.year.values + 1)
        ds_ec_sel = ds_ec_sel.isel(time=ds_ec_sel.time.dt.year >= start_year_trim)
    if ds_ec_sel.time[-1].dt.dayofyear < 365:
        end_year_trim = int(ds_ec_sel.time[-1].dt.year.values - 1)
        ds_ec_sel = ds_ec_sel.sel(time=ds_ec_sel.time.dt.year <= end_year_trim)
    print(
        f"""Selected data (full years): Start = {ds_ec_sel.time[0].values},
            End = {ds_ec_sel.time[-1].values} \n"""
    )

    # Sanity check: make sure that no more water is evaporated than it rains
    # Inspired by check in PLUMBER 2 model intercomparison, see Abramowitz et al. (2024)
    # https://doi.org/10.5194/egusphere-2023-3084
    latent_heat_J_per_kg = lambda x: 1000 * (
        2500.8 - 2.36 * x + 0.0016 * x**2 - 0.00006 * x**3
    )  # https://en.wikipedia.org/wiki/Latent_heat
    LE_kg_per_m2s = ds_ec_sel.Qle / latent_heat_J_per_kg(ds_ec_sel.Tair - 273.15)
    if LE_kg_per_m2s.sum() > ds_ec_sel.Precip.sum():
        logging.warning(
            f"For {site}, the total amount of water evaporated is higher than the rainfall"
        )

    # Write to disk
    ds_ec_sel.to_netcdf(conf_module.ec_pro_dir / (site + ".nc"))

# %%