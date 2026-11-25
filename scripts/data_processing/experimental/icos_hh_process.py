# Experimental: model forcing and evaluation data from the half-hourly ICOS FLUXNET
# product (instead of FluxDataKit), for periods after the FluxDataKit record.
# Used for the reproduction of Kronenberg et al. (2026), https://doi.org/10.5194/gmd-19-5261-2026
# %% Imports
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from conf import conf_module, logger
from utils import knn_gapfill

# %% Settings
site = "DE-Tha"
start_year, end_year = (
    2023,
    2024,
)  # 2023 spin-up, 2024 evaluation as in Kronenberg et al. (2026)
utc_offset_hours = 1.0  # ICOS timestamps are in local standard time (UTC+1 for DE-Tha)

# ICOS FLUXNET name: (name in processed FluxDataKit data, conversion to its unit)
variables = {
    "TA_F": ("Tair", lambda x: x + 273.15),  # °C -> K
    "SW_IN_F": ("SWdown", lambda x: x),
    "LW_IN_F": ("LWdown", lambda x: x),
    "VPD_F": ("VPD", lambda x: x),  # hPa
    "PA_F": ("Psurf", lambda x: x * 1000),  # kPa -> Pa
    "P_F": ("Precip", lambda x: x / 1800),  # mm per 30 min -> kg m⁻² s⁻¹
    "WS_F": ("Wind", lambda x: x),
    "NETRAD": ("Rnet", lambda x: x),
    "G_F_MDS": ("Qg", lambda x: x),
    "LE_F_MDS": ("Qle", lambda x: x),
    "LE_F_MDS_QC": ("Qle_qc", lambda x: x),
    "LE_CORR": ("Qle_cor", lambda x: x),
    "LE_CORR_25": ("Qle_cor_25", lambda x: x),
    "LE_CORR_75": ("Qle_cor_75", lambda x: x),
    "H_F_MDS": ("Qh", lambda x: x),
    "H_F_MDS_QC": ("Qh_qc", lambda x: x),
    "H_CORR": ("Qh_cor", lambda x: x),
    "H_CORR_25": ("Qh_cor_25", lambda x: x),
    "H_CORR_75": ("Qh_cor_75", lambda x: x),
}


def main():
    out_dir = conf_module.ec_pro_dir / "experimental"
    out_dir.mkdir(exist_ok=True, parents=True)

    # %% Read the half-hourly ICOS FLUXNET data
    site_folder = list(conf_module.fluxnet_dir.glob("ICOS_" + site + "_FLUXNET*"))[0]
    site_file = list(site_folder.glob("*_HH_*.csv"))[0]
    logger.info(f"Reading {site_file}")
    df = pd.read_csv(
        site_file,
        usecols=["TIMESTAMP_START", *variables],
        na_values=-9999,
        index_col=0,
        date_format="%Y%m%d%H%M",
    )
    df = df.loc[str(start_year) : str(end_year)]

    # %% Gap fill net radiation as in FluxDataKit (KNN on temperature and radiation)
    # Deviation from FluxDataKit: SW_IN_F instead of PPFD as radiation predictor, since
    # PPFD_IN is not gap-filled in the ICOS FLUXNET product and is missing in many of the
    # NETRAD gaps. SW_IN_F as proxy for PPFD.
    logger.info(
        f"Missing net radiation data fraction: {df['NETRAD'].isna().mean():.3f}"
    )
    df["NETRAD"], netrad_filled = knn_gapfill(df, "NETRAD", ["TA_F", "SW_IN_F"], k=5)
    df["NETRAD_QC"] = netrad_filled.astype(float)  # 0 = measured, 1 = KNN gap-filled

    missing = df.columns[df.isna().any()]
    if len(missing) > 0:
        logger.warning(
            f"Variables with missing data after gap filling: {list(missing)}"
        )

    # %% Convert to the names and units of the processed FluxDataKit data
    ds = xr.Dataset(
        {
            name: ("time", conv(df[var].values))
            for var, (name, conv) in variables.items()
        }
    )
    ds["Rnet_qc"] = ("time", df["NETRAD_QC"].values)
    # Gap-filling metadata: attributes of the FluxDataKit data, updated for the KNN filling
    ds_fdk = xr.open_dataset(conf_module.ec_pro_dir / (site + ".nc"))
    ds["Rnet"].attrs = ds_fdk["Rnet"].attrs | {
        "Missing_%": round(100 * float(ds["Rnet"].isnull().mean()), 1),
        "Gap-filled_%": round(100 * float(netrad_filled.mean()), 1),
        "Gapfilling_method": "KNN (k=5) based on SWdown, Tair (fill_netrad of FluxDataKit, "
        "with SW_IN_F instead of PPFD)",
    }
    ds["Rnet_qc"].attrs = ds_fdk["Rnet_qc"].attrs | {
        "long_name": "Rnet quality control flag (0 = measured, 1 = KNN gap-filled)",
    }
    ds = ds.assign_coords(time=df.index.values)
    ds.time.attrs = {
        "standard_name": "time",
        "info": "Time stamp indicates start time",
        "UTC_offset": utc_offset_hours,
        "time_zone": f"Etc/GMT-{int(utc_offset_hours)}",  # Etc/GMT-x is equal to UTC+x
    }

    # %% LAI: day-of-year climatology of the MODIS LAI in the FluxDataKit data
    # Climatology because outside of FLUXDATAKIT covered period
    lai_daily = ds_fdk.LAI.isel(x=0, y=0).resample(time="1D").mean()
    lai_clim = lai_daily.groupby("time.dayofyear").mean()
    ds["LAI"] = ("time", lai_clim.sel(dayofyear=ds.time.dt.dayofyear).values)
    ds["LAI"].attrs["info"] = (
        f"Day-of-year mean of the MODIS LAI of FluxDataKit, "
        f"{int(ds_fdk.time.dt.year[0])}-{int(ds_fdk.time.dt.year[-1])}"
    )
    # Site metadata from FluxDataKit
    for var in [
        "latitude",
        "longitude",
        "reference_height",
        "canopy_height",
        "elevation",
    ]:
        ds[var] = np.float64(ds_fdk[var].values.flat[0])
    ds.attrs = {
        "site_code": site,
        "source": site_file.name,
    }

    # %% Write to disk
    out_file = out_dir / f"{site}_ICOS_HH_{start_year}-{end_year}.nc"
    # Same time encoding as the FluxDataKit data, which YAXArrays decodes
    time_encoding = {
        "dtype": "float64",
        "units": f"seconds since {start_year}-01-01",
        "calendar": "standard",
    }
    ds.to_netcdf(out_file, encoding={"time": time_encoding})
    logger.info(f"Written to {out_file}")


if __name__ == "__main__":
    main()
