# %% Imports
import logging

import numpy as np
import pandas as pd
import pytz
import xarray as xr
from conf import conf_module

# %% Process soil moisture data
date_format = "%Y%m%d%H%M"
for site in conf_module.sites:
    print(f"site in progress: {site}")

    # FLUXNET2015 data
    site_fluxnet_folder = list(
        conf_module.fluxnet_dir.glob("*" + site + "*FLUXNET2015*")
    )[0]
    df_fluxnet_hh = pd.read_csv(
        list(site_fluxnet_folder.iterdir())[0],
        na_values=-9999,
        index_col=0,
        date_format=date_format,
    )

    # More recent data in FLUXNET format
    # Further info on the format in supplementary material of https://doi.org/10.1038/s41597-020-0534-3
    site_l2_folder = list(conf_module.fluxnet_dir.glob("*" + site + "*ARCHIVE_L2"))[0]
    site_l2_data_file = list(site_l2_folder.glob("*" + site + "_FLUXNET_HH_L2.csv"))[0]
    df_l2_hh = pd.read_csv(
        site_l2_data_file, na_values=-9999, index_col=0, date_format=date_format
    )

    # Combined dataset
    if df_fluxnet_hh.index[-1] < df_l2_hh.index[1]:
        logging.warning(
            f"The Fluxnet dataset (ending on {df_fluxnet_hh.index[-1]})"
            "and the L2 archive dataset (starting on {df_fluxnet_hh.index[1]}"
            "do not overlap in time)"
        )
    df_combined_hh = pd.concat([df_fluxnet_hh, df_l2_hh])
    df_combined_hh["TIMESTAMP_END"] = pd.to_datetime(
        df_combined_hh["TIMESTAMP_END"], format=date_format
    )
    # Sort index and take average if duplicate timestamps
    df_combined_hh = df_combined_hh.sort_index()
    df_combined_hh = df_combined_hh.groupby(df_combined_hh.index).mean()
    # Select columns with SM data
    col_bool = df_combined_hh.columns.str.contains(
        "SWC"
    )  # | df_combined_hh.columns.str.contains("TIMESTAMP")
    df_sel = df_combined_hh.iloc[:, col_bool]
    # Set data where QC > 2 (this means poor quality data) to NaN
    swc_columns = df_sel.columns[
        df_sel.columns.str.contains(r"^SWC_F_MDS_\d+$", regex=True)
    ]
    for swc_col in swc_columns:
        qc_col = swc_col + "_QC"
        df_sel.loc[:, swc_col] = df_sel[swc_col].where(df_sel[qc_col] <= 2)

    # Metadata files
    # Variable metadata
    var_metadata_file = list(
        site_l2_folder.glob("*" + site + "_VARINFO_FLUXNET_HH_L2.csv")
    )[0]
    df_var_meta = pd.read_csv(var_metadata_file)
    df_var_meta_sel = df_var_meta[["GROUP_ID", "VARIABLE", "DATAVALUE"]]
    df_var_meta_wide = df_var_meta_sel.pivot(
        index="GROUP_ID", columns="VARIABLE", values="DATAVALUE"
    )
    df_var_meta_wide_swc = df_var_meta_wide[
        df_var_meta_wide["VAR_INFO_VARNAME"].str.contains(
            r"^SWC_F_MDS_\d+$", regex=True
        )
    ]
    # Site metadata
    time_metadata = xr.open_dataset(conf_module.ec_pro_dir / (site + ".nc")).time.attrs
    site_metadata_file = list(site_l2_folder.glob("*" + site + "_SITEINFO_L2.csv"))[0]
    df_site_meta = pd.read_csv(site_metadata_file)
    utc_offset = df_site_meta[df_site_meta["VARIABLE"].str.contains("UTC_OFFSET")][
        "DATAVALUE"
    ].values[0]
    time_metadata["UTC_offset"] = int(utc_offset)
    # Convert to format suitable for use in pytz.timezone(utc_offset_str)
    # !!CAUTION: Etc/GMT-x is equal to UTC+x!!
    if int(utc_offset) >= 0:
        time_metadata["time_zone"] = "Etc/GMT-" + utc_offset
    else:
        time_metadata["time_zone"] = "Etc/GMT+" - utc_offset

    # Convert to .nc
    df_sel_wide = df_sel.reset_index().melt(
        id_vars=["TIMESTAMP_START"], var_name="depth", value_name="SWC"
    )
    df_sel_wide["variable"] = np.where(
        df_sel_wide["depth"].str.endswith("QC"), "SWC_QC", "SWC"
    )
    df_sel_wide["depth"] = df_sel_wide["depth"].str.strip("_QC")
    depth_series = (
        df_var_meta_wide_swc.reset_index()
        .set_index("VAR_INFO_VARNAME")["VAR_INFO_HEIGHT"]
        .astype(float)
    )  # Links variable name with sensor depth
    depth_series = depth_series * -1  # Define depth below ground positive
    df_sel_wide["depth"] = df_sel_wide["depth"].replace(depth_series.to_dict())
    ds_swc = (
        df_sel_wide.set_index(["TIMESTAMP_START", "depth"])
        .pivot(values="SWC", columns="variable")
        .to_xarray()
    )

# %%
