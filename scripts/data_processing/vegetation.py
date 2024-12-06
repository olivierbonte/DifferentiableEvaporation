# %% Imports
from glob import glob

import numpy as np
import xarray as xr
from conf import conf_module

conf_module.veg_pro_dir.mkdir(exist_ok=True)

# %%
for site in conf_module.sites:
    print(f"site in progress: {site}")
    ds_veg = xr.open_dataset(conf_module.veg_dir / (site + "_horizontal_agg.nc"))
    ds_met = xr.open_dataset(
        glob(str(conf_module.ec_dir / ("*" + site + "*Met.nc")))[0]
    )
    ds_ec = xr.open_dataset(conf_module.ec_pro_dir / f"{site}.nc")
    ## Option 1: Based on Zhong et al. (2022) https://doi.org/10.5194/hess-26-5647-2022
    ## VCF and FPAR dependence
    total_veg_cover = (
        ds_veg.sel(band="Percent_Tree_Cover")
        + ds_veg.sel(band="Percent_NonTree_Vegetation")
    ).to_dataframe()
    total_veg_cover["year"] = total_veg_cover.index.year

    df_fpar = ds_met["FPAR"].isel(x=0, y=0).to_dataframe()
    # df_fpar_climatology = df_fpar.groupby(df_fpar.index.dayofyear).mean()
    df_fpar["FPAR_mean"] = df_fpar.groupby(df_fpar.index.year).transform(
        lambda x: x.mean()
    )["FPAR"]

    df_fpar["year"] = df_fpar.index.year
    df_fpar["cyclic_trend"] = df_fpar["FPAR"] / df_fpar["FPAR_mean"]

    # Join FPAR and veg frac based on year: VCF is yearly ->
    # Applicable for every day of that year
    df_join = df_fpar[["FPAR", "FPAR_mean", "year"]].merge(
        total_veg_cover[["year", "MODIS_VCF"]], how="left", on="year"
    )
    df_join = df_join.set_index(ds_met.isel(x=0, y=0).to_dataframe().index)
    df_join["f_veg"] = df_join["MODIS_VCF"] * df_join["FPAR"] / df_join["FPAR_mean"]
    # Correction 1: limit to max 1
    # Analogous to Trautman et al. (2022): https://doi.org/10.5194/hess-26-1089-2022
    df_join["f_veg"] = np.minimum(df_join["f_veg"], 100)

    ## Option 2: Based on LAI (Lambert-Beer), see van Dijk & Bruijnzeel (2001)
    ## https://doi.org/10.1016/S0022-1694(01)00392-4
    ## k_ext from PML, see https://github.com/gee-hydro/gee_PML/blob/stable/src/pkg_PMLV2_v0.1.5.js#L422C32-L422C134
    k_ext = 0.7
    df_join["f_veg_vdb"] = (
        1 - np.exp(-k_ext * ds_met["LAI"].isel(x=0, y=0).to_dataframe()["LAI"].values)
    ) * 100

    ## To xarray and add metadata
    ds_join = df_join.to_xarray()
    # ds_join["FPAR"].attrs = ds_met["FPAR"].attrs
    ds_join["MODIS_VCF"].attrs = ds_veg["MODIS_VCF"].attrs
    f_veg_dict = {
        "units": "-",
        "long_name": "Fraction of surface covered by vegetation",
    }
    ds_join["MODIS_VCF"].attrs.update(f_veg_dict)
    f_veg_dict.update(
        {
            "Description": """Method used from Zhong et al. (2022), 
            see https://doi.org/10.5194/hess-26-1089-2022"""
        }
    )
    ds_join["f_veg"].attrs = f_veg_dict
    f_veg_dict["Description"] = (
        """Method used from van Dijk & Bruijnzeel (2001) (LAI-based), 
        see https://doi.org/10.1016/S0022-1694(01)00392-4"""
    )
    ds_join["f_veg_vdb"].attrs = f_veg_dict
    ds_join = ds_join.astype(np.float32)

    ## Append to the flux tower datacube
    ds_ec.close()
    ds_merged = xr.merge(
        [ds_ec, ds_join.drop_vars(["year", "FPAR", "FPAR_mean"])], join="left"
    )
    ds_merged.to_netcdf(conf_module.ec_pro_dir / f"{site}.nc", mode="a")


# %%
