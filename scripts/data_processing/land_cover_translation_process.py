# %% Imports
import pandas as pd
from conf import glcc_dir, glcc_pro_dir

glcc_pro_dir.mkdir(parents=True, exist_ok=True)
# %% Read in the land cover translation data per continent
df_list = []
for file in glcc_dir.glob("*.txt"):
    pd_temp = pd.read_table(file, on_bad_lines="warn")
    # pd_temp["PIXELS"].astype(float)
    pd_temp["PIXELS"] = pd.to_numeric(pd_temp["PIXELS"], errors="coerce")
    df_list.append(pd_temp)
df = pd.concat(df_list, ignore_index=True)

# %% Calculate the fraction of each BATS land cover in each IGBP land cover type
df_sel = df[["IGBP V2.0", "IGBP", "BATS", "BATS V2.0", "PIXELS"]]
df_IGBP_to_BATS = (
    pd.DataFrame(df_sel.groupby(["IGBP V2.0", "IGBP", "BATS", "BATS V2.0"]).sum())
    .reset_index()
    .set_index(["IGBP V2.0", "BATS V2.0"])
)
total_pixels_per_IGBP = pd.DataFrame(
    df_IGBP_to_BATS.groupby("IGBP V2.0")["PIXELS"].sum()
)
df_IGBP_to_BATS["fraction_BATS_in_IGBP"] = (
    df_IGBP_to_BATS[["PIXELS"]] / total_pixels_per_IGBP
)
df_IGBP_to_BATS.to_csv(glcc_pro_dir / "glcc_fraction_BATS_in_IGBP.csv")
