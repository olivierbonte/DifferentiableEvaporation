## Vegetation parameters
table_metadata = "Vegetation parameters for Jarvis-Stewart. \
Origin: Table 8.1 of IFS Cy47r3 documentation (Physical processes), \
 see: https://doi.org/10.21957/eyrpir4vj"

df_veg = DataFrame(;
    vegetation_type=[
        "Crops, mixed farming",
        "Short grass",
        "Evergreen needleleaf trees",
        "Deciduous needleleaf trees",
        "Deciduous broadleaf trees",
        "Evergreen broadleaf trees",
        "Tall grass",
        "Desert",
        "Tundra",
        "Irrigated crops",
        "Semidesert",
        "Ice caps and glaciers",
        "Bogs and marshes",
        "Inland water",
        "Ocean",
        "Evergreen shrubs",
        "Deciduous shrubs",
        "Mixed forest/woodland",
        "Interrupted forest",
        "Water and land mixtures",
    ],
    r_smin=[
        100.0,
        100.0,
        250.0,
        250.0,
        175.0,
        240.0,
        100.0,
        250.0,
        80.0,
        190.0,
        150.0,
        missing,
        240.0,
        missing,
        missing,
        225.0,
        225.0,
        250.0,
        175.0,
        150.0,
    ],
    f_veg=[
        0.90,
        0.85,
        0.90,
        0.90,
        0.90,
        0.99,
        0.70,
        0.0,
        0.50,
        0.90,
        0.10,
        missing,
        0.60,
        missing,
        missing,
        0.50,
        0.50,
        0.90,
        0.90,
        0.60,
    ],
    g_d=[
        0.0,
        0.0,
        0.03,
        0.03,
        0.03,
        0.03,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        missing,
        0.0,
        missing,
        missing,
        0.0,
        0.0,
        0.03,
        0.03,
        0.0,
    ],
)
#https://dataframes.juliadata.org/stable/lib/metadata/#Examples 
metadata!(df_veg, "caption", table_metadata; style=:note)
colmetadata!(df_veg, :r_smin, "label", "Minimum stomatal resistance [s/m]"; style=:note)
colmetadata!(df_veg, :f_veg, "label", "Static vegetation fraction cover [-]"; style=:note)
colmetadata!(
    df_veg,
    :g_d,
    "label",
    "Coefficient relating vapour pressure deficit to somatal resistance hPa⁻¹";
    style=:note,
)