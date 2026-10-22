# Synthetic diurnal forcing and parameters for the 1-day reference run.
# Needs `model_parameters` from src/site_data.jl.
# Time `t` is local standard time in seconds since 1970-01-01 (as for the site forcing).
# Units follow EvaporationModel: K, Pa, m/s, W/m², kg m⁻² s⁻¹; R_n is positive downward.

const SYNTHETIC_START = DateTime(2010, 7, 1)
const SYNTHETIC_LON = 4.52 # °E, BE-Bra
const SYNTHETIC_UTC_OFFSET = 1 # h

# Diurnal shape: -1 at midnight, +1 at noon (local standard time)
diurnal_cycle(t::T) where {T} = -cos(T(2π) * mod(t, T(86400)) / T(86400))

# Shower on/off: jumps of the forcing, with `tstops` just before them
const SHOWER_START_HR = 6
const SHOWER_END_HR = 8
P_synthetic(t::T) where {T} = ifelse(
    T(SHOWER_START_HR * 3600) <= mod(t, T(86400)) < T(SHOWER_END_HR * 3600),
    T(2.5e-4),
    zero(T)
) # 0.9 mm/h shower 06-08 h
T_a_synthetic(t::T) where {T} = T(288.0) + T(5.0) * diurnal_cycle(t)
u_a_synthetic(t::T) where {T} = T(3.0)
p_a_synthetic(t::T) where {T} = T(101325.0)
VPD_a_synthetic(t::T) where {T} = T(800.0) + T(600.0) * diurnal_cycle(t)
SW_in_synthetic(t::T) where {T} = max(T(800.0) * diurnal_cycle(t), zero(T))
R_n_synthetic(t::T) where {T} = T(0.7) * SW_in_synthetic(t) - T(40.0)
LAI_synthetic(t::T) where {T} = T(3.0)

synthetic_forcing() = (
    P=P_synthetic,
    T_a=T_a_synthetic,
    u_a=u_a_synthetic,
    p_a=p_a_synthetic,
    VPD_a=VPD_a_synthetic,
    SW_in=SW_in_synthetic,
    R_n=R_n_synthetic,
    LAI=LAI_synthetic,
    lon=SYNTHETIC_LON,
    utc_offset=SYNTHETIC_UTC_OFFSET,
)

# Loam soil texture [%]; all soil properties follow from it
const SYNTHETIC_CLAY = 19.0
const SYNTHETIC_SAND = 40.0

# ISBA pedotransfer functions for the soil moisture thresholds [m³ m⁻³]
# (Noilhan & Lacarrère, 1995, Appendix A.4: continuous formulation of the soil parameters)
w_sat_isba(x_sand) = (-1.08 * x_sand + 494.305) * 1e-3
w_wp_isba(x_clay) = 37.1342e-3 * sqrt(x_clay)
w_fc_isba(x_clay) = 89.0467e-3 * x_clay^0.3496

# Model parameters with BE-Bra heights and vegetation (evergreen needleleaf trees)
synthetic_parameters() = model_parameters(;
    x_clay=SYNTHETIC_CLAY,
    w_sat=w_sat_isba(SYNTHETIC_SAND),
    w_wp=w_wp_isba(SYNTHETIC_CLAY),
    w_fc=w_fc_isba(SYNTHETIC_CLAY),
    # Residual soil moisture: only used by the Martens17 soil evaporation efficiency,
    # not by the default Pielke92 one, so not derived from the texture
    w_res=0.04,
    d_2=1.3,
    h=21.0,
    z_obs=39.0,
)

synthetic_case() = (
    forcings=synthetic_forcing(),
    p=synthetic_parameters(),
    u0=[0.2, 0.2, 0.001],
    tspan=(datetime2unix(SYNTHETIC_START), datetime2unix(SYNTHETIC_START + Day(1))),
    tstops=tstops_before(datetime2unix.(SYNTHETIC_START .+ Hour.([SHOWER_START_HR, SHOWER_END_HR]))),
)
