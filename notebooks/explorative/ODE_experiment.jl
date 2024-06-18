## Idea of defining the ODE model here
using DrWatson
@quickactivate "DifferentiableEvaporation"
using EvaporationModel, Bigleaf

#Define parameters for initial experiment
#Inspired by CLASS: short grass, sandy loam (p.129-131)
#Table 10.6: Short grass
z0m = 0.02 #m
z0h = Bigleaf.roughness_z0h(z0m, log(10)) #for now: use /10
d = 10*z0m*2/3 #assuming displacemengt height is 2/3h and z0m = 1/10h
LAI = 2.0 # [-]
gD = 0.0
rsmin = 110.0 # s/m
fveg = 0.85

#Table 10.7: Sandy Loam
w_sat = 0.472 # [-]
w_fc = 0.323 # [-]
w_wp = 0.171 # [-]
C1sat = 0.132 # [-]
C2ref = 1.8 # [-]
a_ch = 0.219
b_ch = 4.9
p_ch = 4.0

Bigleaf.compute_Ram()

EvaporationModel.jarvis_stewart()

# using ComponentArrays, Parameters
# test = ComponentArray(a = [3,4], b = [5,6])
# @unpack a,b = test