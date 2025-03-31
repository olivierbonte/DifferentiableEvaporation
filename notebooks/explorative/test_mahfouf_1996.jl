## Testing the C3 calculation from Mahfouf (1996) 
## https://doi.org/10.1175/1520-0450(1996)035%3C0987:IOGDIA%3E2.0.CO;2
using Plots

τ = 24 * 60 * 60 # [s] # 24h in seconds
# Values below from Table 1 of the paper
b = [4.05, 4.38, 4.9, 5.3, 5.39, 7.12, 7.75, 8.52, 10.4, 10.4, 11.4] # [-]
X_clay = [3.0, 6.0, 9.0, 14.0, 19.0, 28.0, 34.0, 34.0, 43.0, 49.0, 63.0] # [%]
w_sat = [0.395, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.482, 0.482] # [m³/m³]
K_sat = [176.0, 156.3, 34.1, 7.2, 7.0, 6.3, 1.7, 2.5, 2.2, 1.0, 1.3] .* 1e-6 # [m/s]
w_fc = [0.135, 0.150, 0.195, 0.255, 0.240, 0.255, 0.322, 0.325, 0.310, 0.370, 0.367] # [m³/m³]
b = [4.05, 4.38, 4.90, 5.30, 5.39, 7.12, 7.75, 8.52, 10.4, 10.4, 11.4] # [-]
C_3_mahfouf = [1.705, 1.437, 0.510, 0.183, 0.196, 0.206, 0.098, 0.124, 0.139, 0.112, 0.102] # [m]

"""
Equation 15 multiplied by d2
"""
function calculate_C_3(τ, b, K_sat, w_sat, w_fc)
    w_2_star = w_fc + (w_sat - w_fc) / exp(1.0)
    return τ * (2 * b + 2) * K_sat / (w_sat * (w_2_star / w_sat)^(-2 * b - 2) - 1)
end

C_3 = calculate_C_3.(τ, b, K_sat, w_sat, w_fc)
X_clay_plot = range(0.0, maximum(X_clay), length=100)
scatter(X_clay, C_3_mahfouf, label="C3 from Mahfouf (1996)", color=:blue, marker=:circle, markersize=7)
scatter!(X_clay, C_3, label="C3 calculated", color=:red, marker=:circle, markersize=4)
plot!(X_clay_plot, 5.32 * X_clay_plot .^ (-1.042), label="C3 in function of clay percentage")
xlabel!("Clay content [%]")
ylabel!("C₃ [m]")
ylims!(0, 2.5)