### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 0066388e-a82c-11f0-0374-3b44706faa2f
begin
	using DrWatson
	@quickactivate "DifferentiableEvaporation"
	using Symbolics, Nemo
end

# ╔═╡ d306b779-3ff3-495e-86b7-b6c175cf718b
md"""
The goal of the current notebook is to extend upon the derivations of a series network from:

- "Evaporation from sparse crops-an energy combination theory" by [Shuttleworth and Wallace (1985)](https://doi.org/10.1002/qj.49711146910) (S&W85)
- "Evaporation from Heterogeneous and Sparse Canopies: On the Formulations Related to Multi-Source Representations" by [Lhomme et al. (2012)](https://doi.org/10.1007/s10546-012-9713-x) (Lh12)

Start by defining the different variables
"""

# ╔═╡ a61c6a67-213b-4867-bfa4-42ab0fbe7bb6
@variables begin
    As # Available energy at the surface (soil) (W/m²)
    Ac # Available energy at the canopy (W/m²)
    ρ # Density of air (kg/m³)
    cp # Specific heat of air (J/kg/K)
    Vpd_a # Vapor pressure deficit (Pa) at observation height
    γ # Psychrometric constant (Pa/K)
    Δ # Slope of the saturation vapor pressure curve (Pa/K) 
    r_aa # Aerodynamic resistance (s/m) between source height and observation height
    r_ac # Aerodynamic resistance (s/m) between canopy and source height
    r_as # Aerodynamic resistance (s/m) between soil and source height
    r_ss # Surface resistance (s/m) of soil
    r_sc # Surface resistance (s/m) of canopy
    λE # Latent heat flux from the complete canopy (W/m²)
    # Ra
    # Rs
    # Rc
    f_wet # Fraction of canopy that is wet (-)
end

# ╔═╡ 68e671ec-c356-4871-b2f8-51c9a3e3b949
md"""
## Reference case from Lhomme et al. (2012)

Let's symbolically solve the equations governing the standard case of a two source model, as depicted in Figure 1 of Lh12.

"""

# ╔═╡ 578770e6-5440-4458-9bae-4333c0cf1f93
begin
	A = As + Ac # Total available energy
	Vpd_0 = Vpd_a + (Δ * A - (Δ + γ) * λE) * r_aa / (ρ * cp) # Vapor pressure deficit at source height, (7) of Lh12
	λE_s = (Δ * As + ρ * cp * Vpd_0 / r_as) / (Δ + γ * (1 + r_ss / r_as)) # Latent heat flux from the soil
	λE_c = (Δ * Ac + ρ * cp * Vpd_0 / r_ac) / (Δ + γ * (1 + r_sc / r_ac)) # Latent heat flux from the canopy
	equation = λE - (λE_s + λE_c) ~ 0 # Latent heat flux balance, (1) of Lh12
end

# ╔═╡ 5112f5f3-7b44-41fb-b046-640f21eaca23
begin
	sol = symbolic_solve(equation, λE)
	sol[1];
end

# ╔═╡ c5db62ac-922c-4749-a363-214ce64de67d
md"""
Check if the solution matches with the solution as given by equation 16 of Lh12. Note that f (folliage) replaced by c (canopy). First this reference solution is implemented below.
"""

# ╔═╡ f71769f1-0ec8-4ddf-ab4d-45459c0ff5bc
begin
	Rc = r_sc + (1 + Δ / γ) * r_ac # (20) of Lh12
	Rs = r_ss + (1 + Δ / γ) * r_as # (21) of Lh12
	Ra = (1 + Δ / γ) * r_aa # (22) of Lh12
	Ps = (r_aa * Rc) / (Rc * Rs + Ra * Rc + Rs * Ra) # (18) of Lh12
	Pc = (r_aa * Rs) / (Rc * Rs + Ra * Rc + Rs * Ra) # (19) of Lh12
	λEp = (Δ * A + ρ * cp * Vpd_a / r_aa) / (Δ + γ) # (17) of Lh12
	λE_lhomme = ((Δ + γ) / γ) * (Pc + Ps) * λEp +
	            (Δ / (γ * r_aa)) * (Pc * Ac * r_ac + Ps * As * r_as); # (16) of Lh12
end

# ╔═╡ b4b0df3d-35f8-4a23-ba38-a549d57061f9
md"""
Now symbolically check the equality between the two solutions by subtracting the two equations from each other after simplifying them. 
"""

# ╔═╡ 6812d801-924f-4883-bdb8-a347b10b90b8
begin
	λE_computer_simple = simplify(sol[1])
	λE_lhomme_simple = simplify(λE_lhomme)
	print(isequal(simplify(λE_computer_simple - λE_lhomme_simple), 0))
end

# ╔═╡ 3d3fad59-a0fb-47d8-91fe-b19e3b2efb39
md"""
## Own derivation: adding interception 

Idea is to extend the fluxes so that from a leave, there are now 2 fluxes in parallel:
1. Interception flux from wet fraction of the leaf ($f_{wet}$) with a zero surface resistance
2. Evaporation flux from the dry fraction of the leaf ($1 - f_{wet}$)
"""

# ╔═╡ e507c1c5-78f8-4417-bb8d-42ec07d2fb39
begin
	λE_t = (1 - f_wet) * (Δ * Ac + ρ * cp * Vpd_0 / r_ac) / (Δ + γ * (1 + r_sc / r_ac))
	λE_i = f_wet * (Δ * Ac + ρ * cp * Vpd_0 / r_ac) / (Δ + γ)
	equation_interception = λE - (λE_t + λE_i + λE_s) ~ 0
end

# ╔═╡ 172b5b36-4bc4-483f-87c7-096c4437c693
sol_interception = symbolic_solve(equation_interception, λE);

# ╔═╡ a0552080-bab6-49f3-9bcd-eb44d53eaa7a
md"""
Below the solution from the manual derivation is given.
"""

# ╔═╡ 2ead5dc5-a3c7-4924-b387-f95f0efa5ea9
begin
	Ri = (1 + Δ / γ) * r_ac # In analogy to (20) - (22) of Lh12
	X = Rc * Ri + (1 - f_wet) * Rs * Ri + f_wet * Rs * Rc
	DE = Rs * Rc * Ri + Ra * X
	λE_check = (X * r_aa / DE) * ((Δ + γ) / γ) * λEp + (Δ / γ) * ((Rc * Ri * As * r_as + r_ac * Ac * ((1 - f_wet) * Rs * Ri + f_wet * Rc * Rs)) / DE)
end

# ╔═╡ cf75cc24-bf75-4540-b650-25a89b9eb320
begin
	λE_check_simple = simplify(λE_check)
	λE_interception_simple = simplify(sol_interception[1])
	print(isequal(simplify(λE_check_simple - λE_interception_simple), 0))
end

# ╔═╡ 39e86399-ffa4-4ecc-aa09-033f13ba6785
md"""
The same derivation can also be written in a form closer to the original equations of Lh12 (see (16) and (33)). Also this form will be checked against the computer solution.
"""

# ╔═╡ 52d4b78d-4eb1-455c-8b79-8eaf7b0fc14b
begin
	DE_ = Rs * Rc * Ri + Ra * (Rc * Ri + (1 - f_wet) * Rs * Ri + f_wet * Rs * Rc)
	Pc_ = simplify(r_aa * (1 - f_wet) * Ri * Rs / DE_)
	Pi_ = simplify(r_aa * f_wet * Rc * Rs / DE_)
	Ps_ = simplify(r_aa * Rc * Ri / DE_)
	λE_check_2 = simplify(((Δ + γ) / γ) * (Pc_ + Pi_ + Ps_) * λEp) +
	             simplify((Δ / (γ * r_aa)) * (Pc_ * Ac * r_ac + Pi_ * Ac * r_ac + Ps_ * As * r_as))
	λE_check_2_simple = simplify(λE_check_2)
	print(isequal(simplify(λE_check_2_simple - λE_interception_simple), 0))
end

# ╔═╡ dd9f3b59-fed2-4a33-8be6-bdc36d61b488
md"""
There is also a different approach to solve the problem by hand, where the $R_c$ value as defined in (20) of Lh12 is adapted for the current case (denoted as $R_c^*$). Once this is defined, we can use the same equations as in the reference case (Section 2.2 of Lh12).
"""

# ╔═╡ 904921b6-3dbf-42b1-b455-ea0c1aa2c2b6
begin
	Rcstar = ((1 - f_wet) / Rc + f_wet / Ri)^-1
	Psstar = simplify((r_aa * Rcstar) / (Rcstar * Rs + Ra * Rcstar + Rs * Ra)) # Analogous to (18) of Lh12
	Pcstar = simplify((r_aa * Rs) / (Rcstar * Rs + Ra * Rcstar + Rs * Ra)) # Analogous to (19) of Lh12
	λE_check_alt = simplify(((Δ + γ) / γ) * (Pcstar + Psstar) * λEp) +
	               simplify((Δ / (γ * r_aa)) * (Pcstar * Ac * r_ac + Psstar * As * r_as)) # Analogous to (16) of Lh12
	λE_check_alt_simple = simplify(λE_check_alt)
	print(isequal(simplify(λE_interception_simple - λE_check_alt_simple), 0))
end

# ╔═╡ c916d907-7eb2-4082-9592-d6720c0a7250
md"""
## Add derivation for $\beta$ equivalance
"""

# ╔═╡ f0e696ed-eb14-43f2-958a-995502c3cc9e
@variables β, esatTs, e0


# ╔═╡ fac0dc07-d515-4ae6-8971-ee899b9a3edf
begin
	λE_sp = (Δ * As + ρ * cp * Vpd_0 / r_as) / (Δ + γ)
	equation_β = β* λE_sp - λE_s ~ 0 # Latent heat flux balance, (1) of Lh
end

# ╔═╡ 1e66998d-1187-4cc2-ad76-a6363fefc11c
β_pm = symbolic_solve(equation_β, β)[1]

# ╔═╡ 84fbd0ae-e2c2-4773-8c11-5fc9ca6cf009
λE_sbulk = ((ρ * cp)/ γ) * ((esatTs - e0) / (r_as + r_ss))

# ╔═╡ 51c3fb7f-94a3-4054-a413-8e8cec863d01
λE_sbulk_beta = β * ((ρ * cp)/ γ) * ((esatTs - e0) / (r_as))


# ╔═╡ 874b9e40-3a07-42e8-8989-4470d039e2cd
β_bulk = symbolic_solve(λE_sbulk_beta - λE_sbulk, β)[1]

# ╔═╡ Cell order:
# ╠═0066388e-a82c-11f0-0374-3b44706faa2f
# ╟─d306b779-3ff3-495e-86b7-b6c175cf718b
# ╠═a61c6a67-213b-4867-bfa4-42ab0fbe7bb6
# ╟─68e671ec-c356-4871-b2f8-51c9a3e3b949
# ╠═578770e6-5440-4458-9bae-4333c0cf1f93
# ╠═5112f5f3-7b44-41fb-b046-640f21eaca23
# ╟─c5db62ac-922c-4749-a363-214ce64de67d
# ╠═f71769f1-0ec8-4ddf-ab4d-45459c0ff5bc
# ╟─b4b0df3d-35f8-4a23-ba38-a549d57061f9
# ╠═6812d801-924f-4883-bdb8-a347b10b90b8
# ╟─3d3fad59-a0fb-47d8-91fe-b19e3b2efb39
# ╠═e507c1c5-78f8-4417-bb8d-42ec07d2fb39
# ╠═172b5b36-4bc4-483f-87c7-096c4437c693
# ╟─a0552080-bab6-49f3-9bcd-eb44d53eaa7a
# ╠═2ead5dc5-a3c7-4924-b387-f95f0efa5ea9
# ╠═cf75cc24-bf75-4540-b650-25a89b9eb320
# ╟─39e86399-ffa4-4ecc-aa09-033f13ba6785
# ╠═52d4b78d-4eb1-455c-8b79-8eaf7b0fc14b
# ╟─dd9f3b59-fed2-4a33-8be6-bdc36d61b488
# ╠═904921b6-3dbf-42b1-b455-ea0c1aa2c2b6
# ╟─c916d907-7eb2-4082-9592-d6720c0a7250
# ╠═f0e696ed-eb14-43f2-958a-995502c3cc9e
# ╠═fac0dc07-d515-4ae6-8971-ee899b9a3edf
# ╠═1e66998d-1187-4cc2-ad76-a6363fefc11c
# ╠═84fbd0ae-e2c2-4773-8c11-5fc9ca6cf009
# ╠═51c3fb7f-94a3-4054-a413-8e8cec863d01
# ╠═874b9e40-3a07-42e8-8989-4470d039e2cd
