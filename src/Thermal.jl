module Thermal

using PhysicalConstants.CODATA2014: σ

export planet_eta
export temp_effective, temp_adiabat, temp_melting
export cc_slope, temp_distribution, visc_height, thermal_diff
export lumin_internal, lumin_core, layer_density

"""
    planet_eta(η0::Real, A::Real, T_m::Real, T::Real)

Calculates viscosity. Refer to Stixrude et al. 2021 (eq 15).

# Arguments
- `η0::Real`  - reference viscosity
- `A::Real`   - viscosity exponent
- `T_m::Real` - melting temperature
- `T::Real`   - temperature
"""
planet_eta(η0::Real, A::Real, T_m::Real, T::Real) = η0 * exp(A * (T_m/T - 1))

"""
    temp_effective(T1::Real, B::Real)

Calcultes the effective temperature where it is the temperature the planet would
have in absence of solar luminosity. Refer to Stixrude et al. 2021 (eq 13) or
Guillot et al. 1995.

# Arguments
- `T1::Real` - temperature
- `B::Real`  - constant of normalization
"""
function temp_effective(T1::Real, B::Real)
    (T1 / B)^(1/1.244)
end

"""
    temp_adiabat(P::Real, Tx::Real, Px::Real, ∇::Real)

Calculates the adiabatic temperature. Refer to Stixrude et al. 2021 (eq 12).

# Arguments
- `P::Real`  - pressure
- `Tx::Real` - temperature
- `Px::Real` - reference pressure
- `∇::Real`  - adiabatic gradient
"""
function temp_adiabat(P::Real, Tx::Real, Px::Real, ∇::Real)
    Tx * (P / Px)^∇
end

"""
    temp_melting(P::Real, P0::Real, T0::Real, a::Real, b::Real)

Calculates the melting temperature. Refer to Stixrude et al. 2021 (eq 11).

# Arguments
- `P::Real`   - pressure
- `P0::Real`  - phase transition reference pressure
- `T0::Real`  - phase transition reference temperature
- `a::Real`   - Simon pressure
- `b::Real`   - Simon exponent
"""
function temp_melting(P::Real, P0::Real, T0::Real, a::Real, b::Real)
    T0 * (1 + (P - P0) / a)^b
end

"""
    cc_slope(P::Real, P0::Real, T0::Real, a::Real, b::Real)

Calculates the Clausius-Clapeyron slope of the phase boundary, Γ, where
Γ = d(T_m)/dP.

# Arguments
- `P::Real`   - pressure
- `P0::Real`  - phase transition reference pressure
- `T0::Real`  - phase transition reference temperature
- `a::Real`   - Simon pressure
- `b::Real`   - Simon exponent
"""
function cc_slope(P::Real, P0::Real, T0::Real, a::Real, b::Real)
    T0 * b / a * (1 + (P - P0) / a)^(b - 1)
end

"""
    temp_distribution(P_c::Real, T_c::Real, ∇::Real)

Calculates the slope, Γ_a, of the planetary temperature distribution at `T_c` and
`P_c`.

# Arguments
- `P_c::Real` - pressure at top of the core
- `T_c::Real` - temperature at top of the core
- `∇::Real`   - adiabatic gradient
"""
function temp_distribution(P_c::Real, T_c::Real, ∇::Real)
    (T_c / P_c) * ∇
end

"""
    visc_height(T_c::Real, ρ_c::Real, g_c::Real, A::Real, Γ::Real, Γ_a::Real)

Calculates the viscous scale height, h.

# Arguments
- `T_c::Real` - temperature at top of the core
- `ρ_c::Real` - density at top of the core
- `g_c::Real` - gravity at top of the core
- `A::Real`   - viscosity exponent
- `Γ::Real`   - Clausius-Clapeyron slope of phase boundary
- `Γ_a::Real` - slope of planetary temp distribution
"""
function visc_height(T_c::Real, ρ_c::Real, g_c::Real, A::Real, Γ::Real, Γ_a::Real)
    T_c / (ρ_c * g_c * A * (Γ - Γ_a))
end

"""
    thermal_diff(k::Real, ρ::Real, C_p::Real)

Calculates the thermal diffusivity.

# Arguments
- `k::Real`   - thermal conductivity
- `ρ::Real`   - density
- `C_p::Real` - specific heat
"""
function thermal_diff(k::Real, ρ::Real, C_p::Real)
    k / (ρ * C_p)
end

"""
    lumin_internal(R::Real, T::Real, T_eq::Real)

Calculates the total luminosity of the interior.

# Arguments
- `R::Real`    - radius
- `T_ef::Real` - effective temperature
- `T_eq::Real` - radiative equilibrium temperature
"""
function lumin_internal(R::Real, T_ef::Real, T_eq::Real)
    4π*R^2 * σ.val * (T_ef^4 - T_eq^4)
end

"""
    lumin_core(T1::Real, Ti::Real, c::Real, P_c::Real, T_c::Real, ρ_c::Real,
               g_c::Real, P1::Real, plnt)

Calculates the luminosty due to the core. Refer to Stixrude et al. 2021 (eq 9-16).

# Arguments
- `T1::Real`  - temperature of envelope
- `Ti::Real`  - temperature at top of thermal boundary
- `c::Real`   - radius of core
- `P_c::Real` - pressure at top of the core
- `T_c::Real` - temperature at top of the core
- `ρ_c::Real` - density at top of the core
- `g_c::Real` - gravity at top of the core
- `P1::Real`  - reference pressure
- `plnt`      - Planet object from `Planets.jl`
"""
function lumin_core(T1::Real, Ti::Real, c::Real, P_c::Real, T_c::Real, ρ_c::Real,
                    g_c::Real, P1::Real, plnt)
    if(Ti < 0)
        P_c = plnt.P0 - plnt.a + 1e-6
        Ti = temp_adiabat(P_c, T1, P1, plnt.∇)
    end

    P_ratio = (P_c / P1)^plnt.∇

    ΔT = Ti - T_c
    ΔP = P_c - plnt.P0

    Γ = cc_slope(P_c, plnt.P0, plnt.T0, plnt.a, plnt.b)
    Γ_a = temp_distribution(P_c, T_c, plnt.∇)
    ΔΓ = Γ - Γ_a

    h = visc_height(T_c, ρ_c, g_c, plnt.A, Γ, Γ_a)
    if(h > c) h = c end

    T_m = temp_melting(P_c, plnt.P0, plnt.T0, plnt.a, plnt.b)

    η = planet_eta(plnt.η0, plnt.A, T_m, T_c)

    Ra = ρ_c * plnt.α * abs(ΔT) * g_c * abs(h)^3 / (plnt.κ * η)

    K = ρ_c * plnt.α * g_c / (plnt.Ra * η * plnt.κ)
    L_c = plnt.k * K^(1/3) * ((h / c)^(4))^(1/3) * (ΔT^4)^(1/3) * 4*π*c^2

    # this is value next to ∂c/∂t term in the derivation
    K = - 1 / (ρ_c * g_c * ΔΓ) * P_ratio

    return (L_c=L_c, K=K, ΔT=ΔT, Γ=Γ_a/ΔΓ, Ra=Ra)
end

"""
    layer_density(r::Real, P::Real, Px::Real, ∇::Real, ρ::Real)

Function integrated for luminosity. Refer to Stixrude et al. 2021 (eq 8-9) or
Fortney & Nettelmann 2010.

# Arguments
- `r::Real`  - radius
- `P::Real`  - pressure
- `Px::Real` - reference pressure
- `∇::Real`  - adiabatic gradient
- `ρ::Real`  - density
"""
function layer_density(r::Real, P::Real, Px::Real, ∇::Real, ρ::Real)
    (P/Px)^∇ * ρ * r^2
end

end # module
