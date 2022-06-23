module Thermal

using PhysicalConstants.CODATA2014: σ

using Structure: planet_eta

export Bar_to_Pa, Gyr_to_sec
export temp_effective, temp_adiabat, temp_melting
export cc_slope, phase_boundary, visc_height, thermal_diff
export lumin_internal, lumin_core, layer_density

const Bar_to_Pa = 1e5
const Gyr_to_sec = 1e9*365*24*60*60

function temp_effective(T1::Real, B::Real)
    (T1 / B)^(1/1.244)
end

function temp_adiabat(P::Real, Tx::Real, Px::Real, ∇::Real)
    Tx * (P / Px)^∇
end

function temp_melting(P::Real, P0::Real, T0::Real, a::Real, b::Real)
    T0 * (1 + (P - P0) / a)^b
end

function cc_slope(P::Real, P0::Real, T0::Real, a::Real, b::Real)
    T0 * b / a * (1 + (P - P0) / a)^(b - 1)
end

function phase_boundary(P_c::Real, T_c::Real, ∇::Real)
    (T_c / P_c) * ∇
end

function visc_height(T_c::Real, ρ_c::Real, g_c::Real, A::Real, Γ::Real, Γ_a::Real)
    T_c / (ρ_c * g_c * A * (Γ - Γ_a))
end

function thermal_diff(k::Real, ρ::Real, C_p::Real)
    k / (ρ * C_p)
end

function lumin_internal(R::Real, T::Real, T_eq::Real)
    4π*R^2 * σ.val * (T^4 - T_eq^4)
end

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
    Γ_a = phase_boundary(P_c, T_c, plnt.∇)
    ΔΓ = Γ - Γ_a

    h = visc_height(T_c, ρ_c, g_c, plnt.A, Γ, Γ_a)
    if(h > c) h = c end

    T_m = temp_melting(P_c, plnt.P0, plnt.T0, plnt.a, plnt.b)

    η = planet_eta(plnt.η0, plnt.A, T_m, T_c)

    κ = thermal_diff(plnt.k, ρ_c, plnt.C_p)

    K = ρ_c * plnt.α * g_c / (plnt.Ra * η * κ)
    L_c = plnt.k * K^(1/3) * ((h / c)^(4))^(1/3) * (ΔT^4)^(1/3) * 4*π*c^2

    # this is value next to ∂c/∂t term in the derivation
    K = - 1 / (ρ_c * g_c * ΔΓ) * P_ratio

    return (L_c=L_c, K=K, ΔT=ΔT, Γ=Γ_a/ΔΓ)
end

function layer_density(r::Real, P::Real, Px::Real, ∇::Real, ρ::Real)
    (P/Px)^∇ * ρ * r^2
end

end # module
