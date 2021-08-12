module Thermal

const σ = 5.670374419e-8

L_int(R::Real, Teff::Real, Teq::Real) = 4π*R^2 * σ * (Teff^4 - Teq^4)
T_eff(T1::Real, B::Real) = (B / T1)^-1.244

function dTdt(dT, T, p, t)
    M = p[1]
    R = p[2]
    Cp = p[3]
    Teq = p[4]

    Lf = L_int(R, T, Teq)

    dT = -Lf / (Cp * M)

end

end # module
