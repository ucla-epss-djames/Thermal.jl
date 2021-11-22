module Thermal

using PhysicalConstants.CODATA2018: σ

export L_int, T_eff

L_int(R::Real, Teff::Real, Teq::Real) = 4π*R^2 * σ.val * (Teff^4 - Teq^4)
T_eff(T1::Real, B::Real) = (B / T1)^-1.244

end # module
