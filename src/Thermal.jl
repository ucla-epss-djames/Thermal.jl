module Thermal

using PhysicalConstants.CODATA2018: σ

export BarToPa, GyrToSec
export L_int, T_eff

const BarToPa = 1e5
const GyrToSec = 1e9*365*24*60*60

L_int(R::Real, Teff::Real, Teq::Real) = 4π*R^2 * σ.val * (Teff^4 - Teq^4)
T_eff(T1::Real, B::Real) = (B / T1)^-1.244

end # module
