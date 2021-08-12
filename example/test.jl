include("../src/Thermal.jl")

using .Thermal

using DifferentialEquations

function main()

    GyrToSec = 1e9*365*24*60*60

    Me = 5.972e24
    Re = 6371e3
    T1 = 76

    M = 6.55*Me
    R = 2.68*Re
    Cp = 3000
    Teq = 58.1

    B = 0.47529
    T = Thermal.T_eff(T1, B)

    p = (M, R, Cp, Teq)
    tspan = (0, 0.8*GyrToSec)

    prob = ODEProblem(Thermal.dTdt, T, tspan, p, reltol=1e-10, abstol=1e-10)
    sol = solve(prob)

    display(sol)

end

main()
