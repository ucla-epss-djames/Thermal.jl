using DifferentialEquations
using Plots

const σ = 5.670374419e-8

L_int(R, Teff, Teq) = 4π*R^2 * σ * (Teff^4 - Teq^4)
T_eff(T1, B) = (B / T1)^-1.244


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
    T = T_eff(T1, B)

    p = (M, R, Cp, Teq)
    tspan = (0, 10*GyrToSec)

    function dTdt(T, p, t)

        M = p[1]
        R = p[2]
        Cp = p[3]
        Teq = p[4]

        Lf = 4π*R^2 * σ * (T^4 - Teq^4)

        dT = -Lf / (Cp * M)

        return dT

    end


    prob = ODEProblem(dTdt, T, tspan, p)
    sol = solve(prob, reltol=1e-8, abstol=1e-8)

    T = sol[1,:]
    t = sol.t
    t /= GyrToSec

    plot(t, T, xaxis=("Time (Gyr)"), yaxis=("T_eff (K)", (0,140)))

end


main()
