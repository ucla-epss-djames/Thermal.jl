using DifferentialEquations
using Plots

const σ = 5.670374419e-8

L_f(R, Teff, Teq) = 4π*R^2 * σ * (Teff^4 - Teq^4)
T_eff(T1, B) = (T1 / B)^(1/1.244)
T_1(Teff, B) = B * (Teff + 0*im)^1.244


function fluid()

    GyrToSec = 1e9*365*24*60*60
    GPaToBar = 1e4
    ∇ = 0.2585

    Me = 5.972e24
    Re = 6371e3
    T = 59.1
    T1 = 200
    P0 = 40

    M = 6.55*Me
    R = 2.68*Re
    Cp = 5000
    Teq = 58.1

    B = 0.47529
    Teff = T_eff(T1, B)

    p = (M, R, Cp, Teq, B)
    tspan = (0., 10.0*GyrToSec)

    function dTdt(Teff, p, t)

        M = p[1]
        R = p[2]
        Cp = p[3]
        Teq = p[4]
        B = p[5]

        # Teff = T_eff(T, B)
        # T = T_1(T,B) * (P0 * GPaToBar)^∇
        # T = (T / B * 1 / (P0 * GPaToBar)^∇)^-1.244

        Lf = 4π*R^2 * σ * (Teff^4 - Teq^4)

        dT = -Lf / (Cp * M)


        return dT

    end

    prob = ODEProblem(dTdt, T1, tspan, p)
    sol = solve(prob, reltol=1e-8, abstol=1e-8)

    T = sol[1,:]
    t = sol.t
    t /= GyrToSec

    display(T)

    plot(t, T, color="red", xaxis=("Time (Gyr)", 0:2:10), yaxis=("T_eff (K)", (40,140)),
         size=(400,400))

    savefig("plot.png")

end

function fluidcore()

end

fluid()
