using DifferentialEquations
using Polynomials
using QuadGK
using Plots

const σ = 5.670374419e-8
const GPaToBar = 10000

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

    arr = [4424.50390625, 0.0, 0.0, 0.0, -49167.8, 71977.3,
           -27234.0, 0.0, 0.0]
    ρ = Polynomial(arr)

    B = 0.47529
    Teff = T_eff(T1, B)

    mass(r) = ρ(r/R) * 4π * r^2

    function f(r)

        return (ρ(r/R) * 8.314 * T1 * GPaToBar)^∇ * mass(r) * Cp

    end


    # M = quadgk(mass, 0, R)
    # In = (500*GPaToBar)^∇ * Cp * M[1]

    In = quadgk(f, 0, R)

    p = (In[1], R, Teq)
    tspan = (0., 10.0*GyrToSec)

    function dTdt(Teff, p, t)

        In = p[1]
        R = p[2]
        Teq = p[3]

        Lf = 4π*R^2 * σ * (Teff^4 - Teq^4)

        dT = -Lf / In


        return dT

    end

    prob = ODEProblem(dTdt, T1, tspan, p)
    sol = solve(prob, reltol=1e-8, abstol=1e-8)

    T = sol[1,:]
    t = sol.t
    t /= GyrToSec

    plot(t, T, color="red", xaxis=("Time (Gyr)", 0:2:10), yaxis=("T_eff (K)", (40,140)),
         size=(400,400), label="No core", dpi=300)
    hline!([59.1], color="black", label="59.1 K")

    # savefig("plot.png")

end

function helled_model()

    Re = 6371e3
    R = 2.68*Re
    # 
    arr = [4424.50390625, 0.0, 0.0, 0.0, -49167.8, 71977.3,
           -27234.0, 0.0, 0.0]

    # arr = [4039.37, 4039.5087878704071, 0.0, 20.701, 3.78416, -38675.2,
    #        53208.8, -18597.6, 8231.99, 0.2]

    # arr = arr[end:-1:1]

    P = Polynomial(arr)

    print(P)

    x = range(0, 1, length=100)

    display(P.(x))

    plot(x, P.(x))

end

fluid()
