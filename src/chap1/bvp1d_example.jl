using BVP1D
using PyPlot

const L = 2.0
const γ0 = -1.0
const γL = 2.5
const m = 5.0

const A = γ0 + m
const B = ( γL - A + m*exp(-L)) / L

u(x) = A + B*x - m*exp(-x)
f(x) = m * exp(-x)

P = 6
U, x = solve_bvp(L, f, γ0, γL, P)

xx = range(0; stop=L, length=201)
plot(xx, u.(xx), "-", x[0:P], U[0:P], "o")
xticks(x, [ latexstring("\$x_$p\$") for p = 0:6 ])
legend(("Exact solution", "Finite difference approx"))
grid(true)
#xlabel(L"$x$")
savefig("bvp1d_example.pdf")
~              
