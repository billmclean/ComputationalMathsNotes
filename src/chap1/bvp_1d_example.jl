using PyPlot

L = 1.0
γ0 = 0.5
γL = -0.5
x = range(0, stop=L, length=201)
u(x, c) = ( ( L - x ) * γ0 + x * γL ) / L + (c/2) * x * ( L - x )

figure(1)
cvals = collect(6:-3:-6)
plot(x, u.(x,cvals'))
grid(true)
xlabel(L"x")
legend([latexstring("c=$(cvals[k])") for k = 1:length(cvals)])
savefig("bvp_1d_example.pdf")
