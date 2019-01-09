using Printf
using PyPlot

function inverse_const(r)
    B = [ 1 / (i+j-1) for i = 1:r+1, j=1:r+1 ]
    D = zeros(r+1, r+1)
    for j = 2:r+1
        for i = 2:r+1
            D[i,j] = (i-1)*(j-1) / (i+j-3)
        end
    end
    F = eigen(D, B)
    λmax = F.values[r+1]
    amax = F.vectors[:,r+1]
    return λmax, amax
end

figure(1)
ξ = range(0, stop=1, length=201)
@printf("\n%4s %8s\n\n", "r", "Cr")
for r = 0:3
    Cr, ar = inverse_const(r)
    pr = Poly(ar, :ξ)
    plot(ξ, polyval(pr, ξ))
    @printf("%4d& %8.2f\\\\\n", r, Cr)
end
grid(true)
