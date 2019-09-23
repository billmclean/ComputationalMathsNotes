module BVP1D

using LinearAlgebra
using OffsetArrays

export solve_bvp

function solve_bvp(L, f, γ0, γL, P) 
    Δx = L / P
    c = 1 / Δx^2
    dv = fill(2c, P-1)
    ev = fill(-c, P-1)
    A = SymTridiagonal(dv, ev)
    x_ = range(0; length=P+1, stop=L)
    x = OffsetArray(x_, 0:P)
    b = f.(x[1:P-1])
    b[1]   += c * γ0
    b[P-1] += c * γL
    U = OffsetArray{Float64}(undef, 0:P)
    U[0] = γ0
    U[1:P-1] = A \ b
    U[P] = γL
    return U, x
end

end # module BVP1D

