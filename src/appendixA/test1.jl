using GeneralTridiagonal
using LinearAlgebra: Tridiagonal

const T = Rational{Int64}

β = T[-1, 3, 0, 5]
α = T[0, 8, 2, 1, 8]
γ = T[1, 4, 0, -3]
δ = Array{T}(undef, 3)
p = Array{Bool}(undef, 4)

A = Tridiagonal(β, α, γ)
x = T[-2, 5, 8, 0, -3]
b = A*x

factorize!(β, α, γ, δ, p)
println("Row swaps: ", p)

solve!(b, β, α, γ, δ, p)
println("   Exact solution: x = ", x)
println("Computed solution: x = ", b)
