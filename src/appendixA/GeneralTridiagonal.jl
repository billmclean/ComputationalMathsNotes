module GeneralTridiagonal

export factorize!, solve!

"""
    factorize!(β, α, γ, δ, p)

Computes the LU factorization with partial pivoting for an `n×n` tridiagonal 
matrix `A` with subdiagonal `β`, main diagonal `α` and superdiagonal `γ`.  Thus,
the vector `α` has length `n`, whereas `β` and `γ` each have length `n-1`.
For example, if `n=5` then the matrix is given by

          α_1  γ_1
          β_1  α_2  γ_2
    A =        β_2  α_3  γ_3
                    β_3  α_4  γ_4
                         β_4  α_5

The vector `δ` must have length `n-2` to hold the second superdiagonal of 
the upper triangular factor, arising from the fill-in caused by pivoting.
The boolean vector `p` has length `n-1` and records the sequence of pivot
choices: `p[j]` is `true` if rows `j` and `j+1` were swapped at the `j`th 
step.

The factorization is computed in place, so `β` is overwritten with the
nonzero subdiagonal entries of `L`, and `α`, `γ` and `δ` hold the nonzero
entries of `U`.
"""
function factorize!(β::Vector{T}, α::Vector{T}, γ::Vector{T}, δ::Vector{T}, 
                    p::Vector{Bool}) where T <: Number
    n = length(α)
    if ! ( length(β) == n-1 == length(γ) && length(δ) == n-2 )
        throw(BoundsError("Array sizes to not match"))
    end
    for j = 1:n-2
        δ[j] = zero(T)
        if abs(α[j]) < abs(β[j])
            p[j] = true
            α[j], β[j] = β[j], α[j]
            γ[j], α[j+1] = α[j+1], γ[j]
            δ[j], γ[j+1] = γ[j+1], δ[j]
        else
            p[j] = false
        end
        β[j] = β[j] / α[j]
        α[j+1] -= β[j] * γ[j]
        if p[j]
            γ[j+1] = -β[j] * δ[j]
        end
    end
    if abs(α[n-1]) < abs(β[n-1])
        p[n-1] = true
        α[n-1], β[n-1] = β[n-1], α[n-1]
        γ[n-1], α[n]   = α[n], γ[n-1]
    else
        p[n-1] = false
    end
    β[n-1] = β[n-1] / α[n-1]
    α[n] -= β[n-1] * γ[n-1]
end

"""
    solve!(b, β, α, γ, δ, p)

Solves the tridiagonal linear system `Ax = b` using the LU factorization 
computed by `factorize`.
"""
function solve!(b::Vector{T}, β::Vector{T}, α::Vector{T}, γ::Vector{T}, 
                δ::Vector{T}, p::Vector{Bool}) where T <: Number
    n = length(α)
    if ! ( length(β) == n-1 == length(γ) && length(δ) == n-2 )
        throw(BoundsError("Array sizes to not match"))
    end
    for j = 1:n-1
        if p[j]
            b[j], b[j+1] = b[j+1], b[j]
        end
        b[j+1] -= β[j] * b[j]
    end
    b[n] /= α[n]
    b[n-1] = ( b[n-1] - γ[n-1] * b[n] ) / α[n-1]
    for j = n-2:-1:1
        b[j] = ( b[j] - γ[j] * b[j+1] - δ[j] * b[j+2] ) / α[j]
    end
end

end # module
