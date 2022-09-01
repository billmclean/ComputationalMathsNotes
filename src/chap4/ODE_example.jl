using OffsetArrays
import LinearAlgebra: I, lu
using PyPlot

A = [ 0.5  -20.0
      0.0   20.0 ]
u0 = [ 0.0
       1.0 ]

max_t = 2.0
N = 25

function exact_soln(t::OffsetVector{T}, 
	            u0::Vector{T}) where T <: AbstractFloat
    N = length(t) - 1
    Δt =  t[N] / N
    u = OffsetMatrix{T}(undef, 2, 0:N)
    u[:,0] .= u0
    for n = 1:N
	u[1,n] = (40/39) * ( exp(-t[n]/2) - exp(-20t[n]) )
	u[2,n] = exp(-20t[n])
    end
    return u
end

function explicit_Euler(t::OffsetVector{T}, 
	                u0::Vector{T}) where T <: AbstractFloat
    N =length(t) - 1
    Δt = t[N] / N
    U = OffsetMatrix{T}(undef, 2, 0:N)
    U[:,0] .= u0
    B = I - Δt * A 
    for n = 0:N-1
	U[:,n+1] = B * U[:,n]
    end
    return U
end

function implicit_Euler(t::OffsetVector{T}, 
	                u0::Vector{T}) where T <: AbstractFloat
    N =length(t) - 1
    U = OffsetMatrix{T}(undef, 2, 0:N)
    Δt = t[N] / N
    U[:,0] .= u0
    B = I + Δt * A 
    F = lu(B)
    for n = 1:N
	U[:,n] = F \ U[:,n-1]
    end
    return U
end

function Crank_Nicolson(t::OffsetVector{T}, 
	                u0::Vector{T}) where T <: AbstractFloat
    N =length(t) - 1
    U = OffsetMatrix{T}(undef, 2, 0:N)
    Δt = t[N] / N
    U[:,0] .= u0
    B₊ = I + (Δt/2) * A 
    B₋ = I - (Δt/2) * A
    F = lu(B₊)
    for n = 1:N
	U[:,n] = F \ (B₋ * U[:,n-1] )
    end
    return U
end

t = OffsetVector(range(0, max_t, length=N+1), 0:N)
N_fine = 201
t_fine = OffsetVector(range(0, max_t, length=N_fine+1), 0:N_fine)
u = exact_soln(t_fine, u0)
UEE = explicit_Euler(t, u0)
UIE = implicit_Euler(t, u0)
UCN = Crank_Nicolson(t, u0)

figure(1)
plot(t_fine, u[1,:], ":", t, UEE[1,:], "--", t, UIE[1,:], t, UCN[1,:], )
grid(true)
xlabel(L"$t$")
legend(("Exact solution", "Explicit Euler", "Implicit Euler", "Crank-Nicolson"))

figure(2)
plot(t_fine, u[2,:], ":", t, UEE[2,:], "--", t, UIE[2,:], t, UCN[2,:])
axis([0, max_t, -0.5, 1.0])
grid(true)
xlabel(L"$t$")
legend(("Exact solution", "Explicit Euler", "Implicit Euler", "Crank-Nicolson"))
