module FEM_P2

using SparseArrays

@enum BdryCond Dirichlet Neumann

export Dirichlet, Neumann, Mesh 
export assemble_matrix, assemble_vector, Neumann_bc_vector
export deriv_times_deriv!, func_times_func!
export evaluate_from_ndvals

struct Mesh
    node    :: Vector{Float64}
    elm     :: Matrix{Int64}
    left    :: BdryCond
    right   :: BdryCond
    ndof    :: Int64
    leftnd  :: Int64
    rightnd :: Int64
end

# Element stiffness matrix = MAT1 / (3hm)
const MAT1 = [ 7.0  -8.0   1.0
              -8.0  16.0  -8.0
               1.0  -8.0   7.0 ] 

# Element mass matrix = (hm/30) * MAT2
const MAT2 = [  4.0   2.0  -1.0
                2.0  16.0   2.0
               -1.0   2.0   4.0 ]

Psi1(ξ) = 2 * ( ξ - 1/2 ) * ( ξ - 1 )
Psi2(ξ) = 4 * ξ * ( 1 - ξ )
Psi3(ξ) = 2 * ξ * ( ξ - 1/2 )
const PSI = ( Psi1, Psi2, Psi3 )

dPsi1(ξ) = 4ξ - 3
dPsi2(ξ) = 4 * ( 1 - 2ξ )
dPsi3(ξ) = 4ξ - 1
const DPSI = [ dPsi1(0.0) dPsi2(0.0) dPsi3(0.0) 
               dPsi1(0.5) dPsi2(0.5) dPsi3(0.5) 
               dPsi1(1.0) dPsi2(1.0) dPsi3(1.0) ]

const SIMPSON_WT = [ 1/6, 4/6, 1/6 ]

function Mesh(pts::Vector{Float64}, bc::Tuple{BdryCond,BdryCond})
    M = length(pts) - 1
    h = diff(pts)
    node = Vector{Float64}(undef, 2M+1)
    elm = Matrix{Int64}(undef, 3, M)
    if bc[1] == Dirichlet
        for m = 1:M
            node[2m-1] = pts[m] + h[m]/2
            node[2m]   = pts[m+1]
        end
        node[2M+1] = pts[1]
        leftnd = 2M+1
        rightnd = 2M
        elm[:,1] .= [ 2M+1, 1, 2]
        for m = 2:M
            elm[:,m] .= [ 2m-2, 2m-1, 2m ]
        end
        if bc[2] == Dirichlet
            ndof = 2M-1
        else
            ndof = 2M
        end
    else
        node[1] = pts[1]
        leftnd = 1
        rightnd = 2M+1
        for m = 1:M
            node[2m]   = pts[m] + h[m]/2
            node[2m+1] = pts[m+1]
            elm[:,m] .= [ 2m-1, 2m, 2m+1 ]
        end
        if bc[2] == Dirichlet
            ndof = 2M
        else
            ndof = 2M+1
        end
    end
    return Mesh(node, elm, bc[1], bc[2], ndof, leftnd, rightnd) 
end

function assemble_matrix(mesh::Mesh, elm_mat!::Function, 
                         coef::Union{Float64,Function})
    node, elm, ndof = mesh.node, mesh.elm, mesh.ndof
    M = size(elm, 2)
    Am = zeros(3, 3)
    nodem = zeros(3)
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for m = 1:M
        for p = 1:3
            nodem[p] = node[elm[p,m]]
        end
        elm_mat!(Am, nodem, coef)
        for p = 1:3
            if elm[p,m] > ndof
                continue
            end
            for q = 1:3
                push!(I, elm[p,m])
                push!(J, elm[q,m])
                push!(V, Am[p,q])
            end
        end
    end
    return sparse(I, J, V, ndof, 2M+1)
end

function assemble_vector(mesh::Mesh, f::Function)
    node, elm, ndof = mesh.node, mesh.elm, mesh.ndof
    M = size(mesh.elm, 2)
    F = zeros(ndof)
    Fm = zeros(3)
    nodem = zeros(3)
    I = Int64[]
    V = Float64[]
    for m = 1:M
        for p = 1:3
            nodem[p] = node[elm[p,m]]
        end
        # Use Simpson's rule
        hm = nodem[3] - nodem[1]
        for p = 1:3
            Fm[p] = hm * SIMPSON_WT[p] * f(nodem[p])
        end
        for p = 1:3
            if elm[p,m] > ndof
                continue
            end
            F[elm[p,m]] += Fm[p]
        end
    end
    return F
end

function Neumann_bc_vector(mesh::Mesh, γ::Vector{Float64})
    left, right, ndof = mesh.left, mesh.right, mesh.ndof
    leftnd, rightnd = mesh.leftnd, mesh.rightnd
    count = 0
    if left == Neumann
        count += 1
    end
    if right == Neumann
        count += 1
    end
    if length(γ) != count
        ArgumentError("Wrong number of Neumann boundary values")
    end
    I = zeros(Int64, count)
    V = zeros(Float64, count)
    k = 0
    if left == Neumann
        k += 1
        I[k] = leftnd
        V[k] = γ[k]
    end
    if right == Neumann
        k += 1
        I[k] = rightnd
        V[k] = γ[k]
    end
    return sparsevec(I, V, ndof)
end

function deriv_times_deriv!(A::Matrix{Float64}, node::Vector{Float64},
                            coef::Float64)
    hm = node[3] - node[1]
    A .= MAT1
    A .*= coef/(3hm)
end

function deriv_times_deriv!(A::Matrix{Float64}, node::Vector{Float64},
                            coef::Function)
    hm = node[3] - node[1]
    coefs = coef.(node)
    for q = 1:3
        for p = 1:q-1
            A[p,q] = A[q,p]
        end
        for p = q:3
            s = 0.0
            for j = 1:3
                s += SIMPSON_WT[j] * coefs[j] * DPSI[j,q] * DPSI[j,p]
            end
            A[p,q] = s / hm
        end
    end
end

function func_times_func!(A::Matrix{Float64}, node::Vector{Float64},
                          coef::Float64)
    hm = node[3] - node[1]
    A .= MAT2
    A .*= hm / 30
end

function func_times_func!(A::Matrix{Float64}, node::Vector{Float64},
                          coef::Function)
    hm = node[3] - node[1]
    for p = 1:3
        A[p,p] = hm * SIMPSON_WT[p] * coef(node[p]) 
    end
end

function evaluate_from_ndvals(mesh::Mesh, Und::Vector{Float64}, 
                              r::Int64=10)
    node, elm, ndof = mesh.node, mesh.elm, mesh.ndof
    left, right = mesh.left, mesh.right
    leftnd, rightnd = mesh.leftnd, mesh.rightnd
    M = size(elm, 2)
    x = zeros(r*M+1)
    U = zeros(r*M+1)
    x[1] = node[leftnd]
    U[1] = Und[elm[1,1]]
    for m = 1:M
        lo = node[elm[1,m]]
        hi = node[elm[3,m]]
        hm = hi - lo
        step = hm / r
        for k = 1:r
            j = (m-1) * r + 1 + k
            x[j] = lo + k * step
            s = 0.0
            for p = 1:3
                ξ = k / r
                s += Und[elm[p,m]] * PSI[p](ξ)
            end
            U[j] = s
        end
    end
    return x, U
end

end # module
