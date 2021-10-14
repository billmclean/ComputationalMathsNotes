using PyPlot
using OffsetArrays

function bary_grid(n, a)
    x = OffsetArray{Float64}(undef, 0:n, 0:n)
    y = OffsetArray{Float64}(undef, 0:n, 0:n)
    pts = OffsetArray(range(0, 1, length=n+1), 0:n)
    for i = 0:n
        ξ₁ = pts[i]
        for j = 0:n-i
            ξ₂ = pts[j]
            ξ₃ = 1 - ξ₁ - ξ₂
            x[i,j] = ξ₁ * a[1,1] + ξ₂ * a[1,2] + ξ₃ * a[1,3]
            y[i,j] = ξ₁ * a[2,1] + ξ₂ * a[2,2] + ξ₃ * a[2,3]
        end
        for j = n-i+1:n
            x[i,j] = NaN
            y[i,j] = NaN
        end
    end
    return x, y
end

function bary_func(f, n, a)
    z = OffsetArray{Float64}(undef, 0:n, 0:n)
    pts = OffsetArray(range(0, 1, length=n+1), 0:n)
    ξ = Vector{Float64}(undef, 3)
    for i = 0:n
        ξ[1] = pts[i]
        for j = 0:n-i
            ξ[2] = pts[j]
            ξ[3] = 1 - ξ[1] - ξ[2]
            z[i,j] = f(ξ)
        end
    end
    return z
end

function bary_plot(x, y, z)
    n = size(x, 1) - 1
    xd = [ x[i,n-i] for i = 0:n ]
    yd = [ y[i,n-i] for i = 0:n ]
    zd = [ z[i,n-i] for i = 0:n ]
    idx = [1, 2, 3, 1]
    plot3D(a[1,idx], a[2,idx], zeros(4), "k")
    for i = 0:n
        plot3D(x[i,0:n-i], y[i,0:n-i], z[i,0:n-i], "c")
    end
    for j = 0:n
        plot3D(x[0:n-j,j], y[0:n-j,j], z[0:n-j,j], "c")
    end
    for i = 0:n
        plot3D(xd, yd, zd, "c")
    end
    axis("off")
#    xlabel(L"$x$")
#    ylabel(L"$y$")
end

function midpoints(a)
    prev = [ 3, 1, 2 ]
    succ = [ 2, 3, 1 ]
    m = Matrix{Float64}(undef, 2, 3)
    for i = 1:3
        m[:,i] = ( a[:,prev[i]] + a[:,succ[i]] ) / 2
    end
    return m
end

ψ₁(ξ) = 2 * ξ[1] * ( ξ[1] - 0.5 )
ψ₂(ξ) = 2 * ξ[2] * ( ξ[2] - 0.5 )
ψ₃(ξ) = 2 * ξ[3] * ( ξ[3] - 0.5 )
ψ₄(ξ) = 4 * ξ[2] * ξ[3]
ψ₅(ξ) = 4 * ξ[3] * ξ[1]
ψ₆(ξ) = 4 * ξ[1] * ξ[2]

a = Float64[ 0  1  2
             2  0  3 ]
m = midpoints(a)
n = 20
x, y = bary_grid(n, a)

figure(1)
z = bary_func(ψ₁, n, a)
bary_plot(x, y, z)
plot3D(a[1,1], a[2,1], "ro")
savefig("psi1.pdf")

figure(2)
z = bary_func(ψ₂, n, a)
bary_plot(x, y, z)
plot3D(a[1,2], a[2,2], "ro")
savefig("psi2.pdf")

figure(3)
z = bary_func(ψ₃, n, a)
bary_plot(x, y, z)
plot3D(a[1,3], a[2,3], "ro")
savefig("psi3.pdf")

figure(4)
z = bary_func(ψ₄, n, a)
bary_plot(x, y, z)
plot3D(m[1,1], m[2,1], "ro")
savefig("psi4.pdf")

figure(5)
z = bary_func(ψ₅, n, a)
bary_plot(x, y, z)
plot3D(m[1,2], m[2,2], "ro")
savefig("psi5.pdf")

figure(6)
z = bary_func(ψ₆, n, a)
bary_plot(x, y, z)
plot3D(m[1,3], m[2,3], "ro")
savefig("psi6.pdf")
