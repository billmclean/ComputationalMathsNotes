using PyPlot
using LinearAlgebra

a = [ 0.5  -1.0   0.25
      1.0   0.5  -1.0 ]

function points(λ, offset)
    x = zeros(2,3)
    x[:,1] .= (1-2λ) * a[:,1] +      λ * a[:,2] +      λ * a[:,3]
    x[:,2] .=      λ * a[:,1] + (1-2λ) * a[:,2] +      λ * a[:,3]
    x[:,3] .=      λ * a[:,1] +      λ * a[:,2] + (1-2λ) * a[:,3]
    m = zeros(2,3)
    m[:,1] = 0.5 * a[:,3] + 0.5 * a[:,2]
    m[:,2] = 0.5 * a[:,1] + 0.5 * a[:,3]
    m[:,3] = 0.5 * a[:,2] + 0.5 * a[:,1]
    J = [1, 2, 3, 1]
    plot(a[1,J], a[2,J], "-k")
    c = sum(a, dims=2) / 3
    for j = 1:3
        plot([a[1,j], m[1,j]], [a[2,j], m[2,j]], "--k")
        text(x[1,j] + offset[1,j], x[2,j] + offset[2,j], 
             latexstring("x^{\\langle K\\rangle}_", j), color="r",
             verticalalignment="center", horizontalalignment="center")
        v = a[:,j] - c
        v = 0.1 * v / LinearAlgebra.norm(v)
        text(a[1,j] + v[1] , a[2,j] + v[2], latexstring("a_", j), 
             verticalalignment="center", horizontalalignment="center")
    end
    plot(x[1,:], x[2,:], "or", markeredgecolor="w")
    plot(a[1,:], a[2,:], "ok", markersize=5.0)
end

figure(1)
subplot(1, 2, 1)
offset = [ -0.10  0.10  0.15
            0.05  0.10  0.05  ]
points(1/6, offset)
axis("equal")
axis((-1.2, 0.7, -1.2, 1.2))
axis("off")

subplot(1, 2, 2)
offset = [ -0.12  0.15  0.0
           -0.10  0.0   0.15  ]
points(1/2, offset)
axis("equal")
axis((-1.2, 0.7, -1.2, 1.2))
axis("off")

savefig("quadrature_points.pdf")
