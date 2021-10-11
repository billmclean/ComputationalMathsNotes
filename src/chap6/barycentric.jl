using PyPlot
using LinearAlgebra: norm, dot

a = [ 0.5  -1.0   0.25
      1.0   0.5  -1.0 ]
c = (1/3) * sum(a, dims=2)

figure(1)

function level_sets(a, j, ξ)
    prev = [3, 1, 2]
    succ = [2, 3, 1]
    m = length(ξ)
    start = zeros(2, m)
    finish = zeros(2, m)
    for k = 1:m
        η = -0.3
        start[:,k] = ( a[:,prev[j]] + ξ[k] * ( a[:,j] - a[:,prev[j]] ) 
                      + η * ( a[:,succ[j]] - a[:,prev[j]] ) )
        η = 1.1 - ξ[k] 
        finish[:,k] = ( a[:,prev[j]] + ξ[k] * ( a[:,j] - a[:,prev[j]] ) 
                      + η * ( a[:,succ[j]] - a[:,prev[j]] ) )
    end
    return start, finish
end

function draw_level_sets(a, j)
    offset = [ 0.1   0.0  -0.1 
               0.0   0.1  -0.2 ]
    s = latexstring("\\xi_", j, "=")
    vals = ["0", "1/3", "2/3", "1"]
    ha = ["left", "center", "right"]
    va = ["center", "center", "center"]
    ξ = [0, 1/3, 2/3, 1]
    start, finish = level_sets(a, j, ξ)
    for k = 1:4
        plot([ start[1,k], finish[1,k] ], [ start[2,k], finish[2,k] ], "b",
            linewidth=0.5)
        text(start[1,k]+offset[1,j], start[2,k]+offset[2,j],
             latexstring(s, vals[k]), color="b",
             horizontalalignment=ha[j], verticalalignment=va[j])
    end
end

function draw_triangle(a, c, d)
    J = [1, 2, 3, 1]
    plot(a[1,J], a[2,J], "-k")
    for j = 1:3
        s = latexstring("a_", j)
        offset = a[:,j] - c
        offset *= d / norm(offset)
        text(a[1,j]+offset[1], a[2,j]+offset[2], s, color="r",
            horizontalalignment="center", verticalalignment="center")
    end
    plot(a[1,:], a[2,:], "ro", c[1], c[2], "ko")
end

for j = 1:3
   draw_level_sets(a, j)
end
draw_triangle(a, c, 0.2)
text(c[1]+0.06, c[2]-0.02, L"c",
     horizontalalignment="left", verticalalignment="center")
axis("equal")
axis("off")
savefig("barycentric.pdf")

figure(2)

function draw_vectors(a, b, c)
    prev = [3, 1, 2]
    succ = [2, 3, 1]
    quiver([ c[1], c[1], c[1] ], [ c[2], c[2], c[2] ], b[1,:], b[2,:],
           angles="xy", scale_units="xy", scale=1, 
           width=0.004, headwidth=5, color="g")
    offset = [ -0.1   0.1   -0.1
                0.05  0.07   0.0 ]
    for j = 1:3
        text(c[1]+b[1,j]+offset[1,j], c[2]+b[2,j]+offset[2,j],
                 latexstring("b_", j), color="g",
                 horizontalalignment="center", verticalalignment="center")
        bhat = b[:,j] / norm(b[:,j])
        λ = dot(c - a[:,succ[j]], bhat)
        foot = c - λ * bhat
        plot([c[1], foot[1]], [c[2], foot[2]], "--k")
        side = a[:,succ[j]] - a[:,prev[j]]
        v = side / norm(side)
        plot(foot[1] .+ 0.1*[ bhat[1], bhat[1]+v[1], v[1]],
             foot[2] .+ 0.1*[ bhat[2], bhat[2]+v[2], v[2]], linewidth=0.5, "-k")
    end
end

A = [ (a[:,1]-a[:,3]) (a[:,2]-a[:,3]) ]
B = inv(A')
b = [ B[:,1] B[:,2] -sum(B, dims=2) ]

draw_vectors(a, b, c)
draw_triangle(a, c, 0.15)
axis("equal")
axis([-2.4, 1.4, -2.4, 2.4])
axis("off")
savefig("b_vectors.pdf")

figure(3)

function draw_quadratic_triangle(a, c, d)
    J = [1, 2, 3, 1]
    plot(a[1,J], a[2,J], "-k")
    m = zeros(2, 3)
    m[:,1] = ( a[:,3] + a[:,2] ) / 2
    m[:,2] = ( a[:,1] + a[:,3] ) / 2
    m[:,3] = ( a[:,2] + a[:,1] ) / 2
    for j = 1:3
        s = latexstring("a_", j)
        offset = a[:,j] - c
        offset *= d / norm(offset)
        text(a[1,j]+offset[1], a[2,j]+offset[2], s, color="r",
            horizontalalignment="center", verticalalignment="center")
        s = latexstring("m_", j)
        offset = m[:,j] - c
        offset *= d / norm(offset)
        text(m[1,j]+offset[1], m[2,j]+offset[2], s, color="r",
            horizontalalignment="center", verticalalignment="center")
    end
    plot(a[1,:], a[2,:], "ro", m[1,:], m[2,:], "ro")
end

draw_quadratic_triangle(a, c, 0.15)
axis("equal")
axis("off")
savefig("midpoints.pdf")
