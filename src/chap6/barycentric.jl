using PyPlot
using LinearAlgebra: norm

a = [ 1.0  -2.0   0.5
      2.0   1.0  -2.0 ]

figure(1)
J = [1, 2, 3, 1]
plot(a[1,J], a[2,J], "-k")
#hold(true)
function level_sets(a, j)
    prev = [3, 1, 2]
    succ = [2, 3, 1]
    start = zeros(2, 4)
    finish = zeros(2, 4)
    for k = 1:4
        ξ = (k-1)/3
        η = -0.3
        start[:,k] = ( a[:,prev[j]] + ξ * ( a[:,j] - a[:,prev[j]] ) 
                      + η * ( a[:,succ[j]] - a[:,prev[j]] ) )
        η = 1.1 - ξ 
        finish[:,k] = ( a[:,prev[j]] + ξ * ( a[:,j] - a[:,prev[j]] ) 
                      + η * ( a[:,succ[j]] - a[:,prev[j]] ) )
    end
    return start, finish
end

function draw_level_sets(a, j)
    s = latexstring("a_", j)
    c = sum(a, dims=2) / 3 # centroid
    offset = a[:,j] - c
    offset *= 0.4 / norm(offset)
    text(a[1,j]+offset[1], a[2,j]+offset[2], s, color="r",
        horizontalalignment="center", verticalalignment="center")
    offset = [ 0.1   0.0  -0.1 
               0.0   0.1  -0.2 ]
    s = latexstring("\\xi_", j, "=")
    vals = ["0", "1/3", "2/3", "1"]
    ha = ["left", "center", "right"]
    va = ["center", "center", "center"]
    start, finish = level_sets(a, j)
    for k = 1:4
        plot([ start[1,k], finish[1,k] ], [ start[2,k], finish[2,k] ], "b",
            linewidth=0.5)
        text(start[1,k]+offset[1,j], start[2,k]+offset[2,j],
             latexstring(s, vals[k]), color="b",
             horizontalalignment=ha[j], verticalalignment=va[j])
    end
end
for j = 1:3
   draw_level_sets(a, j)
end
c = (1/3) * sum(a, dims=2)
text(c[1]+0.2, c[2]-0.05, L"c",
    horizontalalignment="left", verticalalignment="center")
plot(a[1,:], a[2,:], "ro", c[1], c[2], "ko")
axis("equal")
axis("off")
savefig("barycentric.pdf")
