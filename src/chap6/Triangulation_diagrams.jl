module Triangulation_diagrams

using PyPlot
using ArgCheck

export Triangulation, draw_triangles, enumerate_triangles, enumerate_vertices,
       enumerate_bdry_edges, local_enumerate_vertices, mark_first_local_vertex

function vertices!(xp, yp, N, T, p)
    for j = 1:3
        xp[j] = N[1, T[j,p]]
        yp[j] = N[2, T[j,p]]
    end
end

function endpoints!(xq, yq, N, E, q)
    for j = 1:2
        xq[j] = N[1, E[j,q]]
        yq[j] = N[2, E[j,q]]
    end
end

function outer_normal!(vecn, N, E, q)
    ax = N[1, E[1,q]]
    ay = N[2, E[1,q]]
    bx = N[1, E[2,q]]
    by = N[2, E[2,q]]
    Δx = bx - ax
    Δy = by - ay
    θ = atan(Δy, Δx)
    vecn[1] = sin(θ)
    vecn[2] = -cos(θ)
end

struct Triangulation
    N :: Matrix{Float64}
    T :: Matrix{Int64}
    E :: Matrix{Int64}
    QN :: Int64
    Centroid :: Matrix{Float64}
    Midpoint :: Matrix{Float64}
    function Triangulation(N, T, E, QN)
        M = size(N, 2)
        P = size(T, 2)
        Q = size(E, 2)
        @argcheck QN ≤ Q
        xp = zeros(3)
        yp = zeros(3)
        Centroid = zeros(2, P)
        for p = 1:P
            vertices!(xp, yp, N, T, p)
            Centroid[1,p] = sum(xp) / 3
            Centroid[2,p] = sum(yp) / 3
        end
        Midpoint = zeros(2, Q)
        xq = zeros(2)
        yq = zeros(2)
        for q = 1:Q
            endpoints!(xq, yq, N, E, q)
            Midpoint[1,q] = sum(xq) / 2
            Midpoint[2,q] = sum(yq) / 2
        end
        return new(N, T, E, QN, Centroid, Midpoint)
    end
end

function draw_triangles(tri::Triangulation, fmt="k")
    N, T, E, QN = tri.N, tri.T, tri.E, tri.QN
    P = size(T, 2)
    xp = zeros(4)
    yp = zeros(4)
    for p = 1:P
        np = T[:,p]
        for j = 1:3
            xp[j] = N[1, np[j]]
            yp[j] = N[2, np[j]]
        end
        xp[4] = xp[1]
        yp[4] = yp[1]
        plot(xp, yp, fmt)
    end
    xq = zeros(2)
    yq = zeros(2)
    for q = QN+1:size(E,2)
       for j = 1:2
          xq[j] = N[1, E[j,q]]
          yq[j] = N[2, E[j,q]]
       end 
       plot(xq, yq, fmt, linewidth=3)
    end
end

function enumerate_triangles(tri::Triangulation)
    Centroid = tri.Centroid
    P = size(Centroid, 2)
    for p = 1:P
        text(Centroid[1,p], Centroid[2,p], string(p), color="b", weight="bold",
            horizontalalignment="center", verticalalignment="center")
    end
end

function enumerate_vertices(tri::Triangulation)
    N = tri.N
    R = size(N, 2)
    for r = 1:R
        x = N[1,r]
        y = N[2,r]
        plot(x, y, "o", markersize=14, markerfacecolor="w",
             markeredgecolor="k")
        text(x, y, string(r), color="r",
            horizontalalignment="center",
            verticalalignment="center")
    end
end

function local_enumerate_vertices(tri::Triangulation, offset)
    N, T, Centroid = tri.N, tri.T, tri.Centroid
    P = size(T, 2)
    xp = zeros(3)
    yp = zeros(3)
    vec = zeros(2)
    for p = 1:P
        vertices!(xp, yp, N, T, p)
        for j = 1:3
            vec[1] = Centroid[1,p] - xp[j]
            vec[2] = Centroid[2,p] - yp[j]
            len = hypot(vec[1], vec[2])
            vec[1] *= offset/len
            vec[2] *= offset/len
            x = xp[j] + vec[1]
            y = yp[j] + vec[2]
            text(x, y, string(j), style="oblique",
#                 bbox=Dict("facecolor"=>"white", "edgecolor"=>"none"),
                 horizontalalignment="center",
                 verticalalignment="center")
        end
    end
end

function mark_first_local_vertex(tri::Triangulation, offset)
    N, T, Centroid = tri.N, tri.T, tri.Centroid
    P = size(T, 2)
    xp = zeros(3)
    yp = zeros(3)
    vec = zeros(2)
    for p = 1:P
        vertices!(xp, yp, N, T, p)
        vec[1] = Centroid[1,p] - xp[1]
        vec[2] = Centroid[2,p] - yp[1]
        len = hypot(vec[1], vec[2])
        vec[1] *= offset/len
        vec[2] *= offset/len
        x = xp[1] + vec[1]
        y = yp[1] + vec[2]
        plot([x], [y], "k*")
    end
end

function enumerate_bdry_edges(tri::Triangulation, offset)
    N, E, Midpoint = tri.N, tri.E, tri.Midpoint
    Q = size(E, 2)
    vecn = zeros(2)
    for q = 1:Q
        outer_normal!(vecn, N, E, q)
        x = Midpoint[1,q] + offset * vecn[1]
        y = Midpoint[2,q] + offset * vecn[2]
        text(x, y, string(q), color="g",
            horizontalalignment="center",
            verticalalignment="center")
    end
end

end # module
