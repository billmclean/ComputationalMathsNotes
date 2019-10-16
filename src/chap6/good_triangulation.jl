using PyPlot
using Triangulations

N = Float64[ -1  1  2  1 -1 -2  0
             -1 -1  0  1  1  0  0 ]

T = [ 1  1  2  3  4  5
      7  2  3  4  5  6
      6  7  7  7  7  7 ]

E = [ 1  2  3  4  5  6
      2  3  4  5  6  1 ]

tri = Triangulation(N, T, E)

figure(1)
draw_triangles(tri)
enumerate_triangles(tri)
enumerate_vertices(tri)
enumerate_bdry_edges(tri, 0.1)
local_enumerate_vertices(tri, 0.3)
axis("off")
axis("equal")
