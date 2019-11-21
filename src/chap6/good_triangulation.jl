using PyPlot
using Triangulations

N = Float64[ 0  -2  2  2  1 -1 -2 
             0  -2 -2  0  2  2  0 ]

T = [ 7  2  3  4  5  6
      2  3  4  5  6  7 
      1  1  1  1  1  1 ]

E = [ 7  2  3  4  5  6 
      2  3  4  5  6  7 ]

tri = Triangulation(N, T, E, 4)

figure(1)
draw_triangles(tri)
enumerate_triangles(tri)
enumerate_vertices(tri)
enumerate_bdry_edges(tri, 0.1)
#local_enumerate_vertices(tri, 0.3)
mark_first_local_vertex(tri, 0.3)
axis("off")
axis("equal")
savefig("good_triangulation.pdf")
