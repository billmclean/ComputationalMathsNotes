using PyPlot
using Triangulate

N = Float64[ 0  -2   0  2  2  1 -1 -2 
             0  -2  -2 -2  0  2  2  0 ]

T = [ 2  3  3  1  1  1
      3  5  4  5  6  7
      8  8  5  6  7  8 ]

E = [ 2 3 4 5 6 7 8
      3 4 5 6 7 8 2 ]

tri = Triangulation(N, T, E, 7)

figure(1)
draw_triangles(tri)
enumerate_triangles(tri)
enumerate_vertices(tri)
#mark_first_local_vertex(tri, 0.3)
axis("off")
axis("equal")
savefig("bad_triangulation.pdf")
