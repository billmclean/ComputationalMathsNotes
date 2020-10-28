using PyPlot
using Triangulate

N = [ -1  -1  -1  0   0  1  0
       2   0  -2  0  -2  0  2 ]

T = [ 6  7  1  2  3  4
      7  2  2  3  5  5
      4  4  7  4  4  6 ]

E = [ 7  1  2  3  5  6
      1  2  3  5  6  7 ]

QN = 4

tri = Triangulation(N, T, E, 4)

figure(1)
draw_triangles(tri)
enumerate_triangles(tri)
enumerate_vertices(tri)
enumerate_bdry_edges(tri, 0.1)
mark_first_local_vertex(tri, 0.3)
axis("off")
savefig("example_triangulation.pdf")
