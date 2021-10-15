using Triangulation_diagrams
using PyPlot

N = [ 0  1  2  4  4  5  3  3  2  1  5  6  6  2  6  4
      2  1  0  0  2  1  1  3  2  3  3  0  4  4  2  4 ]

T = [  1  3  3  5  5
       3  5 12 12 13
      14 14  5 13 14
       9  8  6 15 16
      10  9  7 11  8
       2  7  4  6 11 ]

E = [ 14 10  1  2  3  4 12 15 13 16
      10  1  2  3  4 12 15 13 16 14 ]

QN = 6

tri = Triangulation(N, T, E, QN)

figure(1)
draw_triangles(tri)
enumerate_vertices(tri)
enumerate_triangles(tri)
enumerate_bdry_edges(tri, 0.2)
mark_first_local_vertex(tri, 0.3)
axis("off")
savefig("quad_triangulation.pdf")
