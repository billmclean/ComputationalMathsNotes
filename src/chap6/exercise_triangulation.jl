using PyPlot
using Triangulation_diagrams

N = Float64[ 1  0  1  0  -1   0  -1
             1  1  0  0   0  -1  -1 ]

T = [ 1 4 4 3 4 5
      2 3 2 4 5 7
      3 2 5 6 6 6 ]

E = [ 6 3 1 2 5 7
      3 1 2 5 7 6 ]

tri = Triangulation(N, T, E, 4)

figure(1)
draw_triangles(tri)
enumerate_triangles(tri)
enumerate_vertices(tri)
mark_first_local_vertex(tri, 0.15)
enumerate_bdry_edges(tri, 0.1)
axis("off")
axis("equal")
savefig("ex1_triangulation.pdf")

N = Float64[ 4  2  0  1  3  5  2  4  6  6
             2  2  2  1  1  1  0  0  0  2 ]

T = [ 3  5  4  7  6  1  5   6  10  8
      4  2  7  8  5  2  6   9   1  9
      2  4  5  5  8  5  1  10   6  6 ]

E = [ 10  1  2  3  4  7  8   9
       1  2  3  4  7  8  9  10 ]

tri = Triangulation(N, T, E, 5)

figure(2)
draw_triangles(tri)
enumerate_triangles(tri)
enumerate_vertices(tri)
mark_first_local_vertex(tri, 0.3)
enumerate_bdry_edges(tri, 0.2)
axis("off")
axis("equal")
savefig("ex2_triangulation.pdf")
