using FEM_P2
using PyPlot

const L = 1.0

mesh = Mesh([0.0, 0.25, 0.6, L], (Dirichlet, Dirichlet))
A = assemble_matrix(mesh, deriv_times_deriv!, x -> 1.0 + x/2)
F = assemble_vector(mesh, x -> 2.0)
ufree = A[:,1:mesh.ndof] \ F
ufix = [0.0, 0.0]
Und = [ ufree; ufix ]
x, u = evaluate_from_ndvals(mesh, Und)
xnd, und = evaluate_from_ndvals(mesh, Und, 2)

figure(1)
plot(x, u, xnd, und, "o")
