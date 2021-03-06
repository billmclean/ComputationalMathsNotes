using FEM_P2
using PyPlot

bc = (Dirichlet, Neumann)
L = 2.0
γ0 = 1.0
γL = -1.5
a = 1.0
c = 0.0
f(x) = 2.0*a
exact_u(x) = γ0 + (γL+2L)*x - x^2

mesh = Mesh(L .* [0.0, 0.25, 0.6, 1.0], bc)
A = assemble_matrix(mesh, deriv_times_deriv!, a)
F = assemble_vector(mesh, f)
G = Neumann_bc_vector(mesh, [γL])
F .+= G
ufix = [γ0]
ndof = mesh.ndof
F .-= A[:,ndof+1:end] * ufix

ufree = A[:,1:ndof] \ F
Und = [ ufree; ufix ]
x, u = evaluate_from_ndvals(mesh, Und)
xnd, und = evaluate_from_ndvals(mesh, Und, 2)

figure(1)
plot(x, u, "y-" , xnd, und, "o", x, exact_u.(x), "r:")
legend(("FEM", "Nodal values", "Exact solution"), 
       loc="lower center")
grid(true)
xlabel(L"x")
title("Quadratic Finite Elements")
