using FEM_P2
using PyPlot

bc = (Dirichlet, Neumann)
L = 3.0
a(x) = exp(x)
c(x) = exp(-x)
f(x) = c(x) * sin(x) - exp(x) * ( cos(x) - sin(x) ) 
exact_u(x) = sin(x)
γ0 = 0.0
γL = a(L) * cos(L)

#mesh = Mesh(L .* [0.0, 0.25, 0.6, 1.0], bc)
mesh = Mesh(L .* collect(0.0:0.25:1.0), bc)
A = assemble_matrix(mesh, deriv_times_deriv!, a)
C = assemble_matrix(mesh, func_times_func!, c)
A = A + C
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
