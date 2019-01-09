using FEM_P2
using PyPlot

L = 2.0
γ0 = 2.0
γL = 1.0
a = 1.0
c = 1.0

mesh = Mesh([0.0, 0.25, 0.6, 1.0], (Dirichlet, Neumann))

A_free, A_fix = assemble_matrix(mesh, deriv_times_deriv!, a)
C_free, C_fix = assemble_matrix(mesh, func_times_func!, c)
F = assemble_vector(mesh, source_times_func!, x -> 1.0)
G = Neumann_bc_vector(mesh, [γL])
A = A_free + C_free
b = ( F + G ) - ( A_fix + C_fix ) * u_fix
u_free = A \ b
u_fix = [ γ0 ]
x, u = evaluate(mesh, 10, u_free, u_fix)

figure(1)
plot(x, u)
