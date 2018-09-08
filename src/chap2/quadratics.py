# Run in isympy

x = symbols('x')

psi = [2*(x-Rational(1,2))*(x-1),
       4*x*(1-x),
       2*x*(x-Rational(1,2))]

dpsi = [ diff(psi[k], x) for k in range(3) ]

def a(p, q):
    return integrate(dpsi[q]*dpsi[p], (x, 0, 1))

def c(p,q):
    return integrate(psi[q]*psi[p], (x, 0, 1))

A = Matrix(3, 3, a)
C = Matrix(3, 3, c)

print('Quadratic element stiffness matrix is 1/(3h) times')
display(3*A)
print('Quadratic element mass matrix is h/30 times')
display(30*C)
