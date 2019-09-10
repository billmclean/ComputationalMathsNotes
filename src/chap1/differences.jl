using Printf

f(x)  = x^2 + sin(x)
df(x) = 2x + cos(x)
a = 1.0

for k = 1:10
    h = 10.0^(-k)
    central_diff = ( f(a+h) - f(a-h) ) / (2h)
    complex_diff = imag( f(complex(a, h)) ) / h
    deriv = df(a)
    @printf("\$10^{-%d}\$& %10.2e& %10.2e\\\\\n", 
            k, central_diff-deriv, complex_diff-deriv)
end
