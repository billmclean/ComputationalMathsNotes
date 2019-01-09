using PyPlot

"""
    y = ℓ(c, j, x)

Computes `j`th Lagrange polynomial with respect to the centres `c[k]`.
Thus, if `1 ≤ k ≤ n` then `ℓ` has degree `n-1` and `1 ≤ j ≤ n`.

"""
function ℓ(c, j, x)
    n = length(c)
    @assert 1 ≤ j ≤ n
    y = 1.0
    for k = 1:j-1
        y *= ( x - c[k] ) / ( c[j] - c[k] )
    end
    for k = j+1:n
        y *= ( x - c[k] ) / ( c[j] - c[k] )
    end
    return y
end

function ℓ(c, j, x::AbstractArray)
    n = length(x)
    y = fill(1.0, n)
    for k = 1:n
        y[k] = ℓ(c, j, x[k])
    end
    return y
end

n = 6
x = range(-1; stop=1, length=201)
c = range(-1; stop=1, length=n)

figure(1)
for j = 1:n
    plot(x, ℓ(c, j, x))
end
for j = 1:n
    plot(c[j:j], [1.0], "ok", c[j:j], [0.0], "ok")
end
grid(true)
title("Lagrange polynomials of degree $(n-1)")
