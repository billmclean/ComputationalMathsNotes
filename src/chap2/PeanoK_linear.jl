using PyPlot

pospart(x) = ( x > 0 ? x : zero(x) )

function K1(x, y, a, b)
    ℓ2 = ( x - a ) / ( b - a )
    Qπ = (b-y) * ℓ2
    return (b-y)*ℓ2 - pospart(x-y)
end

a = 0.0
b = 1.0
N = 200
x = LinRange(a, b, N+1)
y = LinRange(a, b, N+1)
z = Float64[ K1(x[k], y[j], a, b) for j=1:N+1, k=1:N+1 ]
fig = figure(1)
surf(x, y, z, cstride=4, rstride=4, cmap="jet", linewidth=0.5)
xlabel(L"x")
ylabel(L"y")
ax = gca()
ax[:view_init](elev=40,azim=-100)
savefig("Linear_PeanoK_3d.pdf")


figure(2)
contour(x, y, z)
colorbar()
xlabel(L"x")
ylabel(L"y")
grid(true)
savefig("Linear_PeanoK_contour.pdf")
