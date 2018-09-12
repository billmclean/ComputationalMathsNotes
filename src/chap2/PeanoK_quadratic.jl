using PyPlot

pospart(x) = ( x > 0 ? x : zero(x) )

function K2(x, y, a, b)
    pt = [a, (a+b)/2, b]
    ℓ2 = (x-pt[1])*(x-pt[3])/((pt[2]-pt[1])*(pt[2]-pt[3]))
    ℓ3 = (x-pt[1])*(x-pt[2])/((pt[3]-pt[1])*(pt[3]-pt[2]))
    Qπ = pospart(pt[2]-y)^2 * ℓ2 / 2 + (b-y)^2 * ℓ3 / 2
    return Qπ - pospart(x-y)^2/2
end

a = -1.0
b = 1.0
N = 200
x = LinRange(a, b, N+1)
y = LinRange(a, b, N+1)
z = Float64[ K2(x[k], y[j], a, b) for j=1:N+1, k=1:N+1 ]
fig = figure(1)
surf(x, y, z, cstride=4, rstride=4, cmap="jet", linewidth=0.25)
xlabel(L"x")
ylabel(L"y")
ax = gca()
ax[:view_init](elev=40,azim=-100)
savefig("Quadratic_PeanoK_3d.pdf")

figure(2)
contour(x, y, z, 12)
colorbar()
xlabel(L"x")
ylabel(L"y")
grid(true)
savefig("Quadratic_PeanoK_contour.pdf")
