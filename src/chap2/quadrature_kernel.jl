using PyPlot

pi5(x, y) = x > y ? (x-y)^5/120 : zero(x)
w1 = 1/6
w2 = 5/6
a = 1/sqrt(5)

K(y) = (1-y)^6/720 - ( w2*pi5(-a,y) + w2*pi5(a,y) + w1*pi5(1,y) )

figure(1)
y = range(-1, stop=1, length=201)
plot(y, K.(y))
grid(true)
xlabel(L"$y$")
ylabel(L"$K(y)$")
savefig("quadrature_kernel.pdf")

