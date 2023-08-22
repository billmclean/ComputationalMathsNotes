from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

x = np.array([-1, 0, 1])
y = np.array([-1, 0, 1])
X, Y = np.meshgrid(x, y)
x, y = X.flatten(), Y.flatten()
mygrid = Triangulation(x, y)

Z = np.zeros_like(X)
Z[1,1] = 1.0
z = Z.flatten()

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.triplot(mygrid, color='k')
ax.plot_trisurf(mygrid, z, color='c', alpha=0.5)
#plt.axis(azimuth=-70, elevation=30)
ax.view_init(azim=-70, elev=30)
plt.xlabel('x')
plt.ylabel('y')

plt.show()
plt.savefig('tent_func.pdf')
