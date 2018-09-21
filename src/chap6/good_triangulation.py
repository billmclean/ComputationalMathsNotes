from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt
import numpy as np

x = np.array([-2,  2, 2, 1, -1, -2, 0])
y = np.array([-2, -2, 0, 2,  2,  0, 0])
triangles = np.array([[0, 6, 5],
                      [0, 1, 6],
                      [1, 2, 6],
                      [2, 3, 6],
                      [3, 4, 6],
                      [4, 5, 6]])
mypoly = Triangulation(x, y, triangles)
n_pts = len(x)
n_tri = triangles.shape[0]

def centroid(j, tri):
    x = tri.x[tri.triangles[j,:]]
    y = tri.y[tri.triangles[j,:]]
    return sum(x)/3, sum(y)/3

plt.triplot(mypoly, color='k')
plt.plot(x, y, 'or')
x_offset = 0.1 * np.array([-1,  1, 1, 1, -1, -1,  0])
y_offset = 0.1 * np.array([-1, -1, 0, 1,  1,  0, -2])

for j in range(n_pts):
    plt.text(x[j]+x_offset[j], y[j]+y_offset[j], j+1, color='r',
             verticalalignment='center', horizontalalignment='center')

for j in range(n_tri):
    xc, yc = centroid(j, mypoly)
    plt.text(xc, yc, j+1, color='b',
             verticalalignment='center', horizontalalignment='center')
  

ax = plt.gca()
ax.axis([-2.2, 2.2, -2.2, 2.2])
ax.axis('off')
plt.show()
plt.savefig('good_triangulation.pdf')