from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt
import numpy as np

x = np.array([0, -2,  2, 2,  1, -1, -2])
y = np.array([0, -2, -2, 0,  2,  2,  0])
triangles = np.array([[0, 6, 1],
                      [0, 1, 2],
                      [0, 2, 3],
                      [0, 3, 4],
                      [0, 4, 5],
                      [0, 5, 6]])
mypoly = Triangulation(x, y, triangles)
n_pts = len(x)
n_tri = triangles.shape[0]

def centroid(j, tri):
    x = tri.x[tri.triangles[j,:]]
    y = tri.y[tri.triangles[j,:]]
    return sum(x)/3, sum(y)/3

plt.triplot(mypoly, color='k')
plt.plot(x, y, 'or')
x_offset = 0.1 * np.array([ 0, -1,  1, 1, 1, -1, -1])
y_offset = 0.1 * np.array([-2, -1, -1, 0, 1,  1,  0])

for j in range(n_pts):
    plt.text(x[j]+x_offset[j], y[j]+y_offset[j], j+1, color='r',
             verticalalignment='center', horizontalalignment='center')

for j in range(n_tri):
    xc, yc = centroid(j, mypoly)
    plt.text(xc, yc, j+1, color='b',
             verticalalignment='center', horizontalalignment='center')
  
edges = np.array([[6,1],
                  [1,2],
                  [2,3],
                  [3,4],
                  [4,5],
                  [5,6]])
n_edg = edges.shape[0]
x_offset = 0.1 * np.array([-1,  0, 1, 1, 0, -1])
y_offset = 0.1 * np.array([ 0, -1, 0, 1, 1,  1])
for j in range(n_edg):
    midx = 0.5 * ( x[edges[j,0]] + x[edges[j,1]] )
    midy = 0.5 * ( y[edges[j,0]] + y[edges[j,1]] )
    plt.text(midx+x_offset[j], midy+y_offset[j], j+1, color='g',
             verticalalignment='center', horizontalalignment='center')

ax = plt.gca()
ax.axis([-2.2, 2.2, -2.2, 2.2])
ax.axis('off')
ax.axis('equal')
plt.show()
plt.savefig('good_triangulation.pdf')
