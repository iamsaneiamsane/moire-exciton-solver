from helper_funcs import *
from const import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


Lx = WSe2.Lm*3
Ly = WSe2.Lm*3
Nx = 240
Ny = 240
dx = Lx/Nx
dy = Ly/Ny
meV_si = (10**21)*hbar**2/(2*WSe2.M*eV)

X,Y,V = triangular_bounds(WSe2, Nx, Ny, dx, dy)


fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, V, cmap='viridis', edgecolor='none')
ax.set_title("3D Surface Plot of Triangular Potential")
ax.set_zlim(-100, 100) 
ax.set_xlabel('x(nm)')
ax.set_ylabel('y(nm)')
ax.set_zlabel('mV(x, y)')
plt.tight_layout()
plt.show()
