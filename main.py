from scipy.sparse.linalg import eigsh
from helper_funcs import *
from const import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

cells = 3
Lx = WSe2.Lm*cells
Ly = WSe2.Lm*cells
Nx = 120*cells
Ny = 120*cells
dx = Lx/Nx
dy = Ly/Ny
meV_si = (10**21)*hbar**2/(2*WSe2.M*eV)

neigh_x = neigh_y = -meV_si/dx**2
center = -2*(neigh_x + neigh_y)
phase_x = np.exp(1j*0.0)
phase_y = np.exp(1j*0.0)


X,Y,V = triangular_bounds(WSe2, Nx, Ny, dx, dy)
V_vec = V.ravel()

H = build_hamiltonian(V_vec, Nx,Ny,neigh_x, neigh_y, center, phase_x, phase_y)

E,psi = eigsh(H,k=5, which="SM")



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

plt.figure(figsize=(6,5))
plt.pcolormesh(V, cmap='viridis', shading='auto')
plt.colorbar(label='Potential Energy (eV)')
plt.title('Moiré Potential Landscape')
plt.xlabel('x-grid index')
plt.ylabel('y-grid index')
plt.show()


# Reshape the nth eigenvector
n = 0  # choose eigenstate index
psi_n = np.real(psi[:, n]).reshape((Nx, Ny))

plt.figure(figsize=(6,5))
plt.pcolormesh(psi_n, cmap='RdBu_r', shading='auto')
plt.colorbar(label='Wavefunction amplitude')
plt.title(f'Eigenstate {n}, E = {E[n]:.3f} eV')
plt.xlabel('x-grid index')
plt.ylabel('y-grid index')
plt.show()

plt.figure(figsize=(6,5))
plt.contourf(V, cmap='Greys', alpha=0.6)
plt.contour(psi_n**2, levels=10, cmap='plasma')
plt.colorbar(label='|ψ|² probability density')
plt.title(f'Wavefunction Density over Moiré Potential, E = {E[n]:.3f} eV')
plt.show()

plt.figure(figsize=(6,5))
plt.imshow(np.angle(psi_n.reshape((Nx, Ny))), cmap='twilight')
plt.colorbar(label='Phase (radians)')
plt.title('Phase distribution of eigenstate')
plt.show()
