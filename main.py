from scipy.sparse.linalg import eigsh, lobpcg
from helper_funcs import *
from const import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

os.environ["MKL_NUM_THREADS"] = "8" 
os.environ["OMP_NUM_THREADS"] = "8"

cells = 10
Lx = WSe2.Lm*cells
Ly = WSe2.Lm*cells
Nx = 60*cells
Ny = 60*cells
dx = Lx/Nx
dy = Ly/Ny
meV_si = (10**21)*hbar**2/(2*WSe2.M*eV)

neigh_x = neigh_y = -meV_si/dx**2
center = -2*(neigh_x + neigh_y)
phase_x = np.exp(1j*0.0)
phase_y = np.exp(1j*0.0)


X,Y,V = triangular_bounds(WSe2, Nx, Ny, dx, dy)
V_vec = V.ravel()
print("bounds done")

H = build_hamiltonian(V_vec, Nx,Ny,neigh_x, neigh_y, center, phase_x, phase_y)
print("hamiltonian done")

#E,psi = eigsh(H,k=5, which="SM", mode='normal')
k = 5
eigX = np.random.rand(H.shape[0], k)
E, psi = lobpcg(H, eigX, largest=False, tol=1e-8)


print("eig done")

n = 0  # Index of eigenstate
psi_n = np.real(psi[:, n]).reshape((Nx, Ny))
prob_density = psi_n**2
phase = np.angle(psi[:, n].reshape((Nx, Ny)))


#3D Surface Plot
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

fig = plt.figure(figsize=(18, 12))
grid_spec = fig.add_gridspec(2, 2)

#2D Potential Landscape
ax2 = fig.add_subplot(grid_spec[0, 0])
im2 = ax2.pcolormesh(V, cmap='viridis', shading='auto')
fig.colorbar(im2, ax=ax2, label='Potential Energy (eV)')
ax2.set_title('Moiré Potential Landscape')
ax2.set_xlabel('x-grid index')
ax2.set_ylabel('y-grid index')

#Wavefunction Amplitude
ax3 = fig.add_subplot(grid_spec[0, 1])
im3 = ax3.pcolormesh(psi_n, cmap='RdBu_r', shading='auto')
fig.colorbar(im3, ax=ax3, label='Wavefunction amplitude')
ax3.set_title(f'Eigenstate {n}, E = {E[n]:.3f} eV')
ax3.set_xlabel('x-grid index')
ax3.set_ylabel('y-grid index')

#Wavefunction Density over Potential
ax4 = fig.add_subplot(grid_spec[1, 0])
contour_base = ax4.contourf(V, cmap='Greys', alpha=0.6)
contour_overlay = ax4.contour(prob_density, levels=10, cmap='plasma')
fig.colorbar(contour_overlay, ax=ax4, label='|ψ|²')
ax4.set_title(f'Wavefunction Density over Moiré Potential\nE = {E[n]:.3f} eV')
ax4.set_xlabel('x-grid index')
ax4.set_ylabel('y-grid index')

#Phase Distribution
ax5 = fig.add_subplot(grid_spec[1, 1])
im5 = ax5.imshow(phase, cmap='twilight', origin='lower')
fig.colorbar(im5, ax=ax5, label='Phase (radians)')
ax5.set_title('Phase distribution of eigenstate')
ax5.set_xlabel('x-grid index')
ax5.set_ylabel('y-grid index')

#plt.tight_layout()
plt.show()
