from scipy.sparse.linalg import eigsh, lobpcg
from helper_funcs import *
from const import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
from qutip import Bloch,basis,sigmax,sigmay,sigmaz,sigmam,qeye,mesolve,expect,Qobj


os.environ["MKL_NUM_THREADS"] = "8" 
os.environ["OMP_NUM_THREADS"] = "8"

cells = 1
Lx = WSe2.Lm*cells
Ly = WSe2.Lm*cells
Nx = 60*cells
Ny = 60*cells
dx = Lx/Nx
dy = Ly/Ny
meV_si = (10**18)*hbar**2/(2*WSe2.M*eV)

neigh_x = neigh_y = -meV_si/dx**2
center = -2*(neigh_x + neigh_y)
phase_x = np.exp(1j*0.0)
phase_y = np.exp(1j*0.0)


X,Y,V = triangular_bounds(WSe2, Nx, Ny, dx, dy)
V_vec = V.ravel()
print("bounds done")

kpath, labels = build_kpath(WSe2, Nk=15)
bands = build_bandstructure(V_vec, Nx,Ny,dx,dy,WSe2,kpath)

visualize_bandstructure(bands,labels,Nk=15)
visualize_bilayer_3d(X,Y,V,WSe2, WS2, 2, cells=cells)
'''
H = build_hamiltonian(V_vec, Nx,Ny,neigh_x, neigh_y, center, WSe2, phase_x, phase_y)
print("hamiltonian done")

E,psi = eigsh(H,k=5, which="SM", mode='normal')

#eigX = np.random.rand(H.shape[0], 4)
#E, psi = lobpcg(H, eigX, largest=False, tol=1e-8)


print("eig done")


n = 0  # Index of eigenstate
psi_n = np.real(psi[:, n]).reshape((Nx, Ny))
prob_density = psi_n**2
phase = np.angle(psi[:, n].reshape((Nx, Ny)))


#moire_qubit(E,psi)

visualize_twisted_bilayer(WSe2, WS2, 2, cells=10)

visualize_bilayer_3d(X,Y,V,WSe2, WS2, 2, cells=cells)

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

plt.tight_layout()
plt.show()
'''

'''
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(X, Y, prob_density, facecolors=plt.cm.plasma(prob_density / prob_density.max()),
                       rstride=1, cstride=1, linewidth=0, antialiased=False, alpha=0.9)
ax.plot_surface(X, Y, V, cmap='Greys', alpha=0.3)
ax.set_title(f'3D Wavefunction Density over Moiré Potential\nE = {E[n]:.3f} eV')
ax.set_xlabel('x-grid index')
ax.set_ylabel('y-grid index')
ax.set_zlabel('|ψ|² (Wavefunction Density)')
m = plt.cm.ScalarMappable(cmap='plasma')
m.set_array(prob_density)
fig.colorbar(m, ax=ax, shrink=0.5, aspect=10, label='|ψ|²')

plt.tight_layout()
plt.show()
'''

'''

t_x = neigh_x
t_y = neigh_y

# Brillouin zone path or grid
kx = np.linspace(-np.pi, np.pi, 400)
ky = np.linspace(-np.pi, np.pi, 400)

# Generate energy dispersion
E = np.zeros((len(kx), len(ky)))

for i, kxi in enumerate(kx):
    for j, kyi in enumerate(ky):
        E[i, j] = np.real(build_Hk(kxi, kyi, t_x, t_y, center))  # scalar here

# Plot
plt.figure(figsize=(6,5))
plt.contourf(kx, ky, E.T, levels=100, cmap='viridis')
plt.colorbar(label='Energy')
plt.xlabel('$k_x$')
plt.ylabel('$k_y$')
plt.title('Band structure: $E(k_x, k_y)$')
plt.show()
'''
