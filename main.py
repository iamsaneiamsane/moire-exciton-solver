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
Nx = 20*cells
Ny = 20*cells
dx = 10/Nx
dy = 10/Ny
meV_si = (10**18)*hbar**2/(2*WSe2.M*m0*eV)

neigh_x = neigh_y = -meV_si/dx**2
center = -2*(neigh_x + neigh_y)
phase_x = np.exp(1j*0.0)
phase_y = np.exp(1j*0.0)


for t in [0.05,0.1, 1]:
    H = build_h_bilayer(Nx, Ny, WSe2, WS2, t, dx, dy, phase_x, phase_y)
    E, psi = exciton_solver(**H, n_eigs=5)
    Ee = spla.eigsh(H['He'], k=1, which='SA', return_eigenvectors=False)[0]
    Ev = spla.eigsh(H['Hh'], k=1, which='SA', return_eigenvectors=False)[0]
    print("single-particle gap =", Ee - Ev, "eV") 
    
    psi_n = psi[:, 0]
    N_layer = Nx * Ny
    w1 = np.sum(np.abs(psi_n[:N_layer])**2)
    w2 = np.sum(np.abs(psi_n[N_layer:])**2)
    print(f"t = {t:.3f} eV: layer weights: {w1:.3f}, {w2:.3f}")
    we,wh = extract_exciton_layer_weights(psi, Nx,Ny)

    for i, (e,h) in enumerate(zip(we,wh)):
        print(f"State:{i}: E={E[i]:.3f} eV, w_e={e:.3f}, w_h={h:.3f}, sum={e+h:.3f}")

print("eig done")
'''
n = 0  # Index of eigenstate
psi_n = np.real(psi[:, n])
psi1 = psi_n[:Nx*Ny].reshape((Nx**2,Ny**2))
psi2 = psi_n[Nx*Ny:].reshape((Nx**2,Ny**2))
prob_density = psi_n**2
phase1 = np.angle(psi1[:, n])
phase2 = np.angle(psi2[:, n])
prob_layer1 = np.sum(np.abs(psi_n[:Nx*Ny])**2)
prob_layer2 = np.sum(np.abs(psi_n[Nx*Ny:])**2)
print(f"Layer 1 weight: {prob_layer1:.3f}, Layer 2 weight: {prob_layer2:.3f}")

plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plt.imshow(psi1, cmap='RdBu_r')
plt.title('Layer 1 wavefunction')

plt.subplot(1,2,2)
plt.imshow(psi2, cmap='RdBu_r')
plt.title('Layer 2 wavefunction')
plt.show()
'''








'''
X,Y,V = triangular_bounds(WSe2, Nx, Ny, dx, dy)
V_vec = V.ravel()
print("bounds done")
print("neigh_x, neigh_y:", neigh_x, neigh_y)
print("center term:", center)
print("Potential range:", np.min(V_vec), np.max(V_vec))


kpath, labels = build_kpath(WSe2, Nk=20)
bands = build_bandstructure(V_vec, Nx,Ny,dx,dy,WSe2,kpath)
valence_idx = 1   
conduction_idx = 2
val_max = bands[:, valence_idx].max()
cond_min = bands[:, conduction_idx].min()
print("valence max:", val_max, "cond min:", cond_min, "gap (eV):", cond_min - val_max)

# bandwidth of band n
for n in range(bands.shape[1]):
    print(f"band {n} bandwidth (eV):", bands[:,n].max() - bands[:,n].min())



visualize_bandstructure(bands,labels,Nk=20)
visualize_bilayer_3d(WSe2, WS2, 3.46, cells=2)


#H = build_h_bilayer( Nx,Ny,WSe2, WS2, .04, dx, dy, phase_x, phase_y)
print("hamiltonian done")

#X,Y,V = triangular_bounds(WSe2, Nx, Ny, dx, dy)
#E,psi = eigsh(H,k=5, which="SM", mode='normal', tol = 1e-4)

#eigX = np.random.rand(H.shape[0], 4)
#E, psi = lobpcg(H, eigX, largest=False, tol=1e-8)






#moire_qubit(E,psi)

visualize_twisted_bilayer(WSe2, WS2, 2, cells=10)

fig = plt.figure(figsize=(18, 12))
grid_spec = fig.add_gridspec(2, 2)

#2D Potential Landscape
ax2 = fig.add_subplot(grid_spec[0, 0])
im2 = ax2.imshow(phase2, cmap='twilight', origin='lower')
fig.colorbar(im2, ax=ax2, label='Phase (radians)')
ax2.set_title('Phase distribution of eigenstate')
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
im5 = ax5.imshow(phase1, cmap='twilight', origin='lower')
fig.colorbar(im5, ax=ax5, label='Phase (radians)')
ax5.set_title('Phase distribution of eigenstate')
ax5.set_xlabel('x-grid index')
ax5.set_ylabel('y-grid index')

plt.tight_layout()
plt.show()



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
