from scipy.sparse.linalg import eigsh, lobpcg
from helper_funcs import *
from const import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from qutip import Bloch,basis,sigmax,sigmay,sigmaz,sigmam,qeye,mesolve,expect,Qobj
from qc import *

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



H = build_h_bilayer(Nx,Ny,WS2,WSe2,dx,dy,phase_x,phase_y,3)
E,psi = exciton_solver(Nx,Ny,H['He'], H['Hh'], H['Uk'], H['Vk'])




H_q, psi_K, psi_Kp, b = moire_qubit(E, psi, Nx, Ny)

print("Valley qubit Hamiltonian (meV):")
print(H_q * 1e3)
overlap = np.vdot(psi_K.full().flatten(), psi_Kp.full().flatten())
print(f"⟨K|K′⟩ = {overlap:.3f}") 
b.show()
plt.show(block=True)

qubit = MoireQubit(H_q[0,0].real, H_q[1,1].real)
t,r = qubit.rabi(2,0,5)
qubit.plot_pop(t,r)
qubit.plot_bloch(r)
