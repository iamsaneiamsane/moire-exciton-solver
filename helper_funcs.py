import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse import diags
import numpy as np
from qutip import Bloch, basis,sigmax,sigmay,sigmaz,sigmam,qeye,mesolve,expect,Qobj
from const import *
from numba import njit, prange
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import RegularPolygon


hb2mb = lambda m_eff: (hbar**2 / (2 * m_eff)) / (eV * 10e-18)


def triangular_bounds(material, Nx,Ny,dx,dy): #builds 2d triangular bounds
    xs = (np.arange(Nx) + .5)*dx
    ys = (np.arange(Ny)+ .5)*dy

    X,Y = np.meshgrid(xs,ys,indexing="ij")

    Gmag = 4*np.pi/(np.sqrt(3)*material.a)
    theta = np.array([0,2*np.pi/3, 4*np.pi/3])
    V = np.zeros_like(X)
    G_vec = np.array([[Gmag, 0.0], [Gmag*np.cos(2*np.pi/3), Gmag*np.sin(2*np.pi/3)], [Gmag*np.cos(4*np.pi/3), Gmag*np.sin(4*np.pi/3)]])
    
    for i,j in G_vec:
        V+= np.cos(i*X + j*Y)

    V = material.V0 * V/3
    #V = V.reshape(Nx*Ny)
    return X,Y,V

def hex_lattice(a, cells, angle=0, offset=(0, 0)):
    ax = np.sqrt(3) * a
    ay = 3 * a / 2

    xs, ys = [], []
    for i in range(-cells, cells):
        for j in range(-cells, cells):
            x = i * ax + (j % 2) * (ax / 2)
            y = j * ay
            x_rot = x * np.cos(np.radians(angle)) - y * np.sin(np.radians(angle))
            y_rot = x * np.sin(np.radians(angle)) + y * np.cos(np.radians(angle))
            xs.append(x_rot + offset[0])
            ys.append(y_rot + offset[1])
    return np.array(xs), np.array(ys)




@njit(parallel=True) 
def build_hamiltonian_numba(V_vec,Nx,Ny, neigh_x, neigh_y, center, hb2M, phase_x=1., phase_y=1.):
    N = Nx*Ny
    neighbours = 5
    nbt = neighbours*N
    cols, rows, data = np.zeros(nbt, dtype = np.int32), np.zeros(nbt, dtype = np.int32), np.zeros(nbt, dtype = np.complex128)
     #KE operator
    
    for x in prange(Nx):
        for y in range(Ny):
            idx = x*Ny + y
            p = idx*neighbours

            #diag
            rows[p] = idx
            cols[p] = idx
            data[p] = center + V_vec[idx]

            #-x neighbour
            xp = (x+1)%Nx
            dxp = xp*Ny + y
            phase = phase_x if x == 0 else 1.0
            rows[p+1] = idx
            cols[p+1] = dxp
            data[p+1] = neigh_x*phase
    
            #+x neighbour
            xm = (x-1)%Nx
            dxm = xm*Nx+y
            phase = np.conj(phase_x) if x == 0 else 1.0
            rows[p+2] = idx
            cols[p+2] = dxm
            data[p+2] = neigh_x*phase
        
            #+y
            yp = (y+1)%Nx
            dyp = x*Ny+yp
            phase = phase_y if y == Ny-1 else 1.0
            rows[p+3] = idx
            cols[p+3] = dyp
            data[p+3] = neigh_y*phase
        
            #-y
            ym = (y-1)%Ny
            dym = x*Ny + ym
            phase = np.conj(phase_y) if y == 0 else 1.0
            rows[p+4] = idx
            cols[p+4] = dym
            data[p+4] = neigh_y*phase
    return rows, cols, data
        



def build_hamiltonian(V_vec,Nx,Ny, neigh_x, neigh_y, center, material = WSe2, phase_x=1., phase_y=1.): #real space hamiltonian
    hb2m = hb2mb(material.M)
    rows, cols, data = build_hamiltonian_numba(V_vec,Nx,Ny, neigh_x, neigh_y, center, hb2m, phase_x=1., phase_y=1.)
    H = sp.csr_matrix((data,(rows,cols)),shape=(Nx*Ny,Nx*Ny))
    return H

def build_Hk_bloch(V_vec, Nx,Ny,dx,dy,material,kx=0.0,ky=0.0):
    hb2m = hb2mb(material.M)

    neigh_x = -hb2m/dx**2
    neigh_y = -hb2m/dy**2
    center = -2*(neigh_x + neigh_y)

    phase_x = np.exp(1j*kx*Nx*dx)
    phase_y = np.exp(1j*ky*Ny*dy)

    rows,cols,data = build_hamiltonian_numba(V_vec, Nx,Ny,neigh_x,neigh_y,center,hb2m,phase_x,phase_y)

    H = sp.csr_matrix((data,(rows,cols)),shape = (Nx*Ny, Nx*Ny))
    return H

def build_bandstructure(V_vec, Nx,Ny,dx,dy,material,k_path,n_eigs = 5):
    bands = np.zeros((len(k_path), n_eigs))

    for i,(x,y) in enumerate(k_path):
        Hk = build_Hk_bloch(V_vec,Nx,Ny,dx,dy,material,x,y)
        E, psi = spla.eigsh(Hk,k=n_eigs, which='SM', maxiter=10000, tol = 1e-9)
        bands[i,:] = np.sort(E.real)
    
    return bands

def segment(p1,p2,N): #helper for kpath
    return [tuple(p1+(p2-p1)*t) for t in np.linspace(0,1,N,endpoint=False)]

def build_kpath(material, Nk=40):
    a = material.a
    G = np.array([0,0]) #gamma point
    K = np.array([4*np.pi/(3*a),0])
    M = np.array([np.pi/a, np.pi/(a*np.sqrt(3))])
    path = segment(G,K,Nk) + segment(K,M,Nk) + segment(M,G,Nk) + [tuple(G)]
    labels = ['Γ', 'K', 'M', 'Γ']
    return path,labels


'''
def miniband(kx_list, ky=0, eig_count = 5, H_kwargs={}):
    bands = np.zeros((len(kx_list),eig_count))

    for i,kx in enumerate(kx_list):
        Hk = build_hamiltonian(Hk, kx=kx, ky=ky, **H_kwargs)
'''
#k space
'''
def build_Hk(kx,ky, neigh_x, neigh_y, center):
    Hk = center + neigh_x*np.exp(-1j*kx) + np.conj(neigh_x)*np.exp(1j*kx) + neigh_y*np.exp(-1j*ky) + np.conj(neigh_y)*np.exp(1j*ky)
    return Hk
'''





#visualize



def visualize_twisted_bilayer(material1, material2, angle, cells=100):
    xs1, ys1 =hex_lattice(material1.a, cells, angle=0)
    xs2, ys2 = hex_lattice(material2.a, cells, angle=angle)

    fig, ax = plt.subplots(figsize=(8, 8))
    
    ax.scatter(xs1, ys1, s=20, c='blue', alpha=0.6, label=f'material1 (bottom)')
    ax.scatter(xs2, ys2, s=20, c='red', alpha=0.6, label=f'material2 (top, θ={angle}°)')

    ax.set_aspect('equal')
    ax.set_title(f"Twisted Bilayer: material2 / material2 @ {angle}°")
    ax.legend()
    ax.axis('off')
    plt.show()

def visualize_bilayer_3d(X,Y,V,material1,material2,angle_deg,cells=4):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, V, cmap='viridis', edgecolor='none', alpha=0.85)
    a1 = material1.a*.1
    a2 = material2.a*.1  

    x_center = (X.max() + X.min()) / 2
    y_center = (Y.max() + Y.min()) / 2

    
    xs1, ys1 = hex_lattice(a1, cells*10, angle=0, offset=(x_center, y_center))
    xs2, ys2 = hex_lattice(a2, cells*10, angle=angle_deg, offset=(x_center, y_center))


    zs1 = np.zeros_like(xs1)
    zs2 = np.full_like(xs2, 10.0)

    ax.scatter(xs1, ys1, zs1, c='blue', s=10, label=material1.__name__)
    ax.scatter(xs2, ys2, zs2, c='red', s=10, label=f"{material2.__name__} ({angle_deg}°)")
    ax.set_title("3D Surface Plot of Triangular Potential")
    ax.set_zlim(-100, 100) 
    ax.set_xlabel('x(nm)')
    ax.set_ylabel('y(nm)') 
    ax.set_zlabel('mV(x, y)')
    plt.tight_layout()
    plt.show()

def visualize_bandstructure(bands,labels,Nk):
    fig,ax = plt.subplots(figsize = (7,5))
    npts = bands.shape[0]
    x = np.arange(npts)
    for n in range(bands.shape[1]):
        ax.plot(x,bands[:,n],color='blue')
    ax.set_xticks([0,Nk,2*Nk,3*Nk])
    ax.set_xticklabels(labels)
    ax.set_ylabel("Energy (eV)")
    ax.grid(True,ls=":")
    plt.show()


def qutip_state(psi): #0ket = layer1, 1let = layer2, or K,K'
    psi/=np.linalg.norm(psi[0:2])
    return Qobj(psi[0:2], dims=[[2],[1]])

def moire_qubit(E,psi):
    E0,E1 = E[:2]
    psi0,psi1 = qutip_state(psi[:,0]), qutip_state(psi[:,1])

    Hq = Qobj(np.diag([E0,E1]))
    b = Bloch()
    b.add_states([psi0, psi1])
    b.show()
    return(Hq)