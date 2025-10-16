import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse import diags
import numpy as np
from qutip import basis,sigmax,sigmay,sigmaz,sigmam,qeye,mesolve,expect,Qobj
from const import *
from numba import njit, prange

#def bloch_bounds(material, Nx,Ny,dx,dy): #TODO 2d block bounds 





def triangular_bounds(material, Nx,Ny,dx,dy): #builds 2d triangular bounds
    xs = (np.arange(Nx) + .5)*dx
    ys = (np.arange(Ny)+ .5)*dy
    X,Y = np.meshgrid(xs,ys,indexing="ij")
    Gmag = 4*np.pi/(np.sqrt(3)*material.Lm)
    theta = np.array([0,2*np.pi/3, 4*np.pi/3])
    V = np.zeros_like(X)
    for ang in theta:
        Gx = Gmag*np.cos(ang)
        Gy = Gmag*np.sin(ang)
        V+= np.cos((Gx*X + Gy*Y))
    V = material.V0 * V/3
    #V = V.reshape(Nx*Ny)
    return X,Y,V


def hex_lattice(a,cells,z=0):
    ax = np.sqrt(3)*a
    ay = 3*a/2

    xs,ys,zx = [],[],[]
    for i in range(-cells,cells):
        for j in range(-cells,cells):
            x = i*ax + (j%2) * (ax/2)
            y = j*ay
            xs.a

@njit(parallel=True)
def build_hamiltonian_numba(V_vec,Nx,Ny, neigh_x, neigh_y, center, phase_x=1., phase_y=1.):
    N = Nx*Ny
    neighbours = 5
    nbt = neighbours*N
    cols, rows, data = np.zeros(nbt, dtype = np.int32), np.zeros(nbt, dtype = np.int32), np.zeros(nbt, dtype = np.complex128)
    hb2m = (hbar**2 / (2.0 * m0*.6)) * 1e21 / 1.602176634e-19 #KE operator
    


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
        

def build_hamiltonian(V_vec,Nx,Ny, neigh_x, neigh_y, center, phase_x=1., phase_y=1.):
    rows, cols, data = build_hamiltonian_numba(V_vec,Nx,Ny, neigh_x, neigh_y, center, phase_x=1., phase_y=1.)
    H = sp.csr_matrix((data,(rows,cols)),shape=(Nx*Ny,Nx*Ny))
    return H


            







