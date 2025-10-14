import scipy.sparse as sp
import scipy.sparse.linalg as spla
import numpy as np
from const import *

#simulation bounds for base grid
Lx = WSe2.Lm
Ly = WSe2.Lm 
Nx = 80
Ny = 80
dx = Lx/Nx
dy = Ly/Ny
meV_si = (10**21)*hbar**2/(2*WSe2.M*eV)

#implement bloch bounds
"""
def build_bloch(Nx,Ny,dx,dy,kx,ky,meVnm):

    return 0
    """



#builds 2d triangular bounds

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
    return V







