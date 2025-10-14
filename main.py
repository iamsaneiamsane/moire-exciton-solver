from helper_funcs import *
from const import *

Lx = WSe2.Lm
Ly = WSe2.Lm 
Nx = 80
Ny = 80
dx = Lx/Nx
dy = Ly/Ny
meV_si = (10**21)*hbar**2/(2*WSe2.M*eV)

V = triangular_bounds(WSe2, Nx, Ny, dx, dy)

for i in V:
    print(i)