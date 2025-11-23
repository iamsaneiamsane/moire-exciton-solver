from const import *
from helper_funcs import *
import numpy as np
import matplotlib.pyplot as plt
from qutip import *

class MoireQubit:
    def __init__(self,E0,E1):
        self.E0 = E0
        self.E1 = E1
        self.oq = E1-E0

        self.H0 = .5*self.oq*sigmaz()
        self.Hx=sigmax()
    
    def rabi(self, ramp, detune, dur, steps=500):
        hbar = .6582 #meV*ps
        H_rot = .5*(detune*sigmaz() + ramp*sigmax())
        tlist = np.linspace(0,dur,steps)

        psi0 = basis(2,1)

        result = mesolve(H_rot, psi0, tlist, [], [sigmaz(),sigmay(),sigmax()])
        return tlist,result
    
    def plot_bloch(self, result):
        b = Bloch()
        b.add_points([result.expect[2], result.expect[1], result.expect[0]], meth="l")
        b.zlabel = [r'$|1\rangle$', r'$|0\rangle$']
        b.show()
        plt.show(block=True)

    def plot_pop(self, tlist, result):
        P_exc = (result.expect[0]+1)/2

        plt.figure(figsize=(8,5))
        plt.plot(tlist,P_exc,'b-',label=r'Excited State Population $\rho_{11}$')
        plt.xlabel('t (ps)')
        plt.ylabel('population')
        plt.title("Rabi Oscillations")
        plt.grid(True, linestyle='--', alpha=.6)
        plt.legend()
        plt.tight_layout()
        plt.show()