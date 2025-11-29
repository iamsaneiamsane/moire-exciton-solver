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
        self.hbar = 0.6582119569 
        self.H0 = 0.5*self.oq*sigmaz()
        self.Hx=sigmax()

    def get_hamiltonian(self, rabi_amp, detuning, g = 0.0, cav = 5):
        rabi_ = rabi_amp/self.hbar
        det_ = detuning/self.hbar
        g_ = g/self.hbar
        if g > 0:
            
            sz = tensor(sigmaz(), qeye(cav))
            sx = tensor(sigmax(), qeye(cav))
            sm = tensor(sigmam(), qeye(cav))
            a = tensor(qeye(2), destroy(cav))

            H0 = det_*(0.5*sz + a.dag()*a)
            Hi = g_*(a.dag()*sm + a*sm.dag())
            Hd = 0.5*rabi_*sx
            return H0+Hi+Hd
        else:
            return .5*(det_*sigmaz() + rabi_*sigmax())
    
    
    def evolve(self, psi0, rabi_amp, detuning, dur, pulse_width = None, steps=1000, g=0.0, T1=None, T2=None):
        hbar = .6582 #meV*ps
        cav=5

        H = self.get_hamiltonian(rabi_amp, detuning, g, cav)
        if dur <= 0: dur = 0.1
        if pulse_width is None or pulse_width >=dur:
            pulse_width = dur
        s1 = int(steps*pulse_width/dur)
        tlist = np.linspace(0, pulse_width, s1)

        if g > 0:
            if psi0.dims != [[2], [1]]: 
                psi_in = tensor(psi0, basis(cav, 0))
            else:
                psi_in = psi0
            eops = [tensor(sigmax(), qeye(cav)),tensor(sigmay(), qeye(cav)),
                    tensor(sigmaz(), qeye(cav)), tensor(basis(2.0)*basis(2.0).dag(), qeye(cav)),
                    tensor(qeye(2), destroy(cav).dag()*destroy(cav))]
            cops = []
        else:
            if psi0.dims != [[2], [1]]:
                psi_in = basis(2, 1)
            else:
                psi_in = psi0
            eops = [sigmax(), sigmay(), sigmaz(), basis(2,0)*basis(2,0).dag()]
            cops = []

        if T1:
            gamma = 1.0/T1
            op = tensor(sigmam(), qeye(cav)) if g > 0 else sigmam()
            cops.append(np.sqrt(gamma)*op) 
        
        if T2:
            gamma_phi = 1.0/T2
            op = tensor(sigmaz(), qeye(cav)) if g > 0 else sigmaz()
            cops.append(np.sqrt(gamma_phi)*op)

        result = mesolve(H, psi_in, tlist, cops, eops, options=Options(store_states=True))

        if pulse_width < dur:
            s2 = steps-s1
            if s2 < 2: s2 = 2
            tlist2 = np.linspace(pulse_width, dur, s2)
            H2 = self.get_hamiltonian(0.0, detuning, g, cav)
            psi_ = result.states[-1]
            result2 = mesolve(H2, psi_, tlist2, cops, eops, options=Options(store_states=True))

            ft = np.concatenate((result.times, result2.times[1:]))

            fe = []
            fs = result.states + result2.states[1:]
            for i in range(len(eops)):
                fe.append(np.concatenate((result.expect[i], result2.expect[i][1:])))
            class unifiedres: pass
            res = unifiedres()
            res.times=ft
            res.expect=fe
            res.states=fs
        else:
            res = result
            ft = tlist

        fstatef = result.states[-1]
        fstatea = fstatef.ptrace(0) if g > 0 else fstatef
        return ft, res, fstatea
    
    def simulate_emission(self, g, kappa, gamma, dur=50):
        g_ = g/self.hbar
        kap = kappa/self.hbar
        gam = gamma/self.hbar
        cav = 5
        sm = tensor(sigmam(), qeye(cav))
        a = tensor(qeye(2), destroy(cav))
        H = g_*(a.dag()*sm + a*sm.dag())

        cops = [np.sqrt(kap)*a, np.sqrt(gam)*sm]

        psi0 = tensor(basis(2,0), basis(cav,0))
        tlist = np.linspace(0,dur,200)

        eops = [tensor(basis(2,0)*basis(2,0).dag(), qeye(cav)), a.dag()*a]
        result = mesolve(H,psi0,tlist,cops,eops)
        return tlist, result.expect[0], result.expect[1]

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