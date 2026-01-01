from const import *
# from helper_funcs import * # Unused
import numpy as np
# import matplotlib.pyplot as plt # Unused
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
    
    
    def evolve(self, psi0, rabi_amp, detuning, dur, pulse_width = None, steps=100, g=0.0, T1=None, T2=None):
        hbar = .6582 #meV*ps
        cav=5

        H = self.get_hamiltonian(rabi_amp, detuning, g, cav)
        if dur <= 0: dur = 0.1
        if pulse_width is None or pulse_width >=dur:
            pulse_width = dur
        s1 = int(steps*pulse_width/dur)
        if s1 < 2: s1 = 2
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
            H2 = self.get_hamiltonian(0.0, 0.0, g, cav)
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

        fstatef = res.states[-1]
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
        tlist = np.linspace(0,dur,100)

        eops = [tensor(basis(2,0)*basis(2,0).dag(), qeye(cav)), a.dag()*a]
        result = mesolve(H, psi0, tlist, c_ops=cops, e_ops=eops, options=Options(nsteps=5000))
        return tlist, result.expect[0], result.expect[1]

    def get_emission_statistics(self, g, kappa, gamma):
        g_ = g/self.hbar
        kap = kappa/self.hbar
        gam = gamma/self.hbar
        cav = 6 

        sm = tensor(sigmam(), qeye(cav))
        sp = tensor(sigmap(), qeye(cav))
        sz = tensor(sigmaz(), qeye(cav))
        a = tensor(qeye(2), destroy(cav))
        
        H = g_*(a.dag()*sm + a*sm.dag())
        U = g_ * 1.0 
        H_kerr = U * a.dag() * a.dag() * a * a
        H += H_kerr

        Omega = kap * 1.5
        H_drive = Omega * (sm + sp)
        H += H_drive

        cops = [np.sqrt(kap)*a, np.sqrt(gam)*sm]
        
        rho_ss = steadystate(H, cops)
        wlist = np.linspace(-5*g_, 5*g_, 80)

        try:
            spec = spectrum(H, wlist, cops, a.dag(), a) 
        except:
            spec = np.zeros_like(wlist)
        
        # g2(tau)
        tlist_g2 = np.linspace(0, 50, 40)
        g2_tau = coherence_function_g2(H, rho_ss, tlist_g2, cops, a)[0]
        
        # Wigner Function
        rho_cav = rho_ss.ptrace(1)
        xvec = np.linspace(-3,3,30)
        W = wigner(rho_cav, xvec, xvec)
        
        return {
            'spectrum': {'w': wlist, 'S': spec},
            'g2': {'tau': tlist_g2, 'val': np.real(g2_tau)},
            'wigner': {'x': xvec, 'W': W}
        }