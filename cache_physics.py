import numpy as np
import pickle
import os
from const import *
from helper_funcs import *
import time

cache_dir = "cache"

theta_start = 0.5
theta_end = 5
step = .1
ppnm = 10
ppnmqub = 4

def main():
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)
    
    print(f"--- Starting Unified Physics Export ---")
    print(f"Visual Res: {ppnm} pts/nm | Qubit Res: {ppnmqub} pts/nm")
    
    num_steps = int(round((theta_end - theta_start) / step))
    
    for i in range(num_steps + 1):
        theta = round(theta_start + i * step, 2)
        filename = f"data_{theta:.2f}.pkl"
        filepath = os.path.join(cache_dir, filename)
        
        print(f"[{i+1}/{num_steps+1}] Processing θ = {theta:.2f}° ... ", end="")
        
        if os.path.exists(filepath):
            print("Skipping (Cached)")
            continue
            
        start_t = time.time()
        
        try:
            bands_data, labels, _ = build_bilayer_bands(theta, WS2, WSe2, ppnm=ppnm)
            

            a0 = WS2.a; a1 = WSe2.a
            delta = 2 * abs(a0 - a1) / (a1 + a0)
            Lm = a0 / np.sqrt(np.radians(theta)**2 + delta**2)
            

            Lx_u = Lm
            Ly_u = np.sqrt(3) * Lm
            Nx_u = int(Lx_u * ppnm)
            Ny_u = int(Ly_u * ppnm)
            dx_u = Lx_u / Nx_u
            dy_u = Ly_u / Ny_u
            
            _, _, V_unit = triangular_bounds(Lm, WSe2, Nx_u, Ny_u, dx_u, dy_u, offset=0, theta_deg=theta, psi_deg=WSe2.psi)


            Nx_q = int(3 * Lm)

            Nx_q = int(max(18, Lm * ppnmqub))
            Ny_q = int(max(18, Lm * ppnmqub))
            dx_q = Lm / Nx_q 

            Nx_q = int(3 * Lm) 
            Ny_q = int(3 * Lm)
            if Nx_q < 5: Nx_q = 5
            if Ny_q < 5: Ny_q = 5
            dx_q = 10.0 / Nx_q
            dy_q = 10.0 / Ny_q
            
            phase = np.exp(1j*0.0)
            qubit_setup = build_h_bilayer(Nx_q, Ny_q, WS2, WSe2, dx_q, dy_q, phase_x=phase, phase_y=phase, theta=theta)
            
            E_vals, _ = exciton_solver(Nx_q, Ny_q, qubit_setup['He'], qubit_setup['Hh'], qubit_setup['Uk'], qubit_setup['Vk'], n_eigs=2)

            unified_data = {
                "theta": theta,
                "Lm": float(Lm),
                "bands": {
                    "data": bands_data,
                    "labels": labels,
                    "meta": {"n_valence": 6, "n_conduction": 2}
                },
                "potential": {
                    "V_unit": V_unit.real.tolist(),
                    "Lx": float(Lx_u),
                    "Ly": float(Ly_u),
                    "min": float(np.min(V_unit.real)),
                    "max": float(np.max(V_unit.real))
                },
                "qubit": {
                    "E0": float(E_vals[0] * 1000),
                    "E1": float(E_vals[1] * 1000),
                    "hamiltonian": qubit_setup 
                }
            }
            
            with open(filepath, "wb") as f:
                pickle.dump(unified_data, f)
                
            print(f"Done ({time.time() - start_t:.2f}s)")
            
        except Exception as e:
            print(f"FAILED: {e}")

if __name__ == "__main__":
    main()
