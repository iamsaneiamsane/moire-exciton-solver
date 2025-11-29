from flask import Flask, request, jsonify
from flask_cors import CORS
import numpy as np
from qutip import Qobj, sigmax,sigmay,mesolve,basis
from qc import *
from helper_funcs import *
import pickle
import os

app = Flask(__name__)
CORS(app)
h = .6582 #mev*ps
current_qubit = None
CACHE_DIR = "band_cache"
QUBIT_CACHE_DIR = "qubit_cache"


def init_qubit_logic(theta, cells=1):
    """
    Core logic to initialize the qubit.
    Tries to load Hamiltonian from cache first, then solves for eigenvalues.
    """
    global current_qubit
    theta = round(theta, 2)
    print(f"Initializing Qubit at theta={theta}°...")

    # --- 1. TRY LOADING CACHED HAMILTONIAN ---
    cache_path = os.path.join(QUBIT_CACHE_DIR, f"qubit_{theta:.2f}.pkl")
    setup = None
    
    if os.path.exists(cache_path):
        with open(cache_path, "rb") as f:
            setup = pickle.load(f)
    print(f"  > Solving Eigenvalues...")
    # Extract grid size from setup dict to ensure consistency
    Nx = setup['Nx']
    Ny = setup['Ny']
    E, psi = exciton_solver(setup['Nx'], setup['Ny'], setup['He'], setup['Hh'], setup['Uk'], setup['Vk'], n_eigs=2)

    E0 = E[0] * 1000.0
    E1 = E[1] * 1000.0
    qubit_wrapper = MoireQubit(E0,E1)
    current_qubit = qubit_wrapper
    print(f"--- Qubit Ready: E0={current_qubit.E0:.2f} meV, E1={current_qubit.E1:.2f} meV ---")
    return current_qubit


def make_serializable(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.float32, np.float64)):
        return float(obj)
    if isinstance(obj, (np.int32, np.int64)):
        return int(obj)
    if isinstance(obj, dict):
        return {k: make_serializable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [make_serializable(i) for i in obj]
    return obj

@app.route('/status', methods=['GET'])
def status():
    return jsonify({'status':'online', 'backend':'Qutip+Flask'})

@app.route('/initialize', methods=['POST'])
def initialize():
    global current_qubit
    data = request.json
    theta = float(data.get('theta', 1.0))
    cells = int(data.get('cells', 1))
    
    print(f"Initializing Qubit at theta={theta}°...")

    current_qubit = init_qubit_logic(theta)
    
    return jsonify({
        "success": True, 
        "E0": current_qubit.E0, 
        "E1": current_qubit.E1, 
        "omega_q": current_qubit.oq
    })

@app.route('/bandstructure', methods=['POST'])
def get_bands():
    data = request.json
    raw_theta = float(data.get('theta', 3.0))
    theta = round(raw_theta * 10) / 10
    
    cache_key = f"bands_{theta:.2f}.pkl"
    cache_path = os.path.join(CACHE_DIR, cache_key)
    
    if os.path.exists(cache_path):
        print(f"Serving cached bands for {theta}°")
        with open(cache_path, 'rb') as f:
            return jsonify(make_serializable(pickle.load(f)))
    
    print(f"Cache miss for {theta}°. Calculating live...")
    try:
        bands_data, labels, potential_data = build_bilayer_bands(theta, WS2, WSe2, ppnm=2.0)
        
        return jsonify(make_serializable({
            "theta": theta,
            "bands": {
                "data": bands_data,
                "labels": labels,
                "meta": {"n_valence": 6, "n_conduction": 2}
            },
            "potential": potential_data
        }))
    except Exception as e:
        print(e)
        return jsonify({"error": str(e)}), 500

@app.route('/evolve', methods=['POST'])
def handle_evolution():
    if not current_qubit: 
        return jsonify({"error": "System initializing..."}), 503
        
    data = request.json
    
    # State
    pd = data['psi']
    c1 = pd['c1']['r'] + 1j*pd['c1']['i']
    c0 = pd['c0']['r'] + 1j*pd['c0']['i']
    psi0 = c1*basis(2,0) + c0*basis(2,1)
    if psi0.norm() > 1e-6: psi0 = psi0.unit()

    # Params
    omega = float(data['omega'])
    delta = float(data['delta'])
    duration = float(data['duration'])
    g = float(data.get('g', 0.0))
    
    t1_val = float(data.get('T1', 0))
    T1 = t1_val if t1_val > 0 else None
    
    pulse_width = data.get('pulse_width', None)
    if pulse_width is not None: pulse_width = float(pulse_width)

    t, res, final = current_qubit.evolve(
        psi0, omega, delta, duration, 
        pulse_width=pulse_width, 
        g=g, T1=T1
    )
    
    traj = []
    for i, time in enumerate(t):
        pt = {
            'time': time,
            'sx': res.expect[0][i], 'sy': res.expect[1][i], 'sz': res.expect[2][i],
            'population': res.expect[3][i]
        }
        if g > 0: pt['photons'] = res.expect[4][i] 
        traj.append(pt)

    if final.type == 'oper': 
        pop1 = final[0,0].real
        pop0 = final[1,1].real
        c1_mag = np.sqrt(pop1)
        c0_mag = np.sqrt(pop0)
        
        coherence = final[0,1]
        phase = np.angle(coherence) if np.abs(coherence) > 1e-9 else 0
        c1_val = c1_mag * np.exp(1j * phase)
        c0_val = c0_mag 
    else: 
        c1_val = final.full()[0][0]
        c0_val = final.full()[1][0]

    return jsonify(make_serializable({
        'trajectory': traj, 
        'final_psi': {'c1': c1_val, 'c0': c0_val}
    }))

@app.route('/emission', methods=['POST'])
def handle_emission():
    if not current_qubit: return jsonify({"error": "Not initialized"}), 400
    data = request.json
    
    g = float(data.get('g', 0.5))
    kappa = float(data.get('kappa', 2.0))
    gamma = float(data.get('gamma', 0.1))
    
    t, atom, photon = current_qubit.simulate_emission(g, kappa, gamma, dur=50)
    
    return jsonify(make_serializable({
        'trajectory': [{'time': t[i], 'atom': atom[i], 'photon': photon[i]} for i in range(len(t))]
    }))

if __name__ == '__main__':
    print("--- Moiré Qubit Server Running on Port 5000 ---")
    app.run(port=5000, debug=True)

