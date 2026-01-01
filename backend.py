from flask import Flask, request, jsonify
from flask_cors import CORS
import numpy as np
from qutip import basis
from qc import *
from helper_funcs import *
import pickle
import os
import json

import sys
import os

app = Flask(__name__)
CORS(app)
h = .6582 #mev*ps
current_qubit = None

if getattr(sys, 'frozen', False):
    base_path = sys._MEIPASS
else:
    base_path = os.path.dirname(os.path.abspath(__file__))

CACHE_DIR = os.path.join(base_path, "cache")


class ComplexEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, (np.float32, np.float64)): return float(obj)
        if isinstance(obj, (np.int32, np.int64)): return int(obj)
        if isinstance(obj, (complex, np.complex64, np.complex128)):
            return {"r": float(obj.real), "i": float(obj.imag)}
        return super().default(obj)

app.json_encoder = ComplexEncoder
def make_serializable(obj):
    return json.loads(json.dumps(obj, cls=ComplexEncoder))


def load_cache(theta):
    theta = round(theta, 2)
    path = os.path.join(CACHE_DIR, f"data_{theta:.2f}.pkl")
    if os.path.exists(path):
        with open(path, "rb") as f:
            return pickle.load(f)
    return None

def init_qubit_logic(theta):
    global current_qubit
    print(f"Initializing Qubit at theta={theta}°...")
    
    data = load_cache(theta)
    
    if data and 'qubit' in data:
        print("  > Loaded from Unified Cache")
        q_data = data['qubit']
        current_qubit = MoireQubit(q_data['E0'], q_data['E1'])
        
        return current_qubit
    
    print("  > Cache missing, cannot initialize.")
    return None

@app.route('/status', methods=['GET'])
def status():
    return jsonify({'status':'online', 'initialized': current_qubit is not None})

@app.route('/initialize', methods=['POST'])
def initialize():
    data = request.json
    try:
        q = init_qubit_logic(float(data.get('theta', 1.0)))
        if q: return jsonify(make_serializable({"success": True, "E0": q.E0, "E1": q.E1}))
        else: return jsonify({"error": "Cache not found for this angle"}), 404
    except Exception as e: return jsonify({"error": str(e)}), 500

@app.route('/bandstructure', methods=['POST'])
def get_bands():
    data = request.json
    theta = round(float(data.get('theta', 3.0)) * 10) / 10
    
    cached = load_cache(theta)
    if cached:
        return jsonify(make_serializable({
            "theta": cached['theta'],
            "bands": cached['bands'],
            "potential": cached['potential']
        }))
        
    return jsonify({"error": "Data not cached"}), 404

@app.route('/evolve', methods=['POST'])
def handle_evolution():
    if not current_qubit: 
        return jsonify({"error": "System initializing..."}), 503
        
    data = request.json
    
    pd = data['psi']
    r1 = pd['c1']['r'] if isinstance(pd['c1'], dict) else float(pd['c1'])
    i1 = pd['c1']['i'] if isinstance(pd['c1'], dict) else 0.0
    r0 = pd['c0']['r'] if isinstance(pd['c0'], dict) else float(pd['c0'])
    i0 = pd['c0']['i'] if isinstance(pd['c0'], dict) else 0.0
    
    c1 = r1 + 1j*i1
    c0 = r0 + 1j*i0
    
    psi0 = c1*basis(2,0) + c0*basis(2,1)
    if psi0.norm() > 1e-9: psi0 = psi0.unit()

    omega = float(data['omega'])
    delta = float(data['delta'])
    duration = float(data['duration'])
    g = float(data.get('g', 0.0))
    t1 = float(data.get('T1', 0))
    T1 = t1 if t1 > 0 else None
    
    t2 = float(data.get('T2', 0))
    T2 = t2 if t2 > 0 else None
    
    pulse_width = data.get('pulse_width', None)
    if pulse_width: pulse_width = float(pulse_width)

    t, res, final = current_qubit.evolve(psi0, omega, delta, duration, pulse_width=pulse_width, g=g, T1=T1, T2=T2)
    
    traj = []
    for i, ti in enumerate(t):
        pt = {
            'time': float(ti),
            'sx': float(np.real(res.expect[0][i])), 
            'sy': float(np.real(res.expect[1][i])), 
            'sz': float(np.real(res.expect[2][i])),
            'population': float(np.real(res.expect[3][i]))
        }
        traj.append(pt)

    if final.type == 'oper':
        pop1 = final[0,0].real
        pop0 = final[1,1].real
        c1_mag = np.sqrt(pop1)
        c0_mag = np.sqrt(pop0)

        coherence = final[0,1]
        phase = np.angle(coherence) if np.abs(coherence) > 1e-9 else 0
        
        c1_val = c1_mag * np.exp(-1j * phase)
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
    global current_qubit
    if not current_qubit: 
        current_qubit = init_qubit_logic(1.0)
        if not current_qubit:
            return jsonify({"error": "Initialization failed (Cache missing)"}), 500

    data = request.json
    
    g = float(data.get('g', 0.5))
    kappa = float(data.get('kappa', 2.0))
    gamma = float(data.get('gamma', 0.1))

    t, atom, photon = current_qubit.simulate_emission(g, kappa, gamma, dur=20)
    stats = current_qubit.get_emission_statistics(g, kappa, gamma)
    
    return jsonify(make_serializable({
        'trajectory': [{'time': t[i], 'atom': atom[i], 'photon': photon[i]} for i in range(len(t))],
        'stats': stats
    }))

if __name__ == '__main__':
    print("--- Moiré Qubit Server Running on Port 5000 ---")
    app.run(port=5000, debug=True)

