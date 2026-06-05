from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
import os
# Import our new python solver
import circuit_solver 

app = Flask(__name__, static_folder='.', static_url_path='')
CORS(app)

# ---------------- CONFIGURE ----------------
NETLIST_FILE = "circuit_netlist.txt"
RESULTS_FILE = "Results.txt"
FIGURE_FILES = ['Voltages_graph.png', 'Currents_graph.png', 'Cap_curr.png', 'Res_C.png']

@app.route('/')
def index():
    return send_file('index.html')

@app.route('/api/health', methods=['GET'])
def health_check():
    return jsonify({
        "status": "online",
        "message": "Python-Sim2 Backend is running",
        "matlab_available": False # We don't need MATLAB anymore
    })

@app.route('/api/simulate', methods=['POST'])
def simulate_circuit():
    try:
        data = request.json
        netlist = data.get('netlist', [])
        tf = float(data.get('tf', 1.0))

        if not netlist:
            return jsonify({"error": "No netlist provided"}), 400

        # 1. Write netlist file
        with open(NETLIST_FILE, 'w') as f:
            for line in netlist:
                # Ensure line components are strings
                netlist_line = ' '.join(str(x) for x in line)
                f.write(netlist_line + '\n')

        # 2. Run the Python Solver (Replacing MATLAB)
        # This function generates the PNGs and Results.txt directly
        success, msg = circuit_solver.solve_circuit(NETLIST_FILE, tf)

        if not success:
             return jsonify({
                 "error": "Simulation failed",
                 "matlab_output": msg
             }), 500

        # 3. Read results
        results = ""
        if os.path.exists(RESULTS_FILE):
            with open(RESULTS_FILE, 'r') as f:
                results = f.read()

        # 4. Check for generated images
        figures = [fig for fig in FIGURE_FILES if os.path.exists(fig)]

        return jsonify({
            "success": True,
            "results": results,
            "figures": figures,
            "matlab_output": msg,
            "matlab_error": ""
        })

    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('/api/results', methods=['GET'])
def get_results():
    if not os.path.exists(RESULTS_FILE):
        return jsonify({"error": "No results available"}), 404
    with open(RESULTS_FILE, 'r') as f:
        return jsonify({"results": f.read()})

@app.route('/api/figures/<filename>', methods=['GET'])
def get_figure(filename):
    if filename not in FIGURE_FILES:
        return jsonify({"error": "File not allowed"}), 403
    if os.path.exists(filename):
        return send_file(filename, mimetype='image/png')
    else:
        return jsonify({"error": "File not found"}), 404

    
if __name__ == '__main__':
    print("\n" + "="*60)
    print("Circuit Simulator")
    print("Running with sim2.py calculation logic + Matplotlib")
    print("="*60)
    app.run(debug=True, port=5000, host='0.0.0.0')