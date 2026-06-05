#  Circuit Simulator

A full-stack web app that solves linear **RLC circuits symbolically in the Laplace (s) domain** — not just numerically — using two independent network-graph-theory methods and reconstructing the time-domain response via inverse Laplace transforms.

**Stack:** Python · Flask · SymPy · NumPy · Matplotlib · React · Tailwind CSS

---

## Overview

You describe a circuit as a SPICE-style netlist (resistors, inductors, capacitors, and independent voltage/current sources) in a React frontend, choose a source waveform and a simulation time, and run it. A Flask REST API passes the netlist to a Python engine that:

1. Converts every element into its s-domain impedance/admittance.
2. Solves the circuit **two ways** — loop analysis (tie-set matrix) and nodal analysis (incidence matrix) — straight from network graph theory.
3. Inverts the s-domain solution back to the time domain and plots node voltages and branch/loop currents.

It also writes a full symbolic dump of every intermediate matrix, so the math is transparent rather than a black box.

## Features

- **Symbolic Laplace-domain solving** — exact transfer functions in `s`, not finite-difference approximations.
- **Dual analysis methods** — fundamental-loop (tie-set) *and* nodal (incidence) formulations, computed independently.
- **Graph-theoretic core** — spanning tree built with a union-find (disjoint-set) algorithm; fundamental loops found via BFS.
- **Source handling** — DC step, sine, and cosine sources; automatic source transformation (Thévenin → Norton) for the nodal solver; current-source loop constraints.
- **Time-domain reconstruction** — symbolic inverse Laplace transform, then `lambdify` to NumPy for fast evaluation.
- **Visual + traceable output** — Matplotlib plots plus a `matrices.txt` dump of the tie-set/incidence, impedance/admittance, and s-domain solution matrices.
- **Simple REST API** — clean JSON endpoints, CORS-enabled, with the UI served by the same Flask app.

## How It Works

**1. s-domain modeling.** Each element is mapped to its Laplace-domain form: `R → R`, `L → Ls`, `C → 1/(Cs)`. A DC source of value `V` becomes `V/s`; a sine source `A·sin(ωt)` becomes `Aω/(s² + ω²)`; a cosine source `A·cos(ωt)` becomes `As/(s² + ω²)`.

**2. Loop analysis (tie-set matrix B).** Branches are sorted (sources → passives → current sources), a spanning tree is built with union-find, and each link defines one fundamental loop (its path through the tree is found by BFS). This yields the tie-set matrix **B**, partitioned into passive (`B_fp`) and source (`B_fg`) columns. The loop-impedance system is formed as `Z_L = B_fp · Z_p · B_fpᵀ` and solved for the loop currents `I_l(s)`; current sources are imposed as known loop currents and the remaining unknowns solved by submatrix reduction.

**3. Nodal analysis (incidence matrix A).** Voltage sources in series with passives are converted to equivalent current sources, then a reduced incidence matrix **A** (ground = node 0) is partitioned into passive (`A_p`) and source (`A_g`) columns. The nodal-admittance system `Y_n = A_p · Y_p · A_pᵀ` is solved for the node voltages `V_n(s)`.

**4. Time domain.** The s-domain node voltages and branch/loop currents are run through SymPy's `inverse_laplace_transform`, evaluated over a time grid, and plotted.

## Tech Stack

| Layer | Technology |
| --- | --- |
| Frontend | React 18, Tailwind CSS (single-page `index.html`) |
| Backend | Flask + Flask-CORS (REST API) |
| Math engine | SymPy (symbolic algebra, Laplace transforms), NumPy |
| Plotting | Matplotlib (Agg backend) |

## Project Structure

```
Circuit Simulator/
├── app.py                 # Flask server + REST API
├── circuit_solver.py      # Core symbolic solver (tie-set + incidence analysis)
├── index.html             # React + Tailwind single-page frontend
├── requirements.txt       # Python dependencies (see below)
├── .gitignore
└── README.md
```

> Generated at runtime (and git-ignored): `circuit_netlist.txt`, `Results.txt`, `matrices.txt`, and the four output PNGs.

## Getting Started

### Prerequisites
- Python 3.10 or newer

### Installation

```bash
# 1. Clone the repository
git clone <your-repo-url>
cd "Circuit Simulator"

# 2. Create and activate a virtual environment
python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS / Linux
source .venv/bin/activate

# 3. Install dependencies
pip install flask flask-cors sympy numpy matplotlib
# or, if you keep a requirements file:
pip install -r requirements.txt
```

A minimal `requirements.txt`:

```
flask
flask-cors
sympy
numpy
matplotlib
```

### Running

```bash
python app.py
```

Then open **http://localhost:5000** in your browser.

## Netlist Format

Each line describes one element: `<Name> <node1> <node2> <value>`. The **first letter of the name** sets the element type; the rest is just a label. **Node 0 is ground**, and for sources the **first node listed is the positive terminal**.

| Element | Syntax | Example | Meaning |
| --- | --- | --- | --- |
| Resistor | `R… n1 n2 <ohms>` | `R1 2 1 10` | 10 Ω between nodes 2 and 1 |
| Inductor | `L… n1 n2 <henries>` | `L1 1 0 50` | 50 H between nodes 1 and 0 |
| Capacitor | `C… n1 n2 <farads>` | `C1 1 0 1e-6` | 1 µF between nodes 1 and 0 |
| DC voltage source | `V… n+ n- <volts>` | `V1 2 0 5` | 5 V step, + at node 2 |
| Sine voltage source | `V… n+ n- sine <amp> <ω>` | `V2 2 0 sine 10 100` | 10·sin(100·t) V |
| Cosine voltage source | `V… n+ n- cosine <amp> <ω>` | `V3 2 0 cosine 10 100` | 10·cos(100·t) V |
| Current source | `I… n1 n2 <amps>` | `I1 1 0 2` | 2 A from node 1 to 0 |

> `ω` is **angular frequency in rad/s** (not Hz).

**Example** — a series R–L circuit driven by a 10·sin(100t) V source:

```
R1 2 1 10
L1 1 0 50
V2 2 0 sine 10 100
```

## API Reference

| Method | Endpoint | Description |
| --- | --- | --- |
| `GET` | `/` | Serves the web UI |
| `GET` | `/api/health` | Backend health check |
| `POST` | `/api/simulate` | Runs a simulation. JSON body: `{ "netlist": [...], "tf": <final_time> }` |
| `GET` | `/api/results` | Returns the latest results as text |
| `GET` | `/api/figures/<filename>` | Returns a generated plot PNG |

## Output Files

| File | Contents |
| --- | --- |
| `Voltages_graph.png` | Node voltages vs. time |
| `Currents_graph.png` | Component (branch) currents vs. time |
| `Res_C.png` | Component voltages vs. time |
| `Cap_curr.png` | Loop currents vs. time |
| `Results.txt` | Final node voltages, loop currents, and branch currents |
| `matrices.txt` | Full symbolic dump: tie-set/incidence matrices, impedance/admittance matrices, and s-domain solutions |

## Demo / Screenshots

To showcase results in this README, drop a few example plots into a `screenshots/` folder (this folder is **not** git-ignored, unlike the live outputs) and embed them, e.g.:

```markdown
![Node voltages](screenshots/voltages.png)
```

## Limitations

- Handles **linear, time-invariant** circuits with **independent** sources (R, L, C, V, I).
- Source waveforms supported: DC step, sine, cosine.
- No dependent/controlled sources, nonlinear elements (diodes, transistors), or mutual inductance.
- The solver reads/writes fixed filenames in the working directory, so it processes one simulation at a time.

## Possible Improvements

- Dependent sources and op-amp models.
- Frequency-response (Bode) plots from the s-domain transfer functions.
- Per-request working directories to support concurrent simulations.
- Export of results to CSV/JSON.

## License

Released under the MIT License. Add a `LICENSE` file if you intend to open-source it.
