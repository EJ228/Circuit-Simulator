import numpy as np
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import os
import copy

# ==============================================================================
#  HELPER: PARSE VALUE
# ==============================================================================
def get_magnitude(val_str):
    try:
        nums = re.findall(r"[-+]?\d*\.\d+|\d+", val_str)
        if nums: return float(nums[0])
    except: pass
    return 0.0

# ==============================================================================
#  HELPER 1: TIE SET MATRIX (B) - FOR LOOP ANALYSIS
# ==============================================================================
def get_tie_set_matrix(num_nodes, branches):
    # 1. SORT BRANCHES
    # Priority: V (0) -> Passive (1) -> I (2)
    sorted_branches = sorted(branches, key=lambda x: 0 if x['type'] == 'V' else (2 if x['type'] == 'I' else 1))
    
    # 2. RE-INDEX BRANCHES (CRITICAL FIX)
    # Update IDs to match the new sorted order so Matrix Columns align
    for i, b in enumerate(sorted_branches):
        b['id'] = i

    # 3. BUILD TREE
    parent = list(range(num_nodes + 1))
    def find(i):
        if parent[i] == i: return i
        parent[i] = find(parent[i])
        return parent[i]
    def union(i, j):
        root_i, root_j = find(i), find(j)
        if root_i != root_j:
            parent[root_i] = root_j
            return True
        return False

    tree_branches = []
    links = []
    
    for br in sorted_branches:
        if union(br['n1'], br['n2']):
            tree_branches.append(br)
        else:
            links.append(br)

    # 4. BUILD B MATRIX
    num_loops = len(links)
    num_branches = len(sorted_branches)
    B = sp.zeros(num_loops, num_branches)
    loop_names = []
    loop_is_fixed = [] 
    fixed_current_vals = []

    adj = {i: [] for i in range(num_nodes + 1)}
    for br in tree_branches:
        adj[br['n1']].append({'to': br['n2'], 'id': br['id'], 'dir': 1})
        adj[br['n2']].append({'to': br['n1'], 'id': br['id'], 'dir': -1})

    def get_path(start, end):
        queue = [(start, [])]
        visited = {start}
        while queue:
            curr, path = queue.pop(0)
            if curr == end: return path
            for neighbor in adj[curr]:
                if neighbor['to'] not in visited:
                    visited.add(neighbor['to'])
                    new_path = path + [{'id': neighbor['id'], 'val': neighbor['dir']}]
                    queue.append((neighbor['to'], new_path))
        return []

    for i, link in enumerate(links):
        loop_names.append(f"Loop_{link['name']}")
        
        # Temporary storage to calculate net voltage drop
        row_vals = {} 
        row_vals[link['id']] = 1 
        
        path = get_path(link['n2'], link['n1'])
        for step in path:
            row_vals[step['id']] = step['val']
            
        # SMART ORIENTATION: Flip loop if it opposes the Voltage Source
        # We want "Clockwise" / Positive current flow out of Source+
        net_drop = 0.0
        
        # Check Link
        mag = get_magnitude(link['val_orig'])
        if link['type'] == 'V': net_drop += mag
        elif link['type'] == 'I': net_drop -= mag * 1e6 # Strong align with I

        # Check Path
        for step in path:
            br = sorted_branches[step['id']]
            mag = get_magnitude(br['val_orig'])
            if br['type'] == 'V':
                net_drop += step['val'] * mag
            elif br['type'] == 'I':
                net_drop -= step['val'] * mag * 1e6

        # If Net Drop is Positive (Opposes Source), Flip the loop
        flip = -1 if net_drop > 0 else 1
        
        for bid, val in row_vals.items():
            B[i, bid] = val * flip
            
        # Set Constraints
        if link['type'] == 'I':
            loop_is_fixed.append(True)
            fixed_current_vals.append(link['Z'] * flip)
        else:
            loop_is_fixed.append(False)
            fixed_current_vals.append(0)

    return B, loop_names, loop_is_fixed, fixed_current_vals, sorted_branches

# ==============================================================================
#  HELPER 2: INCIDENCE MATRIX (A) - FOR NODAL ANALYSIS
# ==============================================================================
def get_incidence_matrix(num_nodes, branches, ground_node=0):
    active_nodes = sorted(list(set([b['n1'] for b in branches] + [b['n2'] for b in branches])))
    if ground_node in active_nodes:
        active_nodes.remove(ground_node)
        
    num_unknowns = len(active_nodes)
    num_branches = len(branches)
    
    A = sp.zeros(num_unknowns, num_branches)
    node_map = {n: i for i, n in enumerate(active_nodes)}
    
    for j, b in enumerate(branches):
        if b['n1'] in node_map: A[node_map[b['n1']], j] = 1 
        if b['n2'] in node_map: A[node_map[b['n2']], j] = -1 
            
    return A, active_nodes

# ==============================================================================
#  MAIN SOLVER
# ==============================================================================
def solve_circuit(netlist_file, tf=1.0):
    # Cleanup
    for f in ['Voltages_graph.png', 'Currents_graph.png', 'Res_C.png', 'Cap_curr.png', 'Results.txt', 'matrices.txt']:
        if os.path.exists(f):
            try: os.remove(f)
            except: pass
    plt.close('all')

    s = sp.symbols('s')
    t = sp.symbols('t', real=True, positive=True)

    # --- 1. PARSE NETLIST ---
    raw_branches = []
    nodes = set()
    with open(netlist_file, 'r') as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        parts = line.strip().split()
        if not parts: continue
        name = parts[0]; type_char = name[0].upper()
        n1 = int(parts[1]); n2 = int(parts[2])
        val_str = " ".join(parts[3:])
        nodes.add(n1); nodes.add(n2)
        
        val_s = 0
        try:
            if 'sine' in val_str:
                nums = re.findall(r"[-+]?\d*\.\d+|\d+", val_str)
                val_s = (float(nums[0]) * float(nums[1])) / (s**2 + float(nums[1])**2)
            elif 'cosine' in val_str:
                nums = re.findall(r"[-+]?\d*\.\d+|\d+", val_str)
                val_s = (float(nums[0]) * s) / (s**2 + float(nums[1])**2)
            else:
                val = float(val_str)
                if type_char == 'R': val_s = val
                elif type_char == 'L': val_s = val * s
                elif type_char == 'C': val_s = 1 / (val * s)
                elif type_char in ['V', 'I']: val_s = val / s
        except: val_s = 0

        raw_branches.append({
            'id': idx, 'name': name, 'type': type_char,
            'n1': n1, 'n2': n2, 'Z': val_s, 'val_orig': val_str
        })

    num_nodes = max(nodes) if nodes else 0

    # =========================================================
    # METHOD A: LOOP ANALYSIS
    # =========================================================
    loop_branches = copy.deepcopy(raw_branches)
    
    # Get B (sorted_loop_branches has Re-indexed IDs)
    B, loop_names, loop_is_fixed, fixed_vals, sorted_loop_branches = get_tie_set_matrix(num_nodes, loop_branches)
    
    pass_idx = [b['id'] for b in sorted_loop_branches if b['type'] in ['R','L','C']]
    src_idx = [b['id'] for b in sorted_loop_branches if b['type'] == 'V']
    
    B_fp = B[:, pass_idx]
    B_fg = B[:, src_idx]
    
    Z_p = sp.zeros(len(pass_idx), len(pass_idx))
    for i, pid in enumerate(pass_idx): Z_p[i,i] = sorted_loop_branches[pid]['Z']
        
    V_g = sp.zeros(len(src_idx), 1)
    for i, sid in enumerate(src_idx): V_g[i,0] = sorted_loop_branches[sid]['Z']

    Z_L = B_fp * Z_p * B_fp.T
    RHS_Loop = -1 * B_fg * V_g 
    
    u_indices = [i for i, fixed in enumerate(loop_is_fixed) if not fixed]
    k_indices = [i for i, fixed in enumerate(loop_is_fixed) if fixed]
    
    I_loop_s = sp.zeros(len(loop_names), 1)
    I_k = sp.Matrix([fixed_vals[i] for i in k_indices])
    for idx in k_indices: I_loop_s[idx, 0] = fixed_vals[idx]
        
    if u_indices:
        Z_uu = Z_L[u_indices, :][:, u_indices]
        Z_uk = Z_L[u_indices, :][:, k_indices]
        RHS_u = RHS_Loop[u_indices, :]
        RHS_eff = RHS_u - Z_uk * I_k if k_indices else RHS_u
        try:
            I_u = Z_uu.inv() * RHS_eff
            for local_i, global_i in enumerate(u_indices):
                I_loop_s[global_i, 0] = I_u[local_i, 0]
        except: pass

    I_branches_s = B.T * I_loop_s

    # =========================================================
    # METHOD B: NODAL ANALYSIS
    # =========================================================
    nodal_branches = copy.deepcopy(raw_branches)
    v_sources = [b for b in nodal_branches if b['type'] == 'V']
    
    for v_src in v_sources:
        for br in nodal_branches:
            if br == v_src: continue
            if br['type'] not in ['R','L','C']: continue
            
            new_I = v_src['Z'] / br['Z']
            i_name = f"I_{v_src['name']}"
            converted = False
            
            # Logic: I points towards V_pos (n1)
            if br['n1'] == v_src['n2']:
                br['n1'] = v_src['n1'] 
                nodal_branches.append({'id':-1, 'name':i_name, 'type':'I', 'n1':br['n2'], 'n2':v_src['n1'], 'Z':new_I})
                converted = True
            elif br['n2'] == v_src['n2']:
                br['n2'] = v_src['n1']
                nodal_branches.append({'id':-1, 'name':i_name, 'type':'I', 'n1':br['n1'], 'n2':v_src['n1'], 'Z':new_I})
                converted = True
            elif br['n1'] == v_src['n1']:
                br['n1'] = v_src['n2']
                nodal_branches.append({'id':-1, 'name':i_name, 'type':'I', 'n1':v_src['n2'], 'n2':br['n2'], 'Z':new_I})
                converted = True
            elif br['n2'] == v_src['n1']:
                br['n2'] = v_src['n2']
                nodal_branches.append({'id':-1, 'name':i_name, 'type':'I', 'n1':v_src['n2'], 'n2':br['n1'], 'Z':new_I})
                converted = True
                
            if converted:
                nodal_branches.remove(v_src)
                break

    nodal_branches = sorted(nodal_branches, key=lambda x: 0 if x['type'] in ['R','L','C'] else 1)
    for i, b in enumerate(nodal_branches): b['id'] = i

    p_idx = [b['id'] for b in nodal_branches if b['type'] in ['R','L','C']]
    g_idx = [b['id'] for b in nodal_branches if b['type'] == 'I']
    
    A, active_nodes = get_incidence_matrix(num_nodes, nodal_branches)
    A_p = A[:, p_idx]; A_g = A[:, g_idx]
    Y_p = sp.zeros(len(p_idx), len(p_idx))
    for i, pid in enumerate(p_idx): Y_p[i,i] = 1 / nodal_branches[pid]['Z']
    I_g = sp.zeros(len(g_idx), 1)
    for i, gid in enumerate(g_idx): I_g[i,0] = nodal_branches[gid]['Z']
        
    Y_n = A_p * Y_p * A_p.T
    RHS_Node = -1 * A_g * I_g
    V_node_s = sp.zeros(len(active_nodes), 1)
    try:
        if Y_n.rows > 0: V_node_s = Y_n.inv() * RHS_Node
    except: pass

    # =========================================================
    # TIME DOMAIN & PLOTTING
    # =========================================================
    t_vals = np.linspace(0, tf, 200)
    
    # Data Processing
    v_node_t = {0: np.zeros_like(t_vals)}
    for i, n in enumerate(active_nodes):
        try:
            ft = sp.lambdify(t, sp.inverse_laplace_transform(V_node_s[i], s, t).replace(sp.Heaviside, lambda x:1), modules='numpy')
            vals = ft(t_vals)
            if np.isscalar(vals): vals = np.full_like(t_vals, vals)
            v_node_t[n] = vals
        except: v_node_t[n] = np.zeros_like(t_vals)

    i_branch_t = {}
    for i, b in enumerate(sorted_loop_branches):
        try:
            ft = sp.lambdify(t, sp.inverse_laplace_transform(I_branches_s[i], s, t).replace(sp.Heaviside, lambda x:1), modules='numpy')
            vals = ft(t_vals)
            if np.isscalar(vals): vals = np.full_like(t_vals, vals)
            i_branch_t[b['name']] = vals
        except: i_branch_t[b['name']] = np.zeros_like(t_vals)

    v_branch_t = {}
    for b in raw_branches:
        v1 = v_node_t.get(b['n1'], np.zeros_like(t_vals))
        v2 = v_node_t.get(b['n2'], np.zeros_like(t_vals))
        v_branch_t[b['name']] = v1 - v2

    i_loop_t = []
    for i in range(len(loop_names)):
        try:
            ft = sp.lambdify(t, sp.inverse_laplace_transform(I_loop_s[i], s, t).replace(sp.Heaviside, lambda x:1), modules='numpy')
            vals = ft(t_vals)
            if np.isscalar(vals): vals = np.full_like(t_vals, vals)
            i_loop_t.append(vals)
        except: i_loop_t.append(np.zeros_like(t_vals))

    # Plots
    plt.figure(figsize=(10,6)); plt.style.use('dark_background')
    for n in active_nodes: plt.plot(t_vals, v_node_t[n], label=f'Node {n}')
    plt.title('Node Voltages'); plt.legend(); plt.savefig('Voltages_graph.png'); plt.close()

    plt.figure(figsize=(10,6)); plt.style.use('dark_background')
    for name, data in i_branch_t.items(): plt.plot(t_vals, data, label=name)
    plt.title('Component Currents'); plt.legend(); plt.savefig('Currents_graph.png'); plt.close()

    plt.figure(figsize=(10,6)); plt.style.use('dark_background')
    for name, data in v_branch_t.items(): plt.plot(t_vals, data, label=name)
    plt.title('Component Voltages'); plt.legend(); plt.savefig('Res_C.png'); plt.close()

    plt.figure(figsize=(10,6)); plt.style.use('dark_background')
    for i, data in enumerate(i_loop_t): plt.plot(t_vals, data, label=loop_names[i])
    plt.title('Loop Currents'); plt.legend(); plt.savefig('Cap_curr.png'); plt.close()

    # --- OUTPUT ---
    def mat_to_str(M, title, row_labels=None, col_labels=None):
        s = f"\n{title} ({M.rows}x{M.cols}):\n" + "-"*40 + "\n"
        if col_labels: s += "Cols: " + ", ".join(col_labels) + "\n"
        for r in range(M.rows):
            r_lab = f"{row_labels[r]} | " if row_labels else ""
            s += r_lab + str(M.row(r)) + "\n"
        return s

    # Labels
    loop_pass_names = [sorted_loop_branches[i]['name'] for i in pass_idx]
    loop_src_names = [sorted_loop_branches[i]['name'] for i in src_idx]
    node_pass_names = [nodal_branches[i]['name'] for i in p_idx]
    node_src_names = [nodal_branches[i]['name'] for i in g_idx]
    node_labels = [f"Node_{n}" for n in active_nodes]

    with open('matrices.txt', 'w') as f:
        f.write(f"MATHEMATICAL ANALYSIS DUMP\nFile: {netlist_file}\n")
        f.write("========================================\n")
        f.write("\n--- LOOP ANALYSIS (Tie Set) ---\n")
        f.write(mat_to_str(B, "Tie Set Matrix (B)", row_labels=loop_names, col_labels=[b['name'] for b in sorted_loop_branches]))
        f.write(mat_to_str(B_fp, "Tie Set - Passive (B_fp)", row_labels=loop_names, col_labels=loop_pass_names))
        f.write(mat_to_str(B_fg, "Tie Set - Sources (B_fg)", row_labels=loop_names, col_labels=loop_src_names))
        f.write(mat_to_str(Z_p, "Impedance Matrix (Z_p)", row_labels=loop_pass_names, col_labels=loop_pass_names))
        f.write(mat_to_str(Z_L, "Loop Impedance (Z_L)"))
        f.write(mat_to_str(I_loop_s, "Loop Currents I_l(s)"))

        f.write("\n\n--- NODAL ANALYSIS (Incidence) ---\n")
        f.write(mat_to_str(A, "Incidence Matrix (A)", row_labels=node_labels, col_labels=[b['name'] for b in nodal_branches]))
        f.write(mat_to_str(A_p, "Incidence - Passive (A_p)", row_labels=node_labels, col_labels=node_pass_names))
        f.write(mat_to_str(A_g, "Incidence - Sources (A_g)", row_labels=node_labels, col_labels=node_src_names))
        f.write(mat_to_str(Y_p, "Admittance Matrix (Y_p)", row_labels=node_pass_names, col_labels=node_pass_names))
        f.write(mat_to_str(Y_n, "Nodal Admittance (Y_n)"))
        f.write(mat_to_str(V_node_s, "Node Voltages V_n(s)"))

    with open('Results.txt', 'w') as f:
        f.write("SIMULATION RESULTS\n------------------\n")
        f.write("Final Node Voltages:\n")
        for n, v in v_node_t.items(): f.write(f"  v_{n} = {v[-1] if len(v)>0 else 0:.4f} V\n")
        f.write("\nFinal Loop Currents:\n")
        for i, l in enumerate(loop_names): f.write(f"  {l} = {i_loop_t[i][-1] if len(i_loop_t)>i else 0:.4f} A\n")
        f.write("\nFinal Component Currents:\n")
        for n, v in i_branch_t.items(): f.write(f"  I_{n} = {v[-1] if len(v)>0 else 0:.4f} A\n")

    return True, "Simulation Complete"