import json
from sage.interfaces.gap import gap
from sage.all import GF, matrix, Matrix
import numpy as np

# ==========================================
# 1. PARAMETERS & CONFIGURATION
# ==========================================
l = 384
m = 512
s = 27
a = [1, 9, 87]
b = [1, 21, 23]

save_output = True # Set to True to generate the .json file

num_trials = 1000 # Number of trials for QDistRnd Algorithm
min_d_target = 0

print("===== PARAMETERS LOADED =====")
print(f"l={l}, m={m}, s={s}")
print(f"a={a}, type: {type(a)}")
print(f"b={b}, type: {type(b)}\n")

# ==========================================
# 2. GAP GROUP SETUP & SUBGROUP ANALYSIS
# ==========================================
gap.eval(f"""
G := SmallGroup({l}, {m});
desc := StructureDescription(G);
subs := AllSubgroups(G);;
non_normal := Filtered(subs, H -> not IsNormal(G,H));;
""")

G = gap("G")
G_order = G.Size()
desc = gap("desc")

subs_list = gap("subs")
non_normal_list = gap("non_normal")

all_subgroups = [gap(H) for H in subs_list]
non_normal_groups = [gap(H) for H in non_normal_list]

print("===== GROUP G INFO =====")
print(f"Description of G: {desc}")
print(f"Order of G: {G_order}")
print(f"Total subgroups: {len(all_subgroups)}")
print(f"Non-normal subgroups: {len(non_normal_groups)}")

gap.eval(f"H := non_normal[{s}];")
order_H = gap.eval("Size(H)")
gap.eval("descH := StructureDescription(H)")
descH = gap("descH")

print("\n===== SUBGROUP H INFO =====")
print(f"Description of H: {descH}")
print(f"Order of H: {order_H}")
print(f"Label of H: non-normal, index {s}\n")

# Cosets G/H (Action Space)
gap.eval(f"Cs := LeftCosets(G, H);")
num_leftCosets = int(gap.eval("Length(Cs);")) 
print(f"[G:H] = {num_leftCosets}") 

# Compute Normalizer N and Core
gap.eval(f"N := Normalizer(G, H);")
gap.eval("Cs_NH_H := LeftCosets(N,H);")
num_N_modH = int(gap.eval("Length(Cs_NH_H);")) 
print(f"Index |N/H| = {num_N_modH}") 

gap.eval(f"Core_H := Core(G,H);")
gap.eval(f"Cs_G_coreH := LeftCosets(G,Core_H);")
num_G_CoreH = int(gap.eval("Length(Cs_G_coreH);")) 
print(f"Index |G/Core(H)| = {num_G_CoreH}") 

print("-" * 40)

# ==========================================
# 3. MATRIX PRE-CALCULATION PHASE
# ==========================================
# This step pre-calculates all unique left and right action matrices so they do not need to be repeatedly built while searching for codes.
F2 = GF(2)

def calculate_permutation_matrix(local_gap, coset_list_G_H, coset_list_G_X, coset_index_G_X, is_left_action):
    local_gap.eval(f"g := Elements({coset_list_G_X}[{coset_index_G_X}])[1];")
    n_cosets = int(local_gap.eval(f"Length({coset_list_G_H});")) 
    sigma = []

    for i in range(1, n_cosets + 1):
        local_gap.eval(f"rep := Elements({coset_list_G_H}[{i}])[1];")
        if is_left_action:
            local_gap.eval("prod := g * rep;")
        else:
            local_gap.eval("prod := rep * g;")
        
        coset_index = int(local_gap.eval("""
        i := 1;;
        while not prod in Cs[i] do
            i := i + 1;;
        od;;
        i;
        """)) - 1 
        sigma.append(coset_index)

    P = matrix(F2, n_cosets, n_cosets)
    for i, j in enumerate(sigma):
        P[j, i] = 1
    return P

print(f"Pre-calculating {num_G_CoreH} L matrices (A_i)...")
L_mat_lookup = {i: calculate_permutation_matrix(gap, "Cs", "Cs_G_coreH", i, True) for i in range(1, num_G_CoreH + 1)}
print("L matrices complete.")

print(f"Pre-calculating {num_N_modH} R matrices (B_j)...")
R_mat_lookup = {j: calculate_permutation_matrix(gap, "Cs", "Cs_NH_H", j, False) for j in range(1, num_N_modH + 1)}
print("R matrices complete.")
print("-" * 40)

# ==========================================
# 4. CONSTRUCT PARITY CHECK MATRICES
# ==========================================
def A_mat(a_list):
    M = matrix(F2, num_leftCosets, num_leftCosets) 
    for k in a_list:
        M += L_mat_lookup[k]
    return M

def B_mat(b_list):
    M = matrix(F2, num_leftCosets, num_leftCosets) 
    for k in b_list:
        M += R_mat_lookup[k]
    return M 

def get_nonzero_coords(mat):
    mat_np = np.array(mat, dtype=int)
    rows, cols = np.nonzero(mat_np)
    return {"rows": rows.tolist(), "cols": cols.tolist()}

print("Loading QDistRnd Package...")
gap.eval('LoadPackage("QDistRnd");')
print("Package loaded.\n")
print("-" * 40)

A = A_mat(a)
B = B_mat(b)

H_x = np.hstack((A.numpy(), B.numpy()))
H_z = np.vstack((B.numpy(), -A.numpy())).T

Hx = Matrix(F2, H_x)
Hz = Matrix(F2, H_z)

Hx_gap = Hx._gap_()
Hz_gap = Hz._gap_()

n_code = Hx.ncols()
k_code = n_code - Hx.rank() - Hz.rank()

print(f"Testing Code [n={n_code}, k={k_code}]")


d_z = int(gap.eval(f"DistRandCSS({Hx_gap}, {Hz_gap}, {num_trials}, {min_d_target});"))
d_x = int(gap.eval(f"DistRandCSS({Hz_gap}, {Hx_gap}, {num_trials}, {min_d_target});"))
d_ub = min(d_x, d_z)
print(f"d = {d_ub} (d_x={d_x}, d_z={d_z})")

# ==========================================
# 5. JSON EXPORT
# ==========================================
if save_output:
    a_str = "_".join(map(str, a))
    b_str = "_".join(map(str, b))

    filename = f"code_n{n_code}_k{k_code}_d{d_ub}_l{l}_m{m}_s{s}_a{a_str}_b{b_str}.json"

    export_data = {
        "metadata": {
            "n": int(n_code),
            "k": int(k_code),
            "d": int(d_ub),
            "is_exact_distance": False, 
            "l": int(l),
            "m": int(m),
            "s": int(s),
            "a": a,
            "b": b
        },
        "building_blocks": {}
    }

    for i, idx_a in enumerate(a, start=1):
        export_data["building_blocks"][f"A{i}"] = get_nonzero_coords(L_mat_lookup[idx_a])
    for j, idx_b in enumerate(b, start=1):
        export_data["building_blocks"][f"B{j}"] = get_nonzero_coords(R_mat_lookup[idx_b])
        
    with open(filename, "w") as f:
        json.dump(export_data, f, indent=2)
        
    print(f"\nSaved: {filename}")