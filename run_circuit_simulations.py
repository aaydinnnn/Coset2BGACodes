import numpy as np
import galois
import json
import sys
import os

GF = galois.GF2
from syndrome_extraction_circuits import get_memory_experiment_circuit 
from stimbposd import SinterDecoder_BPOSD
import sinter
from relay_bp.stim import sinter_decoders


# --- OPTIONAL BEAM SEARCH IMPORT ---
sys.path.append(os.path.join(os.path.dirname(__file__), 'BeamSearchDecoder'))
try:
    from sinter_beamsearch import SinterDecoder_BeamSearch
    BEAM_SEARCH_AVAILABLE = True
except ModuleNotFoundError:
    BEAM_SEARCH_AVAILABLE = False


# ==========================================
# --- INPUT CONFIGURATION ---
# ==========================================

# 1. Choose table: 'Table 2', 'Table 4_Base', 'Table 4_Cover', or 'Table 6'
source_table = 'Table 2'

# 2. Define standard code parameters and GAP identifiers l,m
n = 48
k = 8
d = 6
l = 384
m = 512

# Used for Table 1, Table 3_Cover, and Table 5
s = 53  

# 3. Define a and b
a = [1, 34, 48]
b = [1, 6, 12]

# --- Filename Construction Logic ---
a_str = "a" + "_".join(map(str, a))
b_str = "b" + "_".join(map(str, b))

if source_table in ['Table 2', 'Table 6']:
    prefix = 'code'
    input_filename = f"{prefix}_n{n}_k{k}_d{d}_l{l}_m{m}_s{s}_{a_str}_{b_str}.json"

elif source_table == 'Table 4_Base':
    prefix = 'base'
    input_filename = f"{prefix}_n{n}_k{k}_d{d}_l{l}_m{m}_{a_str}_{b_str}.json"

elif source_table == 'Table 4_Cover':
    prefix = 'cover'
    input_filename = f"{prefix}_n{n}_k{k}_d{d}_l{l}_m{m}_s{s}_{a_str}_{b_str}.json"

else:
    raise ValueError(f"Unknown source_table: {source_table}")

# --- FOLDER PATHING ---
target_folder = 'code_dict'
input_filepath = os.path.join(target_folder, input_filename)

# Auto-extract code name for metadata/saving
code_name = input_filename.replace('.json', '')

# Simulation parameters
basis = 'X'
error_rates = [0.002]
num_cycle = d
num_workers = 4
max_shots = 100000000
max_errors = 50

# ==========================================
# --- MAIN EXECUTION ---
# ==========================================

if __name__ == "__main__":

    print(f"Target File Constructed: {input_filename}")
    
    def parse_matrices_from_json(file_path):
        """
        Parses A and B permutation matrices from a structured JSON file.
        Also returns the metadata dictionary.
        """
        with open(file_path, 'r') as f:
            data = json.load(f)
            
        A_mat_dict = {}
        B_mat_dict = {}
        metadata = data.get('metadata', {}) # <-- Extract metadata
        
        building_blocks = data.get('building_blocks', {})
        
        for key, block in building_blocks.items():
            rows = block['rows']
            cols = block['cols']
            
            size = len(rows)
            mat = np.zeros((size, size), dtype=np.int8)
            
            for r, c in zip(rows, cols):
                mat[r, c] = 1
                
            mat_type = key[0]        # 'A' or 'B'
            idx = int(key[1:])       # 1, 2, 3, etc.
            
            if mat_type == 'A':
                A_mat_dict[idx] = mat
            elif mat_type == 'B':
                B_mat_dict[idx] = mat
                
        return A_mat_dict, B_mat_dict, metadata

    try:
        A_mat_dict, B_mat_dict, metadata = parse_matrices_from_json(input_filepath)
    except FileNotFoundError:
        print(f"\n[ERROR] Could not find the file: {input_filepath}")
        sys.exit(1)

    A_mat = sum(A_mat_dict.values()) % 2
    B_mat = sum(B_mat_dict.values()) % 2

    H_x=np.hstack((A_mat,B_mat))
    H_z=np.hstack((B_mat.T, A_mat.T))

    n_half = A_mat_dict[1].shape[0]      
    n_data = 2 * n_half      
    num_qubits = 2 * n_data 
    q_map = {}
    for i in range(1, n_half + 1):
        q_map[('X', i)] = i - 1
    for i in range(1, n_half + 1):
        q_map[('L', i)] = n_half + (i - 1)
    for i in range(1, n_half + 1):
        q_map[('R', i)] = 2 * n_half + (i - 1)
    for i in range(1, n_half + 1):
        q_map[('Z', i)] = 3 * n_half + (i - 1)

    # --- CODE SANITY CHECK PRINT ---
    meta_n = metadata.get('n', '?')
    meta_k = metadata.get('k', '?')
    meta_d = metadata.get('d', '?')
    exact = metadata.get('is_exact_distance', True)
    
    d_str = f"{meta_d}" if exact else f"<={meta_d}"
    
    # 1. Calculate the rank of the check matrices over GF(2)
    rank_X = np.linalg.matrix_rank(GF(H_x))
    rank_Z = np.linalg.matrix_rank(GF(H_z))
    
    # 2. k = n - rank(H_X) - rank(H_Z)
    calculated_k = H_x.shape[1] - rank_X - rank_Z
    
    # 3. Check for a match
    k_match = "MATCH" if calculated_k == meta_k else f"MISMATCH (Calculated: {calculated_k})"

    print(f"--- Code Sanity Check ---")
    print(f"Code parameters: [[{meta_n}, {meta_k}, {d_str}]]")
    print(f"data qubits = {H_x.shape[1]}")
    print(f"X checks = {H_x.shape[0]} (Rank: {rank_X})")
    print(f"Z checks = {H_z.shape[0]} (Rank: {rank_Z})")
    print(f"logical qubits (meta) = {meta_k} [{k_match}]")
    print(f"-------------------------")


    bposd_decoder = SinterDecoder_BPOSD(
        max_bp_iters=10000,     # 10,000 iterations
        bp_method="ms",         # Min-Sum
        osd_order=10,           # Combination Sweep order
        osd_method="osd_cs"     # Combination Sweep
    )

    relay_bp_decoder = sinter_decoders(
    gamma0=0.125,
    pre_iter=80,
    num_sets=60,
    set_max_iter=60,
    gamma_dist_interval=(-0.24, 0.66),
    stop_nconv=1,
    )['relay-bp']

    # Initialize the dictionary with standard decoders
    custom_decoders_dict = {
        'bposd': bposd_decoder, 
        'relay-bp': relay_bp_decoder
    }

    # Only add Beam Search if it was successfully imported
    if BEAM_SEARCH_AVAILABLE:
        custom_decoders_dict['beam32_340iters'] = SinterDecoder_BeamSearch(
            max_rounds=10,
            beam_width=32, 
            initial_iters=40, 
            iters_per_round=30,
            num_results=1
        )

    checkpoint_file = f"checkpoint_{code_name}_{basis}.csv"
    final_file = f"results_{code_name}_{basis}.csv"

    print(f"Generating tasks for {code_name} (Basis: {basis})...")

    tasks = []
    for p in error_rates:
        tasks.append(
            sinter.Task(
                circuit=get_memory_experiment_circuit(d=num_cycle,
                                                    basis=basis,
                                                    A_mat_dict=A_mat_dict,
                                                    B_mat_dict=B_mat_dict,
                                                    q_map=q_map,
                                                    p=p),
                json_metadata={'Code Name': code_name, 'weight': len(A_mat_dict)+len(B_mat_dict), 'basis': basis, 'p': p, 'NumberofCycle': num_cycle}
            )
        )

    print(f"Starting Sinter simulation for {code_name}...")
    print(f"-> Checkpoints will be saved to: {checkpoint_file}")

    # CHOOSE DECODER:
    # Options: 'bposd', 'relay-bp', 'beam32_340iters'
    active_decoders = ['bposd']
    results = sinter.collect(
        num_workers=num_workers,
        max_shots=max_shots,
        max_errors=max_errors,
        tasks=tasks,
        decoders=active_decoders,
        custom_decoders=custom_decoders_dict,
        print_progress=True,
        save_resume_filepath=checkpoint_file 
    )

# --- Save & Print Results ---
    print("\n--- Simulation Complete ---")
    print(f"Saving final results to '{final_file}'...")
    
    with open(final_file, "w") as f: # <-- Uses dynamic name
        f.write(sinter.CSV_HEADER + "\n")
        for sample in results:
            f.write(sample.to_csv_line() + "\n")
            
    print("\nResults Summary:")
    print(sinter.CSV_HEADER)
    for sample in results:
        print(sample.to_csv_line())