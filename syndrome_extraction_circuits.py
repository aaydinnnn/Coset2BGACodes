from typing import Dict, Tuple, Any
import numpy as np
import stim
import qldpc


def FindQubit(P: np.ndarray, j: int) -> int:
    # Get the row j-1
    row = P[j-1, :]
    indices = np.nonzero(row)[0]
    
    if len(indices) == 0:
        raise ValueError(f"Row {j-1} is empty. Matrix P must be a valid permutation matrix.")
        
    return int(indices[0])+1

def get_one_cycle_circuit(A_mat_dict: Dict[int, np.ndarray], B_mat_dict: Dict[int, np.ndarray], q_map: Dict[Tuple[str, int], int]) -> stim.Circuit:
    """
    Generates exactly ONE cycle for weight 6 codes (Rounds 1-8) of the syndrome extraction.
    Does NOT include the initial zeroing of Z-ancillas (that belongs in the outer loop).
    """
    circuit = stim.Circuit()
    n_half = A_mat_dict[1].shape[0]
    x_ancillas = [q_map[('X', i)] for i in range(1, n_half + 1)]
    z_ancillas = [q_map[('Z', i)] for i in range(1, n_half + 1)]

    # --- Round 1 (Init X / Interact) ---
    circuit.append("TICK")
    circuit.append("RX", x_ancillas)
    targets = []
    mat = A_mat_dict[1].T
    for i in range(1, n_half + 1):
        targets.extend([q_map[('R', FindQubit(mat, i))], q_map[('Z', i)]])
    circuit.append("CNOT", targets)
    # --- Round 2 ---
    circuit.append("TICK")
    targets = []
    mat_X, mat_Z = A_mat_dict[2], A_mat_dict[3].T
    for i in range(1, n_half + 1):
        targets.extend([q_map[('X', i)], q_map[('L', FindQubit(mat_X, i))]])
        targets.extend([q_map[('R', FindQubit(mat_Z, i))], q_map[('Z', i)]])
    circuit.append("CNOT", targets)

    # --- Round 3 ---
    circuit.append("TICK")
    targets = []
    mat_X, mat_Z = B_mat_dict[2], B_mat_dict[1].T
    for i in range(1, n_half + 1):
        targets.extend([q_map[('X', i)], q_map[('R', FindQubit(mat_X, i))]])
        targets.extend([q_map[('L', FindQubit(mat_Z, i))], q_map[('Z', i)]])
    circuit.append("CNOT", targets)

    # --- Round 4 ---
    circuit.append("TICK")
    targets = []
    mat_X, mat_Z = B_mat_dict[1], B_mat_dict[2].T
    for i in range(1, n_half + 1):
        targets.extend([q_map[('X', i)], q_map[('R', FindQubit(mat_X, i))]])
        targets.extend([q_map[('L', FindQubit(mat_Z, i))], q_map[('Z', i)]])
    circuit.append("CNOT", targets)

    # --- Round 5 ---
    circuit.append("TICK")
    targets = []
    mat_X, mat_Z = B_mat_dict[3], B_mat_dict[3].T
    for i in range(1, n_half + 1):
        targets.extend([q_map[('X', i)], q_map[('R', FindQubit(mat_X, i))]])
        targets.extend([q_map[('L', FindQubit(mat_Z, i))], q_map[('Z', i)]])
    circuit.append("CNOT", targets)

    # --- Round 6 ---
    circuit.append("TICK")
    targets = []
    mat_X, mat_Z = A_mat_dict[1], A_mat_dict[2].T
    for i in range(1, n_half + 1):
        targets.extend([q_map[('X', i)], q_map[('L', FindQubit(mat_X, i))]])
        targets.extend([q_map[('R', FindQubit(mat_Z, i))], q_map[('Z', i)]])
    circuit.append("CNOT", targets)

    # --- Round 7 (Interact / Measure Z) ---
    circuit.append("TICK")
    targets = []
    mat_X = A_mat_dict[3]
    for i in range(1, n_half + 1):
        targets.extend([q_map[('X', i)], q_map[('L', FindQubit(mat_X, i))]])
    circuit.append("CNOT", targets)
    #circuit.append("TICK") 
    circuit.append("M", z_ancillas)

    # --- Round 8 (Measure X) ---
    circuit.append("TICK")
    circuit.append("MX", x_ancillas)
    circuit.append("R", z_ancillas) # Prepare Z ancillas for next cycle (or end)

    return circuit

def get_one_cycle_circuit_weight_8(A_mat_dict: Dict[int, np.ndarray], B_mat_dict: Dict[int, np.ndarray], q_map: Dict[Tuple[str, int], int]) -> stim.Circuit:
    """
    Generates exactly ONE cycle of the syndrome extraction for weight-8 codes.
    Fully synchronous: Prepares all ancillas at the start, measures all at the end.
    """
    circuit = stim.Circuit()
    n_half = A_mat_dict[1].shape[0]
    x_ancillas = [q_map[('X', i)] for i in range(1, n_half + 1)]
    z_ancillas = [q_map[('Z', i)] for i in range(1, n_half + 1)]

    # --- Initialization (Synchronous Prep) ---
    circuit.append("RX", x_ancillas)
    circuit.append("R", z_ancillas)

    # --- Round 1: CNOT_A2(1, 2) and CNOT_A1(3, 4) ---
    circuit.append("TICK")
    targets = []
    mat_X, mat_Z = A_mat_dict[2], A_mat_dict[1].T
    for i in range(1, n_half + 1):
        targets.extend([q_map[('X', i)], q_map[('L', FindQubit(mat_X, i))]])
        targets.extend([q_map[('R', FindQubit(mat_Z, i))], q_map[('Z', i)]])
    circuit.append("CNOT", targets)

    # --- Round 2: CNOT_A3(1, 2) and CNOT_A4(3, 4) ---
    circuit.append("TICK")
    targets = []
    mat_X, mat_Z = A_mat_dict[3], A_mat_dict[4].T
    for i in range(1, n_half + 1):
        targets.extend([q_map[('X', i)], q_map[('L', FindQubit(mat_X, i))]])
        targets.extend([q_map[('R', FindQubit(mat_Z, i))], q_map[('Z', i)]])
    circuit.append("CNOT", targets)

    # --- Round 3: CNOT_B1(1, 3) and CNOT_B2(2, 4) ---
    circuit.append("TICK")
    targets = []
    mat_X, mat_Z = B_mat_dict[1], B_mat_dict[2].T
    for i in range(1, n_half + 1):
        targets.extend([q_map[('X', i)], q_map[('R', FindQubit(mat_X, i))]])
        targets.extend([q_map[('L', FindQubit(mat_Z, i))], q_map[('Z', i)]])
    circuit.append("CNOT", targets)

    # --- Round 4: CNOT_B2(1, 3) and CNOT_B1(2, 4) ---
    circuit.append("TICK")
    targets = []
    mat_X, mat_Z = B_mat_dict[2], B_mat_dict[1].T
    for i in range(1, n_half + 1):
        targets.extend([q_map[('X', i)], q_map[('R', FindQubit(mat_X, i))]])
        targets.extend([q_map[('L', FindQubit(mat_Z, i))], q_map[('Z', i)]])
    circuit.append("CNOT", targets)

    # --- Round 5: CNOT_B3(1, 3) and CNOT_B4(2, 4) ---
    circuit.append("TICK")
    targets = []
    mat_X, mat_Z = B_mat_dict[3], B_mat_dict[4].T
    for i in range(1, n_half + 1):
        targets.extend([q_map[('X', i)], q_map[('R', FindQubit(mat_X, i))]])
        targets.extend([q_map[('L', FindQubit(mat_Z, i))], q_map[('Z', i)]])
    circuit.append("CNOT", targets)

    # --- Round 6: CNOT_B4(1, 3) and CNOT_B3(2, 4) ---
    circuit.append("TICK")
    targets = []
    mat_X, mat_Z = B_mat_dict[4], B_mat_dict[3].T
    for i in range(1, n_half + 1):
        targets.extend([q_map[('X', i)], q_map[('R', FindQubit(mat_X, i))]])
        targets.extend([q_map[('L', FindQubit(mat_Z, i))], q_map[('Z', i)]])
    circuit.append("CNOT", targets)

    # --- Round 7: CNOT_A1(1, 2) and CNOT_A2(3, 4) ---
    circuit.append("TICK")
    targets = []
    mat_X, mat_Z = A_mat_dict[1], A_mat_dict[2].T
    for i in range(1, n_half + 1):
        targets.extend([q_map[('X', i)], q_map[('L', FindQubit(mat_X, i))]])
        targets.extend([q_map[('R', FindQubit(mat_Z, i))], q_map[('Z', i)]])
    circuit.append("CNOT", targets)

    # --- Round 8: CNOT_A4(1, 2) and CNOT_A3(3, 4) ---
    circuit.append("TICK")
    targets = []
    mat_X, mat_Z = A_mat_dict[4], A_mat_dict[3].T
    for i in range(1, n_half + 1):
        targets.extend([q_map[('X', i)], q_map[('L', FindQubit(mat_X, i))]])
        targets.extend([q_map[('R', FindQubit(mat_Z, i))], q_map[('Z', i)]])
    circuit.append("CNOT", targets)

    # --- Measurements (Synchronous Readout) ---
    circuit.append("TICK")
    circuit.append("M", z_ancillas)
    circuit.append("MX", x_ancillas)

    circuit.append("TICK")


    return circuit

def get_noisy_circuit(noiseless_circuit: stim.Circuit, p: float) -> stim.Circuit:
    """
    Applies gate noise, idle noise, and state preparation/measurement noise.
    
    NOTES:
    1. Measurements (M): Noise is applied BEFORE to flip the result.
       - Uses X_ERROR(p) for Z-basis M so the flip probability is exactly p.
    2. Resets/Gates: Noise is applied AFTER to scramble the resulting state.
    3. Idle Noise: Applied to any qubit not involved in an operation during a TICK.
    """
    if p == 0:
        return noiseless_circuit
        
    noisy_circuit = stim.Circuit()
    
    # 1. Define the Universe of Qubits
    all_qubits_set = set()
    for instruction in noiseless_circuit:
        for t in instruction.targets_copy():
            if t.is_qubit_target:
                all_qubits_set.add(t.value)
                
    # 2. Inject Noise
    busy_qubits = set()
    first_tick = True
    
    for instruction in noiseless_circuit:
        
        # --- A. Handle Time Steps (Idle Noise) ---
        if instruction.name == "TICK":
            # Only skip if the circuit literally started with an empty TICK
            if first_tick and len(busy_qubits) == 0:
                first_tick = False
                noisy_circuit.append("TICK")
                continue
            
            first_tick = False # Make sure we turn the flag off

            # Apply idle noise to any qubit that was NOT busy in the previous moment
            idle_qubits = list(all_qubits_set - busy_qubits)
            if idle_qubits:
                noisy_circuit.append("DEPOLARIZE1", idle_qubits, p)
            
            busy_qubits = set()
            noisy_circuit.append("TICK")
            continue

        # --- B. Handle Instructions ---
        targets = instruction.targets_copy()
        
        # Track busy qubits for the current moment
        for t in targets:
            if t.is_qubit_target:
                busy_qubits.add(t.value)
        
        # CASE 1: Measurements (Noise BEFORE)
        # We must flip the qubit *before* Stim records the measurement result.
        if instruction.name == "M":
            # X_ERROR flips Z-basis eigenstates (|0> <-> |1>)
            noisy_circuit.append("X_ERROR", targets, p)
            noisy_circuit.append(instruction)
            
        elif instruction.name == "MX":
            # Z_ERROR flips X-basis eigenstates (|+> <-> |->)
            noisy_circuit.append("Z_ERROR", targets, p)
            noisy_circuit.append(instruction)
            
        # CASE 2: Resets & Gates (Noise AFTER)
        else:
            noisy_circuit.append(instruction)
            
            # Apply noise immediately after the operation finishes
            if instruction.name in ["R", "RZ"]:
                 # Reset to |0>. Noise is X (bit flip).
                 noisy_circuit.append("X_ERROR", targets, p)

            elif instruction.name == "RX":
                 # Reset to |+>. Noise is Z (phase flip).
                 noisy_circuit.append("Z_ERROR", targets, p)
                
            elif instruction.name in ["H", "S", "X", "Y", "Z", "SQRT_X", "SQRT_Y", "I"]:
                # Single-qubit gate error
                noisy_circuit.append("DEPOLARIZE1", targets, p)
                
            elif instruction.name in ["CNOT", "CX", "CY", "CZ", "SWAP"]:
                # Two-qubit gate error
                noisy_circuit.append("DEPOLARIZE2", targets, p)

# --- C. Handle the Final Round Idle ---
    # Only apply if there are leftover operations that weren't closed by a TICK
    if not first_tick and len(busy_qubits) > 0: 
        idle_qubits = list(all_qubits_set - busy_qubits)
        if idle_qubits:
            noisy_circuit.append("DEPOLARIZE1", idle_qubits, p)

    return noisy_circuit

def get_data_init_block(basis: str, q_map: Dict[Tuple[str, int], int], p: float = 0.0) -> stim.Circuit:
    circuit = stim.Circuit()
    n_half = int(len(q_map)/4)
    
    # Gather Data Qubits in order
    l_data = [q_map[('L', i)] for i in range(1, n_half + 1)]
    r_data = [q_map[('R', i)] for i in range(1, n_half + 1)]
    all_data = l_data + r_data
    
    # 1. Initialize in the correct basis
    if basis == 'X':
        circuit.append("RX", all_data)
        # Text: "qubits initialised in |+⟩ experience single-qubit Z errors"
        if p > 0:
            circuit.append("Z_ERROR", all_data, p)
    else: # basis == 'Z'
        circuit.append("R", all_data)
        # Text: "qubits initialised in |0⟩ experience single-qubit X errors"
        if p > 0:
            circuit.append("X_ERROR", all_data, p)
        
    return circuit

def get_readout_block(d: int, basis: str, A_mat_dict: Dict[int, np.ndarray], B_mat_dict: Dict[int, np.ndarray], q_map: Dict[Tuple[str, int], int], p: float = 0.0) -> stim.Circuit:
    circuit = stim.Circuit()
    n_half = int(len(q_map)/4)
    
    # =========================================================================
    # 1. ESTABLISH THE "MATRIX TO QUBIT" ORDER 
    # =========================================================================
    matrix_col_to_qubit_uid = []
    
    # Left Qubits
    for i in range(1, n_half + 1):
        matrix_col_to_qubit_uid.append(q_map[('L', i)])
    # Right Qubits
    for i in range(1, n_half + 1):
        matrix_col_to_qubit_uid.append(q_map[('R', i)])
        
    N_data = len(matrix_col_to_qubit_uid)

    # =========================================================================
    # 2. PERFORM MEASUREMENTS
    # =========================================================================
    if basis == 'X':
        # Measure in X basis
        # To simulate error 'p', we must flip the result using Z_ERROR before measuring
        if p > 0:
            circuit.append("Z_ERROR", matrix_col_to_qubit_uid, p)
        circuit.append("MX", matrix_col_to_qubit_uid)
        
    else: # basis == 'Z'
        # Measure in Z basis
        # To simulate error 'p', we must flip the result using X_ERROR before measuring
        if p > 0:
            circuit.append("X_ERROR", matrix_col_to_qubit_uid, p)
        circuit.append("M", matrix_col_to_qubit_uid)

    # =========================================================================
    # 3. BUILD LOOKUP
    # =========================================================================
    qubit_to_record_map = {}
    for index, qubit_uid in enumerate(matrix_col_to_qubit_uid):
        record_offset = -N_data + index
        qubit_to_record_map[qubit_uid] = stim.target_rec(record_offset)

    # =========================================================================
    # 4. BOUNDARY DETECTORS 
    # =========================================================================
    
    # Construct Matrices
    A_mat = sum(A_mat_dict.values()) % 2
    B_mat = sum(B_mat_dict.values()) % 2
    
    # Select Relevant Matrix
    if basis == 'X':
        relevant_H = np.hstack((A_mat, B_mat)) # X-check matrix
    else: 
        relevant_H = np.hstack((B_mat.T, A_mat.T)) # Z-check matrix

    for stab_idx in range(n_half):
        targets = []
        
        # A. Data Parity
        row = relevant_H[stab_idx]
        for col_idx, is_involved in enumerate(row):
            if is_involved:
                qubit_uid = matrix_col_to_qubit_uid[col_idx]
                targets.append(qubit_to_record_map[qubit_uid])
        
        # B. Last Ancilla & Coordinate Assignment
        if basis == 'X':
            # -----------------------------------------------------------------
            # Spatial: stab_idx (0 to n_half-1)
            # Time:    d
            # Type:    0
            # -----------------------------------------------------------------
            # Calculate offset for the X-ancilla from the previous cycle
            # (Assuming X-ancillas are the second half of the ancilla list in memory)
            ancilla_offset = -N_data - n_half + stab_idx
            targets.append(stim.target_rec(ancilla_offset))
            
            circuit.append("DETECTOR", targets, [stab_idx, d, 0])
            
        elif basis == 'Z':
            # -----------------------------------------------------------------
            # Spatial: stab_idx + n_half (n_half to 2*n_half-1)
            # Time:    d
            # Type:    1
            # -----------------------------------------------------------------
            # Calculate offset for the Z-ancilla from the previous cycle
            # (Assuming Z-ancillas are the first half of the ancilla list in memory)
            ancilla_offset = -N_data - 2*n_half + stab_idx
            targets.append(stim.target_rec(ancilla_offset))
            
            circuit.append("DETECTOR", targets, [stab_idx + n_half, d, 1])

    # =========================================================================
    # 5. LOGICAL OBSERVABLES (via qldpc)
    # =========================================================================
    css_code = qldpc.codes.CSSCode(np.hstack((A_mat, B_mat)), np.hstack((B_mat.T, A_mat.T)))
    
    if basis == 'X':
        log_ops_gf = css_code.get_logical_ops(qldpc.objects.Pauli.X)
    else:
        log_ops_gf = css_code.get_logical_ops(qldpc.objects.Pauli.Z)
        
    log_ops = np.array(log_ops_gf, dtype=int)
    num_logicals = log_ops.shape[0]

    for k in range(num_logicals):
        targets = []
        log_op_vec = log_ops[k]
        for col_idx, is_involved in enumerate(log_op_vec):
            if is_involved:
                qubit_uid = matrix_col_to_qubit_uid[col_idx]
                targets.append(qubit_to_record_map[qubit_uid])
        circuit.append("OBSERVABLE_INCLUDE", targets, k)
    
    return circuit

def get_memory_experiment_circuit(d: int, basis: str, A_mat_dict: Dict[int, np.ndarray], B_mat_dict: Dict[int, np.ndarray], q_map: Dict[Tuple[str, int], int], p: float = 0.0) -> stim.Circuit:
    
    # 1. Initialize Circuit
    circuit = stim.Circuit()
    n_half = int(len(q_map)/4)
    
    # 2. Add Data Initialization Block
    circuit += get_data_init_block(basis, q_map, p)

    weight = len(A_mat_dict) + len(B_mat_dict)

    if weight==6:
        # This ensures the first cycle's Z-ancillas have the same prep noise as later cycles.
        z_ancillas = [q_map[('Z', i)] for i in range(1, n_half + 1)]
        circuit.append("R", z_ancillas)
        if p > 0:
            # Match the noise model: Reset (R) followed by scrambling
            circuit.append("X_ERROR", z_ancillas, p)
        noiseless_cycle = get_one_cycle_circuit(A_mat_dict, B_mat_dict, q_map)
    elif weight==8:
        noiseless_cycle = get_one_cycle_circuit_weight_8(A_mat_dict, B_mat_dict, q_map)
    else:
        raise ValueError("Unsupported weight. Only 6 and 8 are supported in this function.")


    noisy_cycle = get_noisy_circuit(noiseless_cycle, p)

    noisy_first_cycle = stim.Circuit()
    first_tick_seen = False
    for inst in noisy_cycle:
        # Skip idle noise on data qubits during the first initialization step
        if not first_tick_seen and inst.name == "DEPOLARIZE1":
            continue
        noisy_first_cycle.append(inst)
        if inst.name == "TICK":
            first_tick_seen = True

    noisy_last_cycle = stim.Circuit()
    # If d=1, the last cycle IS the first cycle, so we must strip both ends!
    insts = list(noisy_first_cycle) if d == 1 else list(noisy_cycle)
    
    if len(insts) > 0 and insts[-1].name == "TICK":
        insts = insts[:-1] 
    if len(insts) > 0 and insts[-1].name == "DEPOLARIZE1":
        insts = insts[:-1] 
        
    for inst in insts:
        noisy_last_cycle.append(inst)

    # 4. Loop over 'd' cycles
    for cycle_idx in range(d):
        
        # Determine which cycle variant to use based on the boundary
        if d == 1:
            circuit += noisy_last_cycle     # Strips front AND back
        elif cycle_idx == 0:
            circuit += noisy_first_cycle    # Strips front only
        elif cycle_idx == d - 1:
            circuit += noisy_last_cycle     # Strips back only
        else:
            circuit += noisy_cycle          # Standard middle cycle
        
        # B. Add "Bulk" Detectors
        #    Note: Offsets must still be calculated relative to the *end* of the stack.
        
        # --- Handle X-Stabilizers (Type 0) ---
        for k in range(n_half):
            targets = []
            
            # X-ancillas are measured LAST in the cycle: Offset = -n_half + k
            current_offset = -n_half + k
            targets.append(stim.target_rec(current_offset))
            
            if cycle_idx > 0:
                prev_offset = current_offset - (2 * n_half)
                targets.append(stim.target_rec(prev_offset))
                circuit.append("DETECTOR", targets, [k, cycle_idx, 0])
            
            elif cycle_idx == 0:
                # Boundary condition: Check only if basis is deterministic for X
                if basis == 'X':
                    circuit.append("DETECTOR", targets, [k, cycle_idx, 0])

        # --- Handle Z-Stabilizers (Type 1) ---
        for k in range(n_half):
            targets = []
            
            # Z-ancillas are measured FIRST in the cycle: Offset = -2*n_half + k
            current_offset = -2 * n_half + k
            targets.append(stim.target_rec(current_offset))
            
            if cycle_idx > 0:
                prev_offset = current_offset - (2 * n_half)
                targets.append(stim.target_rec(prev_offset))
                circuit.append("DETECTOR", targets, [k + n_half, cycle_idx, 1])
                
            elif cycle_idx == 0:
                # Boundary condition: Check only if basis is deterministic for Z
                if basis == 'Z':
                    circuit.append("DETECTOR", targets, [k + n_half, cycle_idx, 1])

    # 5. Add Readout Block
    circuit += get_readout_block(d, basis, A_mat_dict, B_mat_dict, q_map, p)

    circuit.append("TICK")
    
    return circuit