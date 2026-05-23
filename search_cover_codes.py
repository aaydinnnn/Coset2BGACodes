import multiprocessing as mp
from sage.all import GF, matrix, Matrix, block_matrix
from sage.interfaces.gap import gap, Gap
import itertools
import ast
import os

# ==========================================
# 1. PARAMETERS & CONFIGURATION
# ==========================================

F2 = GF(2)
NUM_TRIALS = 1000  #Number of trials for QDistRnd algorithm
num_cores = 4 #Number of CPU cores to be used in parallized search

lift_size = 2 # The code will search for lift_size-fold cover codes of the given base code


# Identify the base code (do not include 1 here for the group algebra elements a and b)
base_l = 51
base_m = 1
a_base = [14, 45, 35]
b_base = [27, 32, 50]

# Automatically prepend 1 to the arrays
a_base = [1] + a_base
b_base = [1] + b_base


# ==========================================
# 2. COVER GROUP DISCOVERY VIA GAP
# ==========================================
def find_target_groups(base_l, base_m, lift_size):
    """Dynamically finds all valid cover groups and all valid normal subgroups."""
    target_size = base_l * lift_size 
    print(f"Finding covers of size {target_size} for base group SmallGroup({base_l}, {base_m})...\n")

    gap_script = f"""
    valid_covers := [];
    
    FindValidCovers := function()
        local G_base, base_id, ts, num_groups, i, G_cover, subs, normal, s, H, Q;
        
        G_base := SmallGroup({base_l}, {base_m});
        base_id := IdGroup(G_base);
        ts := {target_size};
        num_groups := NrSmallGroups(ts);

        for i in [1..num_groups] do
            G_cover := SmallGroup(ts, i);
            
            subs := AllSubgroups(G_cover);
            normal := Filtered(subs, obj -> IsNormal(G_cover, obj));
            
            for s in [1..Length(normal)] do
                H := normal[s];
                
                if Size(H) = {lift_size} then
                    Q := G_cover / H;
                    if IdGroup(Q) = base_id then
                        # Collect all valid configurations (no break here!)
                        Add(valid_covers, [i, StructureDescription(G_cover), StructureDescription(H), s]);
                    fi;
                fi;
            od;
        od;
    end;
    
    FindValidCovers();
    """
    gap.eval(gap_script)

    num_valid = int(gap.eval("Length(valid_covers);"))
    targets = []

    if num_valid == 0:
        print("No valid cover groups found for these parameters.")
    else:
        print(f"Found {num_valid} valid cover configuration(s):")
        print("-" * 65)
        
        for i in range(1, num_valid + 1):
            g_id = int(gap.eval(f"valid_covers[{i}][1];"))
            g_desc = gap.eval(f"valid_covers[{i}][2];").strip('"')
            h_desc = gap.eval(f"valid_covers[{i}][3];").strip('"')
            s_idx = int(gap.eval(f"valid_covers[{i}][4];"))
            
            print(f"Cover Group (l={target_size}, m={g_id:<2}) | Structure: {g_desc:10}")
            print(f"Subgroup    (s={s_idx:<2})             | Structure: {h_desc:10}")
            print("-" * 65)
            
            targets.append((target_size, g_id, s_idx))
            
    return targets

def get_group_info(base_l, base_m, l, m, s):
    """Initializes the main GAP process for a specific target group and establishes the isomorphism."""
    gap.eval(f"G_base := SmallGroup({base_l}, {base_m});")
    gap.eval("elems_base := Elements(G_base);")
    
    gap.eval(f"G := SmallGroup({l}, {m});")
    gap.eval(f"H := Filtered(AllSubgroups(G), h -> IsNormal(G,h))[{s}];")
    gap.eval("elemsG := Elements(G);") 
    
    gap.eval("nat := NaturalHomomorphismByNormalSubgroup(G, H);")
    gap.eval("Q := Image(nat);")
    gap.eval("iso := IsomorphismGroups(Q, G_base);")
    
    struct_g = gap.eval("StructureDescription(G)").strip('"')
    struct_h = gap.eval("StructureDescription(H)").strip('"')
    size_h = int(gap.eval("Size(H)"))
    num_elements = int(gap.eval("Length(elemsG);"))
    return struct_g, struct_h, size_h, num_elements

def evaluate_base_code(base_l, base_m, a_base, b_base):
    """Calculates the n, k, d parameters of the base code before lifting."""
    gap.eval(f"G := SmallGroup({base_l}, {base_m});")
    gap.eval("elemsG := Elements(G);")
    base_num_elems = int(gap.eval("Length(elemsG);"))
    
    A_base = sum(calculate_permutation_matrix(gap, k, True, base_num_elems) for k in a_base)
    B_base = sum(calculate_permutation_matrix(gap, k, False, base_num_elems) for k in b_base)
    
    Hx = block_matrix(F2, [[A_base, B_base]])
    Hz = block_matrix(F2, [[B_base.transpose(), A_base.transpose()]])
    
    n_bits = Hx.ncols()
    k_base = n_bits - Hx.rank() - Hz.rank()
    
    d_base = 0
    if k_base > 0:
        gap.eval('LoadPackage("QDistRnd");')
        Hx_list = [[int(val) for val in row] for row in Hx]
        Hz_list = [[int(val) for val in row] for row in Hz]
        gap.eval(f"Hx_gap := {Hx_list} * Z(2)^0;")
        gap.eval(f"Hz_gap := {Hz_list} * Z(2)^0;")
        
        dz = int(gap.eval(f"DistRandCSS(Hx_gap, Hz_gap, {NUM_TRIALS}, 0);"))
        dx = int(gap.eval(f"DistRandCSS(Hz_gap, Hx_gap, {NUM_TRIALS}, 0);"))
        d_base = min(dx, dz)
        
    return n_bits, k_base, d_base

# ==========================================
# 3. MATRIX CONSTRUCTION
# ==========================================
def calculate_permutation_matrix(local_gap, element_index, is_left_action, num_elements):
    local_gap.eval(f"g := elemsG[{element_index}];")
    sigma = []
    for i in range(1, num_elements + 1):
        local_gap.eval(f"rep := elemsG[{i}];")
        if is_left_action:
            local_gap.eval("prod := g * rep;")
        else:
            local_gap.eval("prod := rep * g;")
        
        gap_cmd = "idx_val := 1;; while prod <> elemsG[idx_val] do idx_val := idx_val + 1;; od;; idx_val;"
        res = local_gap.eval(gap_cmd)
        
        try:
            target_index = int(res) - 1
        except ValueError:
            target_index = int(local_gap.eval("idx_val;")) - 1
            
        sigma.append(target_index)

    P = matrix(F2, num_elements, num_elements)
    for i, j in enumerate(sigma):
        P[j, i] = 1
    return P

def get_elements_in_coset(local_gap, base_idx):
    """Finds the exact coset in G that maps to the base element at base_idx via the isomorphism."""
    local_gap.eval(f"g_base := elems_base[{base_idx}];")
    local_gap.eval("q_elem := PreImagesRepresentative(iso, g_base);")
    local_gap.eval("coset_elems := Elements(PreImages(nat, q_elem));")
    
    indices_str = local_gap.eval("List(coset_elems, e -> Position(elemsG, e));")
    clean_str = indices_str.replace('\n', '').strip()
    
    return ast.literal_eval(clean_str)

# ==========================================
# 4. PARALLEL WORKER LOGIC
# ==========================================
def worker_init(w_l, w_m, w_s):
    """Initializes a private GAP instance for each CPU core with current l, m, s."""
    global worker_gap
    worker_gap = Gap() 
    worker_gap.eval(f"G := SmallGroup({w_l}, {w_m});")
    worker_gap.eval(f"H := Filtered(AllSubgroups(G), h -> IsNormal(G,h))[{w_s}];")
    worker_gap.eval("elemsG := Elements(G);")
    worker_gap.eval("Cs := LeftCosets(G, H);")
    worker_gap.eval('LoadPackage("QDistRnd");')

def process_batch(a_batch, all_b_covers, L_lookup, R_lookup):
    results = []
    for a_cover in a_batch:
        A = sum(L_lookup[k] for k in a_cover)
        for b_cover in all_b_covers:
            B = sum(R_lookup[k] for k in b_cover)

            Hx = block_matrix(F2, [[A, B]])
            Hz = block_matrix(F2, [[B.transpose(), A.transpose()]])

            n_bits = Hx.ncols()
            k = n_bits - Hx.rank() - Hz.rank()

            if k > 0:
                Hx_list = [[int(val) for val in row] for row in Hx]
                Hz_list = [[int(val) for val in row] for row in Hz]

                worker_gap.eval(f"Hx_gap := {Hx_list} * Z(2)^0;")
                worker_gap.eval(f"Hz_gap := {Hz_list} * Z(2)^0;")
                
                dz = int(worker_gap.eval(f"DistRandCSS(Hx_gap, Hz_gap, {NUM_TRIALS}, 0);"))
                dx = int(worker_gap.eval(f"DistRandCSS(Hz_gap, Hx_gap, {NUM_TRIALS}, 0);"))
                d = min(dx, dz)
                
                print(f"[Core {os.getpid()}] Found: k={k}, d={d}")
                results.append((k, d, a_cover, b_cover))
    return results

# ==========================================
# 5. MAIN EXECUTION
# ==========================================
if __name__ == "__main__":
    print("=" * 70)
    print(f"PROCESSING JOB: base_l={base_l}, base_m={base_m}, lift_size={lift_size}")
    print("=" * 70)
        
    OUTPUT_FILE = f"cover_codes_results_base_l_{base_l}_base_m_{base_m}_a_base_{a_base}_b_base_{b_base}_lift_size_{lift_size}.txt"

    # 1. Discover Target Groups Automatically
    TARGET_GROUPS_RAW = find_target_groups(base_l, base_m, lift_size)

    # Deduplicate while preserving order (just in case)
    TARGET_GROUPS = list(dict.fromkeys(TARGET_GROUPS_RAW))

    # 2. Evaluate Base Code
    print("Evaluating base code parameters...")
    base_n, base_k, base_d = evaluate_base_code(base_l, base_m, a_base, b_base)
    base_struct = gap.eval(f"StructureDescription(SmallGroup({base_l}, {base_m}))").strip('"')
        
    # 3. Open Output File & Write Base Info
    with open(OUTPUT_FILE, "w") as f:
        f.write("=" * 50 + "\n")
        f.write("QLDPC COVER CODE SEARCH RESULTS\n")
        f.write("=" * 50 + "\n")
        f.write(f"BASE GROUP: SmallGroup({base_l}, {base_m}) ({base_struct})\n")
        f.write(f"LIFT SIZE:  {lift_size}\n")
        f.write("BASE CODE PARAMS:\n")
        f.write(f"  a_base = {a_base}\n")
        f.write(f"  b_base = {b_base}\n")
        f.write(f"  [n={base_n}, k={base_k}, d={base_d}]\n")
        f.write("=" * 50 + "\n\n")

    print(f"\nBase info [n={base_n}, k={base_k}, d={base_d}] written to {OUTPUT_FILE}.")
    print(f"Starting parallel search over {len(TARGET_GROUPS)} unique target configurations...\n")

    # 4. Loop Through All Target Groups for this Run
    for (l, m, s) in TARGET_GROUPS:
        g_struct, h_struct, size_h, num_elements = get_group_info(base_l, base_m, l, m, s)
            
        # Prepare header block
        header_str = (
            "-" * 50 + "\n" +
            f"TARGET: Group SmallGroup({l}, {m}) ({g_struct})\n" +
            f"SUBGROUP: Index {s} ({h_struct})\n" +
            f"LIFTING: Factor {size_h} | Total Qubits: {2 * num_elements}\n" +
            "-" * 50 + "\n"
        )
            
        print(header_str, end="")
        with open(OUTPUT_FILE, "a") as f:
            f.write(header_str)

        print("Generating lookup tables...")
        L_lookup = {i: calculate_permutation_matrix(gap, i, True, num_elements) for i in range(1, num_elements + 1)}
        R_lookup = {i: calculate_permutation_matrix(gap, i, False, num_elements) for i in range(1, num_elements + 1)}
        print("Tables ready.")

        a_opts = [[1] if c == 1 else get_elements_in_coset(gap, c) for c in a_base]
        b_opts = [[1] if c == 1 else get_elements_in_coset(gap, c) for c in b_base]
        all_a = list(itertools.product(*a_opts))
        all_b = list(itertools.product(*b_opts))

        print(f"Distributing {len(all_a)} 'a' batches across {num_cores} cores.")

        chunk_size = max(1, len(all_a) // num_cores)
        batches = [all_a[i:i + chunk_size] for i in range(0, len(all_a), chunk_size)]

        # Multiprocessing Pool
        with mp.Pool(processes=num_cores, initializer=worker_init, initargs=(l, m, s)) as pool:
            task_args = [(batch, all_b, L_lookup, R_lookup) for batch in batches]
            final_results_nested = pool.starmap(process_batch, task_args)

        # 5. Final Output & File Append for this specific (l, m, s) target
        final_results = [item for sublist in final_results_nested for item in sublist]
        final_results.sort(key=lambda x: x[1], reverse=True)
            
        result_summary = f"Found {len(final_results)} codes with k > 0.\n"
        print(result_summary, end="")
            
        with open(OUTPUT_FILE, "a") as f:
            f.write(result_summary)
            for k, d, a_c, b_c in final_results:
                line = f"[n={2*num_elements}, k={k}, d={d}] | a_cover={a_c} | b_cover={b_c}\n"
                print(line, end="")
                f.write(line)
            f.write("\n\n")
        