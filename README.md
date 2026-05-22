This repository contains the quantum error-correcting codes constructed in the paper [Insert Paper Title Here], as well as the simulator used to run circuit-level memory experiments for these codes.

## Repository Structure

### `code_dict/`
This folder contains `.json` files that include the metadata and building blocks for all the codes constructed in the paper. The files are named dynamically based on the specific tables they reference:

* **Table 1 and Table 5 Codes:** `code_n{n}_k{k}_d{d}_l{l}_m{m}_s{s}_a{a1}_{a2}..._b{b1}_{b2}....json`
* **Table 3 Base Codes:** `base_n{n}_k{k}_d{d}_l{l}_m{m}_a{a1}_{a2}..._b{b1}_{b2}....json`
* **Table 3 Cover Codes:** `cover_n{n}_k{k}_d{d}_l{l}_m{m}_s{s}_a{a1}_{a2}..._b{b1}_{b2}....json`

### Parameter Definitions:

* **n, k, d:** The standard quantum code parameters.
* **l, m, s:** The GAP identifiers.
* **a, b:** The group algebra elements corresponding to the respective table.

### `sinter_results/`
This folder contains `.csv` files that include the circuit-level simulation results from `sinter` used to generate the plots in Figures 2 and 3. The files are named as follows:

`results_{codename}_{basis}.csv`

### Parameter Definitions:

* **codename, basis:** The `codename` corresponds to the name of the file from the `code_dict/` folder, and the `basis` (`X` or `Z`) specifies the basis in which the memory experiment simulation is conducted in `stim`.

### `syndrome_extraction_circuits.py`
This script contains all the functions required to generate the memory experiment circuits for a given code.

### `run_circuit_simulations.py`
This is the main script used to run the circuit-level simulations.

## Installation
First clone the GitHub repo
```bash
git clone https://github.com/aaydinnnn/Coset2BGACodes.git
cd Coset2BGACodes
```
Then create the virtual envrionment. You can use your preferred environment manager (e.g., `venv`, `conda`, `mamba`, `micromamba`). 
```bash
mamba create -n qecsim -c conda-forge python=3.11 pip rust
```
Activate your environment:
```bash
mamba activate qecsim
```
Install the required dependencies:
```bash
pip install -r requirements.txt
```

## Execution

The primary execution script is `run_circuit_simulations.py`. This script automatically extracts the required `.json` files based on your input parameters, validates the code properties (verifying the rank of $H_x$ and $H_z$ over GF(2) to confirm logical qubit dimensions), constructs the space-time `stim` tasks, and runs the parallelized `sinter` simulation.

### 1. Configure the Simulation Inputs
Open `run_circuit_simulations.py` in your text editor. Under the `--- INPUT CONFIGURATION ---` section, you must define the target code parameters to match the table from the paper you wish to reproduce.

* **`source_table`**: Selects the naming convention for the `.json` file. 
  * *Options:* `'Table 1'`, `'Table 3_Base'`, `'Table 3_Cover'`, `'Table 5'`.
* **`n, k, d`**: Standard quantum code parameters.
* **`l, m, s`**: The GAP mathematical identifiers used to construct the code block.
* **`a, b`**: The group algebra elements (provided as lists, e.g., `a = [1, 34, 48]`).

*Example configuration for a Table 1 code:*
```python
source_table = 'Table 1'
n = 48
k = 8
d = 6
l = 384
m = 512
s = 53  

a = [1, 34, 48]
b = [1, 6, 12]
```
### 2. Configure Simulation Parameters
Scroll down slightly to configure the runtime parameters for the memory experiment:

* **`basis`**: Set to `'X'` to simulate Z-errors, or `'Z'` to simulate X-errors.
* **`error_rates`**: A list of physical error probabilities $p$ to sweep across (e.g., `[0.002, 0.004]`).
* **`num_cycle`**: The number of syndrome extraction cycles.
* **`num_workers`**: The number of CPU cores to allocate to `sinter`.
* **`max_shots` / `max_errors`**: The early-stopping criteria for the Monte Carlo sampler.

### 3. Choose the Decoder
The script supports multiple advanced decoders. Under the `CHOOSE DECODER:` section near the bottom of the file, modify the `active_decoders` list to select which decoders you want `sinter` to benchmark simultaneously.

```python
# Options: 'bposd', 'relay-bp', 'beam32_340iters'
active_decoders = ['bposd']
```
**Available Decoders:**

* **`bposd`**: Standard Belief Propagation with Ordered Statistics Decoding (defaults to Combination Sweep `osd_cs`, Order 10).
* **`relay-bp`**: A real time decoder from the [paper](https://arxiv.org/abs/2506.01779).
* **`beam32_340iters` (Optional):** Decoder from the [paper](https://arxiv.org/abs/2512.07057) . *(Note: This decoder will only be available if the `BeamSearchDecoder` module is present in the repository. The repository for this decoder is available in the [link](https://github.com/ionq-publications/BeamSearchDecoder)).*

### 4. Output and Checkpoints

* **Sanity Check:** Upon execution, the script will print a sanity check confirming the loaded matrix dimensions match the expected $k$ logical qubits.
* **Checkpoints:** As the simulation progresses, results are continuously backed up to a `checkpoint_{codename}_{basis}.csv` file to prevent data loss.
* **Final Results:** Once the target `max_errors` or `max_shots` are reached for all physical error rates, the final data is saved in the repository as `results_{codename}_{basis}.csv`.



