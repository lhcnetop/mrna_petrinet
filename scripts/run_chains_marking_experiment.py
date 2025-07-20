from mrna_petrinet.mrna2pnet_adapters.mrna_ribo2pnet import mRNARibo2PNetAdapter
from petri_net_core.petrinet.pnet import PNet
from petri_net_core.petrinet.pnet_simulation_helper import simulate_with_history_to_parquet
import os
import copy
import time
from datetime import datetime
from multiprocessing import Pool, cpu_count
import logging

# Base input for mRNA simulation with ribosomes
base_mrna_ribo_input = {
    "chains": [
        {
            "name": "chainA",
            "sequence": "MALWMRLLPLLALLALWGRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
            "polipeptide_name": "preinsulin"
        }
    ],
    "simulation_parameters": {
        "initial_chains_marking": 100,  # This will be varied
        "max_protein_output_goal": 100,  # Fixed at 100
        "ribosome_parameters": {
            "initial_ribosomes": 2  # This will be varied
        }
    }
}

# Initial chains marking values to test (varying the substrate availability)
# Reversed order: higher values first
#initial_chains_marking_values = [500, 200, 100, 50, 20, 10]
initial_chains_marking_values = [400, 350, 300, 250, 200, 170, 150, 130, 100, 70, 50, 20, 10]

# Initial ribosome values to test
#initial_ribosome_values = [1, 5, 10, 50, 100, 200, 1000, 5000]
initial_ribosome_values = [2, 4, 6, 8, 20, 30, 40, 60, 70, 80, 90, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]

# Fixed max protein output goal
max_protein_output_goal = 100

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Create experiment subfolder
experiment_folder = "data/chains_marking_experiment"
os.makedirs(experiment_folder, exist_ok=True)
print(f"Created experiment folder: {experiment_folder}")

# Start total experiment timer
total_experiment_start_time = time.time()
print(f"\n=== Starting total experiment at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ===")

def run_single_simulation(initial_chains_marking, initial_ribosomes):
    """
    Run a single simulation with given parameters.
    
    Args:
        initial_chains_marking (int): Initial chains marking value
        initial_ribosomes (int): Initial ribosome count
    
    Returns:
        tuple: (success, data_filepath, simulation_time, error_message)
    """
    try:
        # Calculate excess mRNA chain factor as ratio (inverse of aminoacid factor)
        excess_mrna_chain_factor = initial_chains_marking / max_protein_output_goal
        
        simulation_start_time = time.time()
        print(f"Running simulation: {initial_chains_marking} chains, {initial_ribosomes} ribosomes (excess factor: {excess_mrna_chain_factor:.2f})")
        
        # Update base input with current parameters
        current_input = copy.deepcopy(base_mrna_ribo_input)
        current_input["simulation_parameters"]["initial_chains_marking"] = initial_chains_marking
        current_input["simulation_parameters"]["ribosome_parameters"]["initial_ribosomes"] = initial_ribosomes
        
        # Calculate max steps based on current chains marking
        chain_length = len(current_input["chains"][0]["sequence"])
        num_chains = current_input["simulation_parameters"]["initial_chains_marking"]
        max_steps = chain_length * (num_chains + 1)
        
        # Create the ribosome adapter
        ribo_adapter = mRNARibo2PNetAdapter(current_input)
        
        # Create PNet object
        pnet = PNet(ribo_adapter.get_pnet_json_dict())
        
        # Save simulation history to Parquet file with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        data_filepath = f"{experiment_folder}/chains_marking_simulation_history_{initial_ribosomes}_ribosomes_{excess_mrna_chain_factor:.2f}_excess_mrna_{timestamp}.parquet"
        simulate_with_history_to_parquet(pnet, num_steps=max_steps, parquet_file=data_filepath)
        
        simulation_time = time.time() - simulation_start_time
        print(f"Simulation completed: {data_filepath} ({simulation_time:.2f}s)")
        
        return True, data_filepath, simulation_time, None
        
    except Exception as e:
        error_msg = f"Error in simulation ({initial_chains_marking} chains, {initial_ribosomes} ribosomes): {str(e)}"
        print(error_msg)
        return False, None, 0, error_msg

if __name__ == '__main__':
    # Create parameter combinations for parallel processing
    parameter_combinations = []
    
    # Number of times to repeat each experiment
    num_repetitions = 10
    
    for initial_chains_marking in initial_chains_marking_values:
        for initial_ribosomes in initial_ribosome_values:
            # Add the same parameter combination multiple times
            for _ in range(num_repetitions):
                parameter_combinations.append((initial_chains_marking, initial_ribosomes))

    print(f"Total simulations to run: {len(parameter_combinations)}")
    print(f"Each parameter combination repeated {num_repetitions} times")

    # Determine number of processes (use 60% of available cores)
    num_processes = max(1, int(cpu_count() * 0.7))
    print(f"Using {num_processes} processes for parallel execution")

    # Run simulations in parallel
    print(f"\nStarting parallel execution of {len(parameter_combinations)} simulations...")
    
    # Track progress
    total_simulations = len(parameter_combinations)
    results = []
    
    # Use a simple counter for progress
    from multiprocessing import Value
    counter = Value('i', 0)
    
    def progress_callback(result):
        with counter.get_lock():
            counter.value += 1
            remaining = total_simulations - counter.value
            print(f"[{datetime.now().strftime('%H:%M:%S')}] Progress: {counter.value}/{total_simulations} completed ({remaining} remaining)")
        results.append(result)
    
    with Pool(processes=num_processes) as pool:
        # Submit all tasks and track progress
        async_results = []
        for params in parameter_combinations:
            async_result = pool.starmap_async(run_single_simulation, [params], callback=progress_callback)
            async_results.append(async_result)
        
        # Wait for all to complete
        for async_result in async_results:
            async_result.wait()

    # Process results
    successful_simulations = 0
    failed_simulations = 0
    total_simulation_time = 0

    for success, data_filepath, simulation_time, error_msg in results:
        if success:
            successful_simulations += 1
            total_simulation_time += simulation_time
        else:
            failed_simulations += 1
            print(f"Failed simulation: {error_msg}")

    print(f"\n=== Parallel execution completed ===")
    print(f"Successful simulations: {successful_simulations}")
    print(f"Failed simulations: {failed_simulations}")
    print(f"Total simulation time: {total_simulation_time:.2f} seconds")

    # End total experiment timer
    total_experiment_time = time.time() - total_experiment_start_time
    print(f"\n=== Total experiment completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ===")
    print(f"=== Total experiment time: {total_experiment_time:.2f} seconds ({total_experiment_time/60:.2f} minutes) ===")

    print("\nAll chains marking experiments completed!") 

