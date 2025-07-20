from mrna_petrinet.mrna2pnet_adapters.mrna_ribo2pnet import mRNARibo2PNetAdapter
from petri_net_core.petrinet.pnet import PNet
from petri_net_core.petrinet.pnet_simulation_helper import simulate_with_history_to_parquet
import os
import copy
import time

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
        "initial_chains_marking": 100,
        "max_protein_output_goal": 100,
        "ribosome_parameters": {
            "initial_ribosomes": 2  # This will be varied
        }
    }
}

chain_length = len(base_mrna_ribo_input["chains"][0]["sequence"])
num_chains = base_mrna_ribo_input["simulation_parameters"]["initial_chains_marking"]
max_steps = chain_length * (num_chains + 1)

# Initial ribosome values to test
initial_ribosome_values = [1, 2, 3, 10, 30, 100, 200, 1000]

# Max protein output goal values to test
max_protein_output_goals = [10,20,30,40,50,60,70,80,90,100,110,120]

print(max_steps)

# Create experiment subfolder
experiment_folder = "data/ribosome_experiment_2"
os.makedirs(experiment_folder, exist_ok=True)
print(f"Created experiment folder: {experiment_folder}")

# Calculate excess aminoacid factor as ratio
initial_chains_marking = base_mrna_ribo_input["simulation_parameters"]["initial_chains_marking"]

print("Running multiple ribosome simulations with different initial ribosome values...")
print(f"Testing ribosome values: {initial_ribosome_values}")
print(f"Testing max protein output goals: {max_protein_output_goals}")

for initial_ribosomes in initial_ribosome_values:
    ribosome_loop_start_time = time.time()
    print(f"\n=== Starting simulations for {initial_ribosomes} initial ribosomes ===")
    
    for max_protein_output_goal in max_protein_output_goals:
        # Calculate excess aminoacid factor as ratio
        excess_aminoacid_factor = max_protein_output_goal / initial_chains_marking
        
        simulation_start_time = time.time()
        print(f"\n--- Running simulation with {initial_ribosomes} initial ribosomes and {max_protein_output_goal} max protein output goal (excess factor: {excess_aminoacid_factor:.2f}) ---")
        
        # Update the input with current ribosome value and max protein output goal
        current_input = copy.deepcopy(base_mrna_ribo_input)
        current_input["simulation_parameters"]["ribosome_parameters"]["initial_ribosomes"] = initial_ribosomes
        current_input["simulation_parameters"]["max_protein_output_goal"] = max_protein_output_goal
        
        # Create the ribosome adapter
        ribo_adapter = mRNARibo2PNetAdapter(current_input)
        
        print(f"Generated Petri net dict with {initial_ribosomes} initial ribosomes and {max_protein_output_goal} max protein output goal:")
        pnet_dict = ribo_adapter.get_pnet_json_dict()
        
        # Create PNet object
        pnet = PNet(pnet_dict)
        
        # Save simulation history to Parquet file in the experiment folder
        data_filepath = f"{experiment_folder}/ribosome_simulation_history_{initial_ribosomes}_ribosomes_{excess_aminoacid_factor:.2f}_excess.parquet"
        simulate_with_history_to_parquet(pnet, num_steps=max_steps, parquet_file=data_filepath)
        
        simulation_time = time.time() - simulation_start_time
        print(f"Ribosome simulation history saved to {data_filepath}")
        print(f"Simulation completed in {simulation_time:.2f} seconds")
    
    ribosome_loop_time = time.time() - ribosome_loop_start_time
    print(f"\n=== Completed all simulations for {initial_ribosomes} initial ribosomes in {ribosome_loop_time:.2f} seconds ===")

print("\nAll simulations completed!") 