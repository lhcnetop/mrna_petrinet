from mrna_petrinet.mrna2pnet_adapters.mrna_ribo2pnet import mRNARibo2PNetAdapter
from petri_net_core.petrinet.pnet import PNet
from petri_net_core.petrinet.pnet_simulation_helper import simulate_with_history_to_parquet

# Input for mRNA simulation with ribosomes
mrna_ribo_input = {
    "chains": [
        {
            "name": "chainA",
            "sequence": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
            "polipeptide_name": "preinsulin"
        }
    ],
    "simulation_parameters": {
        "initial_chains_marking": 2,
        "max_protein_output_goal": 2,
        "excess_aminoacids_factor": 1,
        "ribosome_parameters": {
            "initial_ribosomes": 2
        }
    }
}

# Create the ribosome adapter
ribo_adapter = mRNARibo2PNetAdapter(mrna_ribo_input)

print("Generated Petri net dict with ribosomes:")
pnet_dict = ribo_adapter.get_pnet_json_dict()
#print(pnet_dict)

# Create PNet object
pnet = PNet(pnet_dict)

# Save simulation history to Parquet file in the data folder
data_filepath = "data/ribosome_simulation_history.parquet"
simulate_with_history_to_parquet(pnet, num_steps=1000, parquet_file=data_filepath)
print(f"Ribosome simulation history saved to {data_filepath}") 