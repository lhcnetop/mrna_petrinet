from mrna_petrinet.mrna2pnet_adapters.mrna2pnet import mRNA2PNetAdapter
from petri_net_core.petrinet.pnet import PNet
from petri_net_core.petrinet.pnet_simulation_helper import simulate_with_history_to_parquet

# TODO: Replace this with a real minimal valid mRNA input for your adapter
mrna_input = {
    "chains":[
        {
            "name":"chainA",
            "sequence":"MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
            "polipeptide_name":"preinsulin"
        }
    ],
    "simulation_parameters":{
        "initial_chains_marking":200,
        "max_protein_output_goal":100,
        "excess_aminoacids_factor":1
    }
}

adapter = mRNA2PNetAdapter(mrna_input)


print("Generated Petri net dict:")
pnet_dict = adapter.get_pnet_json_dict()
print(pnet_dict)

pnet = PNet(pnet_dict)

# Save simulation history to Parquet file in the data folder
data_filepath = "data/simulation_history.parquet"
simulate_with_history_to_parquet(pnet, num_steps=5, parquet_file=data_filepath)
print(f"Simulation history saved to {data_filepath}") 