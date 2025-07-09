from .mrna2pnet import mRNA2PNetAdapter

class mRNARibo2PNetAdapter(mRNA2PNetAdapter):
    def __init__(self, json_input: dict):
        # Call parent constructor to get the base Petri net
        super().__init__(json_input)
        # Add ribosome places and modify transitions
        self.places_and_transitions = self.add_ribosome_functionality(self.places_and_transitions, json_input)
    
    def add_ribosome_functionality(self, network_dict: dict, json_input: dict) -> dict:
        """
        Add ribosome functionality to an existing Petri net network.
        
        Args:
            network_dict: Dictionary containing places and transitions from parent class
            json_input: Complete input dictionary with chains and simulation parameters
            
        Returns:
            Updated network dictionary with ribosome places and modified transitions
        """
        places_dict = network_dict["places"]
        transitions_array = network_dict["transitions"]
        simulation_parameters = json_input["simulation_parameters"]
        
        # Add free ribosomes place
        places_dict['p_free_ribosomes'] = simulation_parameters.get('initial_free_ribosomes', 50)
        
        # Modify transitions to include ribosome consumption/production
        chains = json_input["chains"]
        for chain_obj in chains:
            chain_name = chain_obj["name"]
            chain = chain_obj["sequence"]
            
            # Find first transition (consumes ribosome)
            first_transition_name = f't_{chain_name}_t1'
            for transition in transitions_array:
                if transition['name'] == first_transition_name:
                    if 'consume' not in transition:
                        transition['consume'] = {}
                    transition['consume']['p_free_ribosomes'] = 1
                    break
            
            # Find last transition (produces ribosome back)
            last_transition_name = f't_{chain_name}_t{len(chain)}'
            for transition in transitions_array:
                if transition['name'] == last_transition_name:
                    if 'produce' not in transition:
                        transition['produce'] = {}
                    transition['produce']['p_free_ribosomes'] = 1
                    break
        
        return network_dict 