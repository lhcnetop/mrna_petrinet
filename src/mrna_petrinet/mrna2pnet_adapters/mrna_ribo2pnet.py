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
        transitions_dict = network_dict["transitions"]
        simulation_parameters = json_input["simulation_parameters"]
        
        # Get ribosome parameters
        ribosome_params = simulation_parameters.get("ribosome_parameters", {})
        
        # Add free ribosomes place
        places_dict['p_free_ribosomes'] = ribosome_params.get('initial_ribosomes', 50)
        
        # Modify transitions to include ribosome consumption/production
        chains = json_input["chains"]
        for chain_obj in chains:
            chain_name = chain_obj["name"]
            chain = chain_obj["sequence"]
            
            # Find first transition (consumes ribosome)
            first_transition_name = f't_{chain_name}_t1'
            if first_transition_name in transitions_dict:
                if 'consume' not in transitions_dict[first_transition_name]:
                    transitions_dict[first_transition_name]['consume'] = {}
                transitions_dict[first_transition_name]['consume']['p_free_ribosomes'] = 1
            
            # Find last transition (produces ribosome back)
            last_transition_name = f't_{chain_name}_t{len(chain)}'
            if last_transition_name in transitions_dict:
                if 'produce' not in transitions_dict[last_transition_name]:
                    transitions_dict[last_transition_name]['produce'] = {}
                transitions_dict[last_transition_name]['produce']['p_free_ribosomes'] = 1
        
        return network_dict 