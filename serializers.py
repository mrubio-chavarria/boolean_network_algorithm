#!/home/mario/Projects/boolean_2/software/venv/bin/venv python

# Libraries
import ast

# Classes
class Serializer:
    """
    DESCRIPTION:
    A class to convert the networks into a JSON-like format.
    """
    def __init__(self, network):
        """
        DESCRIPTION:
        Class constructor.
        :param network: [dict] the boolean network to serialize.
        """
        self.network = network
        self.serialize_network()

    def serialize_network(self):
        """
        DESCRIPTION:
        A method to prepare the network for being visualised as a JSON file.
        Important, every set/fronzenset will be substituted by a list.
        """
        # Convert the sets to lists
        for node in self.network['pre_pathways'].keys():
            # Prepare pre_pathways
            [pathway.update({'domain': list(pathway['domain'])}) 
                for pathway in self.network['pre_pathways'][node]['activators']]
            [pathway.update({'domain': list(pathway['domain'])}) 
                for pathway in self.network['pre_pathways'][node]['inhibitors']]
            # Prepare pathways
            [pathway.update({'domain': list(pathway['domain'])}) 
                for pathway in self.network['pathways'][node]['activators']]
            [pathway.update({'domain': list(pathway['domain'])}) 
                for pathway in self.network['pathways'][node]['inhibitors']]
        # Convert the conflict strings to dicts again
        for conflict in self.network['conflicts']:
            conflict['activator'] = conflict['activator'].replace('frozenset', '').replace('({', '[').replace('})', ']')
            conflict['activator'] = ast.literal_eval(conflict['activator'])
            conflict['inhibitor'] = conflict['inhibitor'].replace('frozenset', '').replace('({', '[').replace('})', ']')
            conflict['inhibitor'] = ast.literal_eval(conflict['inhibitor'])
            conflict['solution'] = conflict['solution'].replace('frozenset', '').replace('({', '[').replace('})', ']')
            conflict['solution'] = ast.literal_eval(conflict['solution'])


    def get_serialized_network(self):
        return self.network
