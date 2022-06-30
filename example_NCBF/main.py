#!/home/mario/Projects/boolean_2/software/venv/bin/venv python



# Libraries
import os
import sys
import json
import shutil
import pickle
import json
from random import choice
from pytictoc import TicToc
from datetime import datetime
from pytictoc import TicToc
from graph import Graph
from bn_utils import minterms2bnet
# from bn_utils import filter_boolean_networks
# from serializers import Serializer


# Functions
def main():
    """
    DESCRIPTION:
    Reduced version of the main in the outer folder to generate just the NCBF
    combination.
    """

    # Read inference information in the data.json file
    data = json.load(open('data.json'))
    # Delete data that are no longer required
    del data['load_stored_networks']
    del data['n_attractors']
    del data['delete_repeated']
    del data['min_conflicts']

    # Infer a new graph
    # Create the object graph
    print('Create the graph and generate all the networks')
    graph = Graph(**data)

    # The name of the methods describe the steps
    # Compute all the possible combinations of pathways from the graph
    # There are multiple groups of pathways compatible with the same
    # graph
    graph.obtain_pathways_from_graph()
    # Generate a priority matrix couple (act and inh) for each simulation to perform
    graph.generate_priority_matrices()
    # Generate all the possible networks made of NCBFs for each group of pathways
    graph.generate_NCBFs()

    # Print result example
    # Example
    print("Example")
    print("Pairs: ", graph.pathway_groups[0])
    print("--------")
    print("Network: ", graph.pre_networks[0])

    # Print networks to files
    networks_folder = "example_NCBF/printed_networks"
    networks = [
        "\n".join([f"{node}, " + minterms2bnet(graph.nodes, network[1][node]) for node in graph.nodes])
        for network in graph.pre_networks
    ]
    i = 0
    os.system(f"mkdir -p {networks_folder}")  # <------- ESTE COMANDO CREA LA CARPETA DENTRO DE example_NCBF. EN WINDOWS ES MEJOR COMENTAR Y CREAR LA CARPETA A MANO
    for network in networks:
        with open(networks_folder + f"/network_{i}.txt", "w") as file:
            file.write(network)
        i += 1



if __name__ == "__main__":
    # Execute main function
    main()

