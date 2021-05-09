#!/home/mario/Projects/boolean_2/software/venv/bin/venv python

# Libraries
import json
import os, shutil
import pickle
from datetime import datetime
from pytictoc import TicToc
from graph import Graph
from bn_utils import filter_boolean_networks


# Functions
def main():
    """
    DESCRIPTION:
    Main method, the guideline to execute the inference.
    """
    # Assess if it is needed to load a previous session
    load_session = False
    # Assess if the networks are going to be checked for 
    # repeated
    delete_repeated = False
    if not load_session:
        # Infer a new graph
        # Read problem information
        data = json.load(open('data.json'))
        # Create the object graph
        print('Create the graph and generate all the networks')
        graph = Graph(**data)
        graph.obtain_pathways_from_graph()
        graph.generate_priority_matrices()
        graph.generate_NCBFs()
        graph.generate_boolean_networks()
        # Perform inference (it takes time)
        print('Solve the conflicts of all the networks')
        graph.solve_conflicts()
        # Compute relevant features of the converging networks
        # and filter by number of desired attractors (total)
        graph.format_network()
        # Store the resulting networks
        filename = 'networks/networks_' + '_'.join(str(datetime.now()).split(' ')) + '.pickle'
        with open(filename, 'wb') as file:
            pickle.dump(graph.boolean_networks, file)
        # Prepare for the filtering
        boolean_networks = graph.boolean_networks
    else:
        # Load networks from a previous session (graph)
        print('Loading selected session')
        sessions_folder = 'networks/'
        # Load the networks
        boolean_networks = [pickle.load(open(sessions_folder + file, 'rb'))
            for file in os.listdir(sessions_folder)]
        boolean_networks = [it for sb in boolean_networks for it in sb]
        # Delete repeated
        if delete_repeated:
            # Take only those network that are not repeated
            expressions = []
            non_repeated_networks = []
            for boolean_network in boolean_networks:
                expression = ''.join(sorted(boolean_network['network'].values(), key=lambda x: x[0]))
                if expression not in expressions:
                    expressions.append(expression)
                    non_repeated_networks.append(boolean_network)
            # Delete the pickle files of all the networks
            shutil.rmtree(sessions_folder)
            os.mkdir(sessions_folder)
            # Store the filtered networks in one pickle file
            filename = 'networks/networks_' + '_'.join(str(datetime.now()).split(' ')) + '.pickle'
            with open(filename, 'wb') as file:
                pickle.dump(non_repeated_networks, file)
    # Analysis over the graph
    print('Filter the resulting networks')
    # Parameters to filter 
    attractors = ['0000', '1111']
    n_attractors = 16
    partial = False
    # Filtering
    boolean_networks = filter_boolean_networks(boolean_networks,
        attractors=attractors,
        n_attractors=n_attractors,
        partial=partial)
    print()


if __name__ == "__main__":
    # Execute main function
    main()

