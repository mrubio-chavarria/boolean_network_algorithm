#!/home/mario/Projects/boolean_2/software/venv/bin/venv python

# Libraries
import json
import os, shutil
import pickle
import json
from random import choice
from pytictoc import TicToc
from datetime import datetime
from pytictoc import TicToc
from graph import Graph
from bn_utils import filter_boolean_networks
from serializers import Serializer


# Functions
def main():
    """
    DESCRIPTION:
    Main method, the guideline to execute the inference.
    """
    # Read problem information
    data = json.load(open('data.json'))  # Load all the information
    load_session = data['load_session']  # Assess if it is needed to load a previous session
    del data['load_session']
    partial = data['partial']  # Parameter to relax the conditions during network filtering
    del data['partial']
    attractors = data['attractors']  # Parameters attractors searched in the networks
    del data['attractors']
    n_attractors = data['n_attractors']  # Total number of attractors searched in the networks
    del data['n_attractors']
    delete_repeated = data['delete_repeated']  # assess whether to elete repeated networks or not
    del data['delete_repeated']
    if not load_session:
        # Infer a new graph
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
        print(f'Storing networks in: {filename}')
        t = TicToc()
        t.tic()
        with open(filename, 'wb') as file:
            # Transform the network in a format that can be written in 
            # a JSON file
            serialized_networks = [Serializer(net).get_serialized_network() 
                for net in graph.boolean_networks]
            pickle.dump(serialized_networks, file)
        # Prepare for the filtering
        boolean_networks = graph.boolean_networks
        t.toc('Computation time: ')
        print('Networks stored')
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
            # Keep the non repeated and store them
            boolean_networks = non_repeated_networks
            filename = 'networks/networks_' + '_'.join(str(datetime.now()).split(' ')) + '.pickle'
            with open(filename, 'wb') as file:
                pickle.dump(non_repeated_networks, file)
    # Analysis over the graph
    # # Filtering
    print('Filter the resulting networks')
    boolean_networks = filter_boolean_networks(boolean_networks,
                                                attractors=attractors,
                                                n_attractors=n_attractors,
                                                partial=partial)
    # Print a random network
    example_filename = 'example.json'
    with open(example_filename, 'w') as file:
        json.dump(choice(boolean_networks), file)



if __name__ == "__main__":
    # Execute main function
    main()

