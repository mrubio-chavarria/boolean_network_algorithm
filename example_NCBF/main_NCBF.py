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
# from bn_utils import filter_boolean_networks
# from serializers import Serializer


# Functions
def main_NCBF():
    """
    DESCRIPTION:
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
    print("Pairs: ", graph.pathway_groups[0])
    print("--------")
    print("Network: ", graph.pre_networks[0])


def main():
    """
    DESCRIPTION:
    Main method, the guideline to execute the inference.
    """
    # Read the initial data
    # Read the pathways index, the repeat index passed during the execution
    pathways_index = '122'#sys.argv[1]
    pathways_index = int(pathways_index)
    # Read inference information in the data.json file
    data = json.load(open('data.json'))
    load_stored_networks = data['load_stored_networks']     # Assess if it is needed to load a previous session
    partial = data['partial']                               # Parameter to relax the conditions during network filtering
    attractors = data['attractors']                         # Attractors searched in the network
    if data['attractors'] == "None":
        attractors = None
        data['attractors'] = None
    n_attractors = data['n_attractors']                     # Total number of attractors searched in the networks
    if n_attractors == "None":
        n_attractors = None
    delete_repeated = data['delete_repeated']               # Assess whether to delete repeated networks from the stored files
    if data['min_conflicts'] is not None:                   # Set a minimum number of conflicts for the written example, which is randomly chosen
        min_conflicts = data['min_conflicts']
    data['pathways_index'] = pathways_index                 # Store pathway index with the rest of data
    # Delete data that are no longer required
    del data['load_stored_networks']
    del data['n_attractors']
    del data['delete_repeated']
    del data['min_conflicts']
    # Check networks directory
    if not os.path.exists('networks'):
        # Create it if it does not exist
        os.makedirs('networks')
    # Measure the time taken in every step
    t = TicToc()
    t.tic()
    if not load_stored_networks:
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
        # Cross all the priority matrices with all the combinations of NCBF networks and
        # pathway groups to obtain all the inferences to perform. The data to perform
        # each inference is stored in a dict refered to as boolean network. This method
        # creates all these dicts with the data
        graph.generate_boolean_networks()
        t.toc('Computation time: ', restart=True)
        # Solve the conflicts (it takes time) for every combination of:
        # - Priority matrix.
        # - Pathway group.
        # - NCBF network.
        print('Solve network conflicts')
        graph.solve_conflicts()
        t.toc('Computation time: ', restart=True)
        # Select only those networks that has the desired attractors
        print('Select only related networks')
        graph.prefilter()
        t.toc('Computation time: ', restart=True)
        # Compute relevant features of the converging networks
        # In this method, the attractors are computed through the 
        # Tarjan algorithm, a filter based on the attractor number
        # is applied 
        graph.format_network()
        t.toc('Computation time: ', restart=True)
        # Store the resulting networks
        filename = 'networks/networks_' + '_'.join(str(datetime.now()).split(' ')) + '.pickle'
        if graph.boolean_networks:
            print(f'Storing networks in: {filename}')
            with open(filename, 'wb') as file:
                # Transform the network in a format that can be written in 
                # a JSON file
                serialized_networks = [Serializer(net).get_serialized_network() 
                    for net in graph.boolean_networks]
                pickle.dump(serialized_networks, file)
        else:
            print('No network passed the inference!')
        # Prepare for the filtering
        boolean_networks = graph.boolean_networks
        t.toc('Computation time: ')
        print('Networks stored')
    else:
        # Load stored networks
        print('Loading previous networks')
        sessions_folder = f'networks/{data["algorithm"]}/'
        boolean_networks = [pickle.load(open(sessions_folder + file, 'rb'))
            for file in os.listdir(sessions_folder)]
        boolean_networks = [it for sb in boolean_networks for it in sb]
        # Delete repeated
        if delete_repeated:
            # Take only those networks that are not repeated
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
            # Keep the non-repeated and store them
            boolean_networks = non_repeated_networks
            filename = sessions_folder + 'networks_' + '_'.join(str(datetime.now()).split(' ')) + '.pickle'
            with open(filename, 'wb') as file:
                pickle.dump(non_repeated_networks, file)
        t.toc('Computation time: ')
    # Analysis over the graph of the resulting networks
    if boolean_networks:
        # Filtering
        print('Filter the resulting networks')
        boolean_networks = filter_boolean_networks(boolean_networks,
                                                    attractors=attractors,
                                                    n_attractors=n_attractors,
                                                    partial=partial)
        t.toc('Computation time: ')
        # Print a random network
        example_filename = 'example.json'
        try:
            selected_network = choice(list(filter(lambda x: len(x['conflicts']) >= min_conflicts, boolean_networks)))
            with open(example_filename, 'w') as file:
                json.dump(selected_network, file)
            print('Example network printed in:', example_filename)
        except IndexError:
            print('No network has enough conflicts to be printed')
        print('Inference finished')



if __name__ == "__main__":
    # Execute main function
    main_NCBF()
    # main()

