#!/home/mario/Projects/boolean_2/software/venv/bin/venv python

# Libraries
import json
from pytictoc import TicToc
from graph import Graph
import pickle
from datetime import datetime


# Functions
def main():
    """
    DESCRIPTION:
    Main method, the guideline to execute the inference.
    """
    # Assess if it is needed to load a previous session
    load_session = False
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
        # Perform inference (takes time)
        print('Solve the conflicts of all the networks')
        graph.solve_conflicts()
        # Store the graph
        filename = 'session_' + '_'.join(str(datetime.now()).split(' ')) + '.pickle'
        with open(filename, 'wb') as file:
            pickle.dump(graph, file)
    else:
        # Load a previous session (graph)
        filename = 'session_2021-05-01_08:16:16.519277.pickle'
        with open(filename, 'rb') as file:
            graph = pickle.load(file)
    # Analysis over the graph
    print('Filter the resulting networks')
    graph.filter_boolean_networks()
    print()



if __name__ == "__main__":
    # Execute main function
    main()

