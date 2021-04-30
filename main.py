#!/home/mario/Projects/boolean_2/software/venv/bin/venv python

# Libraries
import json
from pytictoc import TicToc
from graph import Graph


# Functions
def main():
    """
    DESCRIPTION:
    Main method, the guideline to execute the inference.
    """
    # Read problem information
    data = json.load(open('data.json'))
    # Create the object graph
    t = TicToc()
    graph = Graph(**data)
    graph.obtain_pathways_from_graph()
    graph.generate_priority_matrices()
    graph.generate_NCBFs()
    graph.generate_boolean_networks()
    print('Solve the conflicts of all the networks')
    t.tic()
    graph.solve_conflicts()
    t.toc()
    print('Filter the resulting networks')
    graph.filter_boolean_networks()
    print()
    # # Parse the graph
    # parsed_graph = graph_parser(graph)
    # # Create all the possible NCBF functions by node
    # ncbf_by_node = ncbf_tree(parsed_graph)



if __name__ == "__main__":
    # Execute main function
    main()

