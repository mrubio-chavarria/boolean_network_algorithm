#!/home/mario/Projects/boolean_2/software/venv/bin/venv python

# Libraries
from graph import Graph
import json



# Functions
def main():
    """
    DESCRIPTION:
    Main method, the guideline to execute the inference.
    """
    # Read problem information
    data = json.load(open("data.json"))
    # Create the object graph
    graph = Graph(**data)
    graph.obtain_pathways_from_graph()
    graph.generate_priority_matrices()
    graph.generate_NCBFs()
    graph.generate_boolean_networks()
    # # Parse the graph
    # parsed_graph = graph_parser(graph)
    # # Create all the possible NCBF functions by node
    # ncbf_by_node = ncbf_tree(parsed_graph)



if __name__ == "__main__":
    # Execute main function
    main()

