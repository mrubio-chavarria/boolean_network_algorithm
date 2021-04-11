#!/home/mario/Projects/boolean_2/software/venv/bin/venv python

# Libraries
import numpy as np
import itertools

##############################################################
# FUNCTIONS
##############################################################

# REVISAR
def introduce_contradictory(parsed_graph, node, target_node, which='inhibitors'):
    """
    DESCRIPTION:
    A function to introduce contradictory node in the graph. It will copy the relationships from its
    equivalent node and will be introduced as input in its target. The nature of the input is denoted
    through the parameter which.
    :param graph: [dict] the parsed_graph from which the nodes are to be extracted.
    :param node: [str] name of the node to be introduced.
    :param target_node: [str] the node to which the action of the contradictory nodes is directed.
    :param which: [str] the indication to take either the activators or the inhibitors.
    """
    # Obtain parameters for the new node
    new_node = node + '_C'
    activators = obtain_inputs(parsed_graph, -1)
    inhibitors = obtain_inputs(parsed_graph, -1, 'inhibitors')
    act_pos = len(parsed_graph['activators'])
    inh_pos = len(parsed_graph['inhibitors'])
    # Introduce the new node in the graph
    parsed_graph['nodes'].append(new_node)
    parsed_graph['activators'] += activators
    parsed_graph['inhibitors'] += inhibitors
    parsed_graph['positions_activators'].append(act_pos)
    parsed_graph['positions_inhibitors'].append(inh_pos)
    parsed_graph['n_activators'].append(len(activators))
    parsed_graph['n_inhibitors'].append(len(inhibitors))
    # Introduce the new node
    target_position = list(filter(lambda x: x[1] == target_node, enumerate(parsed_graph['nodes'])))[0][0]
    if which == 'inhibitors':
        # Add the node
        position = parsed_graph['positions_inhibitors'][target_position]
        number = parsed_graph['n_inhibitors'][target_position]
        parsed_graph['inhibitors'] = parsed_graph['inhibitors'][0:position] + parsed_graph['inhibitors'][position:position + number] + [new_node] + parsed_graph['inhibitors'][number::]
        # Update parameters    
        parsed_graph['n_inhibitors'][target_position] += 1 
        parsed_graph['positions_inhibitors'][target_position+1::] = list(map(lambda x: x + 1, parsed_graph['positions_inhibitors'][target_position+1::]))
    else:
        # Add the node
        position = parsed_graph['positions_activators'][target_position]
        number = parsed_graph['n_activators'][target_position]
        parsed_graph['activators'] = parsed_graph['activators'][0:position] + parsed_graph['activators'][position:position + number] + [new_node] + parsed_graph['activators'][number::]
        # Update parameters    
        parsed_graph['n_activators'][target_position] += 1 
        parsed_graph['positions_activators'][target_position+1::] = list(map(lambda x: x + 1, parsed_graph['positions_activators'][target_position+1::]))


def delete_inputs(parsed_graph, node_index, input_node, which='activators'):
    """
    DESCRIPTION:
    A function to delete activators or inhibitors from a node in the parsed graph.
    :param parsed_graph: [dict] the graph from which the nodes are to be extracted.
    :param which: [str] the indication to take either the activators or the inhibitors.
    :param node_index: [int] the position of the node in the initial list.
    :param input: [str] the node to delete.    
    :return: [list] list of nodes, either activators or inhibitors depending on the which parameter.
    """
    # Auxiliar functions
    def aux_fun(x):
        if x:
            return x - 1
    # Parameters
    elements = parsed_graph[which]
    n_elements = parsed_graph['n_' + which][node_index]
    position_elements = parsed_graph['positions_' + which][node_index]
    # Substitute the values
    new_elements = list(set(obtain_inputs(parsed_graph, node_index, which)) - set(input_node))
    elements[position_elements:position_elements + n_elements] = new_elements
    # Update parameters    
    # Rest one to the upwards positions
    parsed_graph['positions_' + which][node_index + 1::] = list(map(aux_fun, parsed_graph['positions_' + which][node_index + 1::]))
    parsed_graph['n_' + which][node_index] = parsed_graph['n_' + which][node_index] - 1


def obtain_inputs(parsed_graph, node_index, which='activators'):
    """
    DESCRIPTION:
    A function to access the activators and inhibitors in the parser graph in the same way that they are accessed in 
    the original graph.
    :param parsed_graph: [dict] the graph from which the nodes are to be extracted.
    :param which: [str] the indication to take either the activators or the inhibitors.
    :param node_index: [int] the position of the node in the initial list.    
    """
    elements = parsed_graph[which]
    n_elements = parsed_graph['n_' + which][node_index]
    position_elements = parsed_graph['positions_' + which][node_index]
    if not n_elements:
        return []
    else:
        return elements[position_elements:position_elements + n_elements]


def auxiliar_function(groupA, groupB, n, path=None):
    """
    DESCRIPTION:
    A generator to iterate over a dynamic number of for loops. Here are generated all possible structures given a set of
    groups.
    :param groupA: first set of levels in the structure.
    :param groupB: second set of levels in the structure.
    :param n: total weight of the path.
    :param path: path to propagate, to generate.
    :return: all possible structures, paths.
    """
    # Parameters to assess the path
    n_A = groupA[-1]
    n_B = groupB[-1]
    if path is not None:
        if len(path) % 2 != 0:
            # Last layer A, add B
            for el_B2 in groupB[:]:
                pathB2 = path + [el_B2]
                valB = np.sum(pathB2)
                if valB == n:
                    model = np.array(pathB2)
                    vecA = np.zeros([1, len(model)])
                    vecA[0, 0::2] = 1
                    vecB = np.zeros([1, len(model)])
                    vecB[0, 1::2] = 1
                    if (np.sum(vecB * model) == n_B) and (np.sum(vecA * model) == n_A):
                        yield pathB2
                elif valB < n:
                    yield from auxiliar_function(groupA, groupB, n, path=pathB2)
        else:
            # Last layer B, add A
            for el_A2 in groupA[:]:
                pathA2 = path + [el_A2]
                valA = np.sum(pathA2)
                if valA == n:
                    model = np.array(pathA2)
                    vecA = np.zeros([1, len(model)])
                    vecA[0, 0::2] = 1
                    vecB = np.zeros([1, len(model)])
                    vecB[0, 1::2] = 1
                    if (np.sum(vecB * model) == n_B) and (np.sum(vecA * model) == n_A):
                        yield pathA2
                elif valA < n:
                    yield from auxiliar_function(groupA, groupB, n, path=pathA2)
    else:
        # Start with a layer of A
        for el_A1 in groupA[:]:
            yield from auxiliar_function(groupA, groupB, n, path=[el_A1])


def graph_parser(graph):
    """
    DESCRIPTION:
    A function to check the graph, and to adapt it to more convenient format.
    :param graph: [dict] graph, contains all the information to compute the NCBF.
    :return: [dict] transformed graph, without nested structures beyond the lists. Note that a in the 
    return object, a position only makes sense if n_acitvators (or n_inhibitors) is greater than 0.
    """
    # Initialise parameters
    act_last = 0
    inh_last = 0
    total_activators = []
    total_inhibitors = []
    positions_activators = []
    positions_inhibitors = []
    # Compute the parameters from the graph
    for i in range(len(graph['nodes'])):
        positions_activators.append(len(total_activators))
        positions_inhibitors.append(len(total_inhibitors))
        total_activators += graph['activators'][i]
        total_inhibitors += graph['inhibitors'][i]
    # Build the return object
    parsed_graph = {
        'nodes': graph['nodes'],  # nodes of the graph
        'activators': total_activators,  # list with all the activators
        'inhibitors': total_inhibitors,  # list with all the inhibitors
        'positions_activators': positions_activators,  # positions of activators for each node
        'positions_inhibitors': positions_inhibitors,  # positions of activators for each node
        # store the number of activators and inhibitors by node
        'n_activators': [len(graph['activators'][i]) for i in range(len(graph['nodes']))],
        'n_inhibitors': [len(graph['inhibitors'][i]) for i in range(len(graph['nodes']))]
    }
    return parsed_graph


def ncbf_node(graph, node_position):
    """
    DESCRIPTION:
    This function computes all the possible NCBF for a given node in the graph.
    :param graph: [dict] graph, contains all the information to compute the NCBF.
    :param node_position: [int] the position that denotes the node in the graph.
    """
    # Check if there are contradictory nodes
    node = graph['nodes'][node_position]
    activators = set(obtain_inputs(graph, node_position))
    inhibitors = set(obtain_inputs(graph, node_position, 'inhibitors'))
    contradictory_nodes = activators & inhibitors
    # Delete contradictory nodes from inhibitors
    [delete_inputs(graph, node_position, contradictory_node, 'inhibitors')
        for contradictory_node in contradictory_nodes]
    # Create the new nodes based on the contradictory ones
    print(contradictory_nodes)
    print(graph)
    # REVISAR
    [introduce_contradictory(graph, contradictory_node, node) for contradictory_node in contradictory_nodes]
    print(graph)

    
        


    # # Generate all the structure codes for this node
    # n_activators = len(graph["activators"][node_position])
    # n_inhibitors = len(graph["inhibitors"][node_position])
    # groupA = np.linspace(1, n_activators, n_activators)
    # groupB = np.linspace(1, n_inhibitors, n_inhibitors)
    # n = len(activators) + len(inhibitors)

    # structural_codes = list(auxiliar_function(groupA, groupB, n))
    # print(structural_codes)


def ncbf_tree(graph):
    """
    DESCRIPTION:
    The function to build all the possible ncbf given the a series of activators and inhibitors.
    :param graph: [dict] the representation of the problem graph. For every node in the graph
    are returned its associated NCBF. Note that this is not the original, but the parsed graph.
    :return: [dict] the NCBF computed by node. 
    """
    # Iterate for every node
    # for i in range(len(graph["nodes"])):
    i = 0
    ncbf_node(graph, i)
    
    return {}

#############
# NEW
#############
def ncbf_recursive(group1, group2, n_elements, path=[]):
    """
    DESCRIPTION:
    The algorithm that computes all the possible NCBF for a given node.
    :param group1: [list] possible layers to build a NCBF. Initially this
    group was the activators.
    :param group2: [list] possible layers to build a NCBF. Initially this
    group was the inhibitors.
    :param n_elements: [int] number of nodes that shoudl gather all the 
    layers in every NCBF.
    :param path: [list] the layers, that share all the developed NCBF.
    :return: [list] all the NCBFs obtained from this path with these
    activators and inhibitors.
    """
    # Helper functions
    def kernel(layer):
        """
        DESCRIPTION:
        A function to filter all the layers that should be discarded
        because the nodes that they gather have been used.
        :param layer: [str] the layer that is in test.
        :return: [bool] result of the test.
        """
        forbidden_elements = set(group1[i])
        layer_elements = set(list(layer))
        return not bool(layer_elements & forbidden_elements)

    if group1:
        # Recursive case
        ncbfs = []
        if not group2:
            group1 = [layer for layer in group1
                if len(''.join(path)) +  len(layer) == n_elements]
        for i in range(len(group1)):
            # Filter for incompatible types
            new_group1 = list(filter(kernel, group1))
            ncbfs.append(
                ncbf_recursive(group2, new_group1, n_elements, path + [group1[i]])
                )
        ncbfs = [it for sb in ncbfs for it in sb]
    else:
        # Base case
        if group2:
            ncbfs = [path + [layer] for layer in group2
                if len(''.join(path)) +  len(layer) == n_elements]
        else:
            ncbfs = [path]
    return ncbfs


def ncbf_obtain_domain(structure, info, space, first=False):
    """
    DESCRIPTION:
    Recursive algorithm to obtain the domain of a NCBF given the layer 
    structure and the info about the nodes in the layers.
    :param structure: [list] strings that represent the nodes in every layer
    of the NCBF.
    :param info: [dict] the needed information of every node, its domain and
    if it is an activator or not.
    :param first: [bool] a variable that indicates if this instantiation of
    the function is the first or it is a nested one.
    :param space: [set] all the possible terms with the number nodes studied.
    """
    if not structure:
        # Base case
        return space
    else:
        # Recursive case
        # Check for layers with several nodes
        if len(structure[0]) == 1:
            domain_current_layer = space - info[structure[0]][0]
        else:
            domain_current_layer = set().union(*[space - info[node][0] for node in structure[0]])
        domain_downwards_layers = ncbf_obtain_domain(structure[1::], info, space)
        domain = domain_current_layer & domain_downwards_layers
        domain = space - domain if (first and info[structure[0][0]][1]) else domain
        return domain


def ncbf_generator(activators, inhibitors, space):
    """
    DESCRIPTION:
    A function to generate all the NCBF, and handle the algorithm selection 
    and the situation of the contradictory nodes.
    :param activators: [list] activator pathways targeting the selected node.
    :param inhibitors: [list] inhibitor pathways targeting the selected node.
    :param space: [set] all the possible terms with the number nodes studied.
    :return: [list] the domain (set) of every NCBF.
    """
    # Detect contradictory nodes
    activator_nodes = set(pathway['antecedent'] for pathway in activators)
    inhibitor_nodes = set(pathway['antecedent'] for pathway in inhibitors)
    contradictory_nodes = activator_nodes & inhibitor_nodes
    # Assess if there is contradictory nodes and replace then
    if contradictory_nodes:
        # Modify the inhibitors from contradictory nodes
        contradictory_inhibitors = []
        non_contradictory_inhibitors = []
        [
            contradictory_inhibitors.append(pathway)
            if pathway['antecedent'] in contradictory_nodes else
            non_contradictory_inhibitors.append(pathway)
            for pathway in inhibitors
        ]
        relation_original_modified = {}
        modified_inhibitors = []
        # Any modified inhibitor will show a numeric antecedent
        id_modified = 0
        modified_nodes = set()
        for node in contradictory_nodes:
            id_modified += 1
            original_pathway = list(filter(lambda pathway: pathway['antecedent'] == node, contradictory_inhibitors))[0]
            modified_pathway = {'antecedent': str(id_modified),
                                'consequent': original_pathway['consequent'],
                                'activator': False,
                                'domain': original_pathway['domain']}
            relation_original_modified[node] = original_pathway
            modified_inhibitors.append(modified_pathway)
            modified_nodes.add(str(id_modified))
        # Replace the contradictory inhibitors by the modified ones
        inhibitors = non_contradictory_inhibitors + modified_inhibitors
        inhibitor_nodes = inhibitor_nodes - contradictory_nodes | modified_nodes
    # Execute the inference algorithm
    activator_possibilities = [itertools.combinations(activator_nodes, i + 1)
        for i in range(len(activator_nodes))]
    activator_possibilities = [''.join(sorted(it)) for sb in activator_possibilities for it in sb]
    inhibitor_possibilities = [itertools.combinations(inhibitor_nodes, i + 1) 
        for i in range(len(inhibitor_nodes))]
    inhibitor_possibilities = [''.join(sorted(it)) for sb in inhibitor_possibilities for it in sb]
    n_elements = len(list(activator_nodes | inhibitor_nodes))
    ncbfs_1 = ncbf_recursive(activator_possibilities, inhibitor_possibilities, n_elements)
    ncbfs_2 = ncbf_recursive(inhibitor_possibilities, activator_possibilities, n_elements)
    ncbfs = ncbfs_1 + ncbfs_2
    # Create relationship between antecedent and pathway
    antecedent_info = {pathway['antecedent']: (pathway['domain'], pathway['activator'])
        for pathway in activators + inhibitors}
    # Obtain the domain of every NCBF and return
    # CHECK THIS FUNCTION
    return [ncbf_obtain_domain(ncbf, antecedent_info, space, first=True) for ncbf in ncbfs]