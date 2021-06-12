#!/home/mario/Projects/boolean_2/software/venv/bin/venv python

# Libraries
import numpy as np
import itertools
from string import ascii_letters, digits
from random import choice


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
        # Compute the terms in this layer
        current_layer = structure[0]
        layer_domain = space
        for factor in current_layer:
            factor_domain = space - info[factor][0]
            layer_domain = layer_domain & factor_domain
        # Compute the terms in downward layers
        downward_domain = ncbf_obtain_domain(structure[1::], info, space)
        layer_domain = layer_domain & downward_domain
        layer_domain = space - layer_domain
        # Check for the canalised value of the outter layer
        if first:
            if not info[structure[0][0]][1]:
                layer_domain = space - layer_domain
        return layer_domain
        

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
    # Modify the repeated nodes to build the NCBF
    available_variables = set(ascii_letters + digits) - inhibitor_nodes - activator_nodes
    # Activators
    antecedents = []
    repeated_activators = []
    for i in range(len(activators)):
        if activators[i]['antecedent'] not in antecedents:
            antecedents.append(activators[i]['antecedent'])
        else:
            repeated_activators.append(activators[i])
    for repeated in repeated_activators:
        new_symbol = choice(list(available_variables))
        available_variables = available_variables - set(new_symbol)
        repeated['antecedent'] = new_symbol
        activator_nodes = activator_nodes | set(new_symbol)
    # Inhibitors
    antecedents = []
    repeated_inhibitors = []
    for i in range(len(inhibitors)):
        if inhibitors[i]['antecedent'] not in antecedents:
            antecedents.append(inhibitors[i]['antecedent'])
        else:
            repeated_inhibitors.append(inhibitors[i])
    for repeated in repeated_inhibitors:
        new_symbol = choice(list(available_variables))
        available_variables = available_variables - set(new_symbol)
        repeated['antecedent'] = new_symbol
        inhibitor_nodes = inhibitor_nodes | set(new_symbol)
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
    ncbfs = ncbfs_1 if not activator_nodes or not inhibitor_nodes else ncbfs_1 + ncbfs_2
    # Create relationship between antecedent and pathway
    antecedent_info = {pathway['antecedent']: (pathway['domain'], pathway['activator'])
        for pathway in activators + inhibitors}
    # Obtain the domain of every NCBF and return
    domains = [ncbf_obtain_domain(ncbf, antecedent_info, space, first=True) for ncbf in ncbfs]
    return domains