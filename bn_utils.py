#!/home/mario/Projects/boolean_2/software/venv/bin/venv python

# Libraries
import itertools
from random import choice, sample
from string import ascii_uppercase, digits
from sympy import SOPform


# Functions
def function_simplifier(node, nodes, n_nodes, terms, method='SymPy'):
    """
    DESCRIPTION:
    A function to simplify the terms of a given boolean function.
    :param node: [str] node that denotes the boolean function that we
    are simplifying.
    :param nodes: [list] variables of the function to simplify.
    :param n_nodes: [int] the number of nodes computed to save time.
    :param terms: [frozenset] strings the represent the terms.
    :param method: [str] variable to indicate the simplification method to
    use.
    :return: [list] minterms of the simplfied function.
    """
    if method == 'SymPy':
        terms = [[int(digit) for digit in list(term)] for term in terms]
        simplified = str(SOPform(nodes, terms))
        terms = [
            term.replace('(', '').replace(')', '').replace(' ', '').split('&')
            for term in simplified.split('|')
        ]
        results = []
        for term in terms:
            result = dict(zip(nodes, ['*'] * n_nodes))
            for factor in term:
                result[factor.replace('~', '')] = '0' if '~' in factor else '1'
            results.append(''.join(result.values()))  
    else:
        raise AttributeError('Introduced non-valid simplification mode')
    return results


def solver(activator, inhibitor, network, graph_space):
    """
    DESCRIPTION:
    A helper function that will just solve the conflicts between the pathways
    of a given pair.
    :param activator: [dict] pathway that provokes a 1 in the consequent.
    :param inhibitor: [dict] pathway that provokes a 0 in the consequent.
    :param network: [dict] all the information of the network to solve the 
    conflicts.
    :param graph_space: [set] all the terms that build the mathematical space
    in which is developed.
    :return: [list] new pathways obtained from the inference (same format as
    the others).
    """
    # Detect the conflicting region (psi)
    psi = activator['domain'] & inhibitor['domain']
    # If there is no conflicting region no pathways are returned
    if not psi:
        return []
    # Obtain the priority
    activator_priority = network['priority']['activators'][activator['antecedent']][activator['consequent']]
    inhibitor_priority = network['priority']['inhibitors'][inhibitor['antecedent']][inhibitor['consequent']]
    # Store the pathways before any modifications
    stored_activator = str(activator)
    stored_inhibitor = str(inhibitor)
    # Compare the priority of both activator and inhibitor
    if activator_priority >= inhibitor_priority:
        prioritised = 'activator'
        # Prioritised: activator
        # Non prioritised: inhibitor
        # Modify the non-prioritised pathway
        inhibitor['antecedent'] = ''.join(sorted(set(activator['antecedent'] + inhibitor['antecedent'])))
        inhibitor['domain'] = inhibitor['domain'] - activator['domain']
        # Modify the consequent node function
        network['network'][activator['consequent']] = (network['network'][activator['consequent']] | activator['domain']) - inhibitor['domain']
        # Obtain the solution pathway
        simplified_minterms = function_simplifier(activator['consequent'], network['network'].keys(), len(network['network'].keys()), graph_space - network['network'][activator['consequent']])
        # It is just needed to meet the condition imposed by one term
        # but there is a new pathway per node obtained in the solution
        target_nodes = [(node, bool(int(canalised_value))) for node, canalised_value in zip(network['network'].keys(), list(str(choice(simplified_minterms)))) if canalised_value != '*']
        # Generate the solution pathway
        solution_pathways = [
            {'id': ''.join(sample(ascii_uppercase + digits, 4)),
            'antecedent': inhibitor['antecedent'], 'consequent': node, 'activator': canalised_value, 'domain': psi}
            for node, canalised_value in target_nodes
        ]
    else:
        prioritised = 'inhibitor'
        # Prioritised: inhibitor
        # Non prioritised: activator
        # Modify the non-prioritised pathway
        activator['antecedent'] = ''.join(sorted(set(activator['antecedent'] + inhibitor['antecedent'])))
        activator['domain'] = activator['domain'] - inhibitor['domain']
        # Modify the consequent node function
        network['network'][activator['consequent']] = (network['network'][activator['consequent']] | activator['domain']) - inhibitor['domain']
        # Obtain the solution pathway
        simplified_minterms = function_simplifier(activator['consequent'], network['network'].keys(), len(network['network'].keys()), network['network'][activator['consequent']])
        # It is just needed to meet the condition imposed by one term
        # but there is a new pathway per node obtained in the solution
        target_nodes = [(node, bool(int(canalised_value))) for node, canalised_value in zip(network['network'].keys(), list(str(choice(simplified_minterms)))) if canalised_value != '*']
        # Generate the solution pathway
        solution_pathways = [
            {'id': ''.join(sample(ascii_uppercase + digits, 4)),
            'antecedent': activator['antecedent'], 'consequent': node, 'activator': canalised_value, 'domain': psi}
            for node, canalised_value in target_nodes
        ]
    # Store the conflict with its solution. Store the string
    # to keep the current state of the pathways
    network['conflicts'].append(
        {'activator': stored_activator, 'inhibitor': stored_inhibitor, 'solution': str(solution_pathways), 'prioritised': prioritised}
    )
    return solution_pathways


def conflicts_solver(network):
    """
    DESCRIPTION:
    The function to solve the conflicts of a given network. It just modified the
    network according to the conflicts solving.
    :param network: [dict] all the information needed to compute the inference
    of a network.
    """
    # Create network and pathways
    network['network'] = dict(network['pre_network'])
    network['pathways'] = {
        node: {
        'activators': [dict(pathway) for pathway in network['pre_pathways'][node]['activators']],
        'inhibitors': [dict(pathway) for pathway in network['pre_pathways'][node]['inhibitors']]
        }
        for node in network['nodes']
    }
    # Solve the conflicts
    # Divide the pathways in those checked and non-checked. If needed, probably
    # this could be improved
    pathways = {
        'non_checked': dict(network['pathways']),
        'checked': {node: {'activators': [], 'inhibitors': []} for node in network['nodes']}
    }
    # Iterate a maximum of times
    iteration = 0
    node_conditions = {node: False for node in network['nodes']}
    while iteration < network['max_iterations']:
        # Go through all the nodes
        for node in network['nodes']:
            # Obtain all the possible pathway pairs
            pair_group1 = list(itertools.product(pathways['non_checked'][node]['activators'],
                                                pathways['non_checked'][node]['inhibitors']))
            pair_group2 = list(itertools.product(pathways['non_checked'][node]['activators'],
                                                pathways['checked'][node]['inhibitors']))
            pair_group3 = list(itertools.product(pathways['checked'][node]['activators'],
                                                pathways['non_checked'][node]['inhibitors']))
            # Change checked by non-checked
            pathways['checked'][node]['activators'].extend(pathways['non_checked'][node]['activators'])
            pathways['checked'][node]['inhibitors'].extend(pathways['non_checked'][node]['inhibitors'])
            pathways['non_checked'][node] = {'activators': [], 'inhibitors': []}
            # Solve all the pairs 
            solution_pathways = [it for sb in 
                [solver(pair[0], pair[1], network, network['graph_space']) for pair in pair_group1 + pair_group2 + pair_group3]
                for it in sb]
            node_conditions[node] = not solution_pathways
            # Introduce the new pathways with the previous ones
            pathways['non_checked'] = {node: {
                'activators': pathways['non_checked'][node]['activators'] +
                    list(filter(lambda pathway: (pathway['consequent'] == node) and pathway['activator'], solution_pathways)),
                'inhibitors': pathways['non_checked'][node]['inhibitors'] +
                    list(filter(lambda pathway: (pathway['consequent'] == node) and not pathway['activator'], solution_pathways))
                } for node in network['nodes']}
        # Update iteration
        iteration += 1
        # Assess break condition
        if all(node_conditions.values()):
            break
    # Assess condition of max iterations
    network['converge'] = all(node_conditions.values())
    # Save new pathways
    network['pathways'] = pathways['checked']
    # Impose all the pathways in the network
    for node in network['nodes']:
        # Activators
        for pathway in network['pathways'][node]['activators']:
            network['network'][node] = network['network'][node] | pathway['domain']
        # Inhibitors
        for pathway in network['pathways'][node]['inhibitors']:
            network['network'][node] = network['network'][node] - pathway['domain']

