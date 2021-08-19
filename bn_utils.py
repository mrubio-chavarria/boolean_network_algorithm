#!/home/mario/Projects/boolean_2/software/venv/bin/venv python

# Libraries
import itertools
from random import choice, sample
from string import ascii_uppercase, digits
from sympy import SOPform
from PyBoolNet.Attractors import compute_attractors_tarjan
from PyBoolNet.FileExchange import bnet2primes
from PyBoolNet.StateTransitionGraphs import primes2stg
from quine_mccluskey.qm import QuineMcCluskey
from exceptions import NoSolutionException, InputModificationException


# Functions
def left_zfill(word, n_digits):
    """
    DESCRIPTION:
    Given a string a number, the function introduces zeros by the left side 
    until the string the size specified in the number.
    :param word: [string] the string that is to be extended. 
    :param n_digits: [int] the number of digits of the final word. 
    :return: [string] extended word. 
    """
    return '0' * (n_digits - len(word)) + word
    

def function_simplifier(node, nodes, n_nodes, terms, method='Quine-McCluskey'):
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
    if not terms:
        raise NoSolutionException()
    if method == 'Quine-McCluskey':
        qm = QuineMcCluskey(use_xor=False)
        minterms = [int(term, 2) for term in terms]
        # Special cases
        if minterms == [0]:
            results = [left_zfill(str(minterm), n_nodes) for minterm in minterms]
        # General case
        else:
            results = [left_zfill(minterm.replace('-', '*'), n_nodes) for minterm
                in qm.simplify(minterms)]
    elif method == 'SymPy':
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


def solver(activator, inhibitor, network, graph_space, checklist):
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
    :param checklist: [set] a set to check for repeated conflicts.
    :return: [list] new pathways obtained from the inference (same format as
    the others).
    """
    # Detect the conflicting region (psi)
    psi = activator['domain'] & inhibitor['domain']
    # If there is no conflicting region 
    if not psi:
        # Apply both pathways to the network
        network['network'][activator['consequent']] = (network['network'][activator['consequent']] | activator['domain']) - inhibitor['domain']
        # And return no conflict
        return []
    # Obtain the priority
    activator_priority = network['priority']['activators'][activator['antecedent']][activator['consequent']]
    inhibitor_priority = network['priority']['inhibitors'][inhibitor['antecedent']][inhibitor['consequent']]
    # Store the pathways before any modifications
    stored_activator = str(activator)
    activator_domain = activator['domain'].copy()
    activator_antecedent = activator['antecedent']
    stored_inhibitor = str(inhibitor)
    inhibitor_domain = inhibitor['domain'].copy()
    inhibitor_antecedent = inhibitor['antecedent']
    # Compare the priority of both activator and inhibitor
    if activator_priority >= inhibitor_priority:
        prioritised = 'activator'
        # Prioritised: activator
        # Non prioritised: inhibitor
        # Check whether the non-prioritised is an input
        if inhibitor['antecedent'] in network['input_nodes']:
            raise InputModificationException
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
            {'id': ''.join(sample(ascii_uppercase + digits, 5)),
            'antecedent': inhibitor['antecedent'], 'consequent': node, 'activator': canalised_value, 'domain': psi}
            for node, canalised_value in target_nodes
        ]
    else:
        prioritised = 'inhibitor'
        # Prioritised: inhibitor
        # Non prioritised: activator
        # Check whether the non-prioritised is an input
        if activator['antecedent'] in network['input_nodes']:
            raise InputModificationException
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
            {'id': ''.join(sample(ascii_uppercase + digits, 5)),
            'antecedent': activator['antecedent'], 'consequent': node, 'activator': canalised_value, 'domain': psi}
            for node, canalised_value in target_nodes
        ]
    # Check for repeated conflicts
    n_initial_conflicts = len(list(checklist))
    checklist |= set([(activator_domain, activator_antecedent, inhibitor_domain, inhibitor_antecedent)])
    n_final_conflicts = len(list(checklist))
    if n_initial_conflicts == n_final_conflicts or not solution_pathways:
        # If there are repeated conflicts or there is no new pathway, there is no solution
        raise NoSolutionException
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
    # Set the type of algorithm
    if network['algorithm'] == 'original':
        while iteration < network['max_iterations']:
            # Go through all the nodes
            for node in network['nodes']:
                new_pathways = True
                checklist = set()
                try:
                    while new_pathways:
                        # Filter pathways with domain
                        pathways['checked'][node]['activators'] = list(filter(lambda pathway: bool(pathway['domain']),
                            pathways['checked'][node]['activators']))
                        pathways['checked'][node]['inhibitors'] = list(filter(lambda pathway: bool(pathway['domain']),
                            pathways['checked'][node]['inhibitors']))
                        # Obtain all the possible pathway pairs
                        pair_group1 = list(itertools.product(pathways['non_checked'][node]['activators'],
                                                            pathways['non_checked'][node]['inhibitors']))
                        pair_group2 = list(itertools.product(pathways['non_checked'][node]['activators'],
                                                            pathways['checked'][node]['inhibitors']))
                        pair_group3 = list(itertools.product(pathways['checked'][node]['activators'],
                                                            pathways['non_checked'][node]['inhibitors']))
                        # Apply pathways when there are no pairs, when there is pairs they are applied downwards
                        new_activators = bool(pathways['non_checked'][node]['activators'])
                        new_inhibitors = bool(pathways['non_checked'][node]['inhibitors'])
                        if not (pair_group1 + pair_group2 + pair_group3) and (new_activators or new_inhibitors):
                            non_applied_pathways = pathways['non_checked'][node]['activators'] + pathways['non_checked'][node]['inhibitors']
                            for nap in non_applied_pathways:  # nap stands for non applied pathway, used for brevity
                                if nap['activator']:
                                    # Activator
                                    network['network'][nap['consequent']] = network['network'][nap['consequent']] | nap['domain']
                                else:
                                    # Inhibitor
                                    network['network'][nap['consequent']] = network['network'][nap['consequent']] - nap['domain'] 
                        # Change checked by non-checked
                        pathways['checked'][node]['activators'].extend(pathways['non_checked'][node]['activators'])
                        pathways['checked'][node]['inhibitors'].extend(pathways['non_checked'][node]['inhibitors'])
                        pathways['non_checked'][node] = {'activators': [], 'inhibitors': []}
                        # Solve all the pairs
                        solution_pathways = [it for sb in 
                            [solver(activator, inhibitor, network, network['graph_space'], checklist) 
                                for activator, inhibitor  in pair_group1 + pair_group2 + pair_group3]
                            for it in sb]
                        # Update the presence of new pathways
                        new_pathways = bool(solution_pathways)
                        # Introduce the new pathways with the previous ones
                        pathways['non_checked'] = {network_node: {
                            'activators': pathways['non_checked'][network_node]['activators'] +
                                list(filter(lambda pathway: (pathway['consequent'] == network_node) and pathway['activator'], solution_pathways)),
                            'inhibitors': pathways['non_checked'][network_node]['inhibitors'] +
                                list(filter(lambda pathway: (pathway['consequent'] == network_node) and not pathway['activator'], solution_pathways))
                            } for network_node in network['nodes']}
                        # Update the node conditions for the outer loop
                        for network_node in network['nodes']:
                            new_node_pathways = pathways['non_checked'][network_node]['activators'] + \
                                pathways['non_checked'][network_node]['inhibitors']
                            node_conditions[network_node] = not new_node_pathways
                except (NoSolutionException, InputModificationException):
                    # Case in which a pathway could not be solved or the solution requires
                    # to modify the input. There is no solution
                    return None
            # Update iteration
            iteration += 1
            # Assess break condition
            if all(node_conditions.values()):
                break
    elif network['algorithm'] == 'new':
        # CHECK THESE ALGORITHM WITH THE NEW MODIFICATIONS
        while iteration < network['max_iterations']:
            # Go through all the nodes
            for node in network['nodes']:
                new_pathways = True
                checklist = set()
                try:
                    # Filter pathways with domain
                    pathways['checked'][node]['activators'] = list(filter(lambda pathway: bool(pathway['domain']),
                        pathways['checked'][node]['activators']))
                    pathways['checked'][node]['inhibitors'] = list(filter(lambda pathway: bool(pathway['domain']),
                        pathways['checked'][node]['inhibitors']))
                    # Obtain all the possible pathway pairs
                    pair_group1 = list(itertools.product(pathways['non_checked'][node]['activators'],
                                                        pathways['non_checked'][node]['inhibitors']))
                    pair_group2 = list(itertools.product(pathways['non_checked'][node]['activators'],
                                                        pathways['checked'][node]['inhibitors']))
                    pair_group3 = list(itertools.product(pathways['checked'][node]['activators'],
                                                        pathways['non_checked'][node]['inhibitors']))
                    # Apply pathways when there are no pairs, when there is pairs they are applied downwards
                    new_activators = bool(pathways['non_checked'][node]['activators'])
                    new_inhibitors = bool(pathways['non_checked'][node]['inhibitors'])
                    if not (pair_group1 + pair_group2 + pair_group3) and (new_activators or new_inhibitors):
                        non_applied_pathways = pathways['non_checked'][node]['activators'] + pathways['non_checked'][node]['inhibitors']
                        for nap in non_applied_pathways:  # nap stands for non applied pathway, used for brevity
                            if nap['activator']:
                                # Activator
                                network['network'][nap['consequent']] = network['network'][nap['consequent']] | nap['domain']
                            else:
                                # Inhibitor
                                network['network'][nap['consequent']] = network['network'][nap['consequent']] - nap['domain']
                    # Change checked by non-checked
                    pathways['checked'][node]['activators'].extend(pathways['non_checked'][node]['activators'])
                    pathways['checked'][node]['inhibitors'].extend(pathways['non_checked'][node]['inhibitors'])
                    pathways['non_checked'][node] = {'activators': [], 'inhibitors': []}
                    # Solve all the pairs
                    solution_pathways = [it for sb in 
                        [solver(activator, inhibitor, network, network['graph_space'], checklist) 
                            for activator, inhibitor  in pair_group1 + pair_group2 + pair_group3]
                        for it in sb]
                    # Update the presence of new pathways
                    new_pathways = bool(solution_pathways)
                    node_conditions[node] = not new_pathways
                    # Introduce the new pathways with the previous ones
                    pathways['non_checked'] = {node: {
                        'activators': pathways['non_checked'][node]['activators'] +
                            list(filter(lambda pathway: (pathway['consequent'] == node) and pathway['activator'], solution_pathways)),
                        'inhibitors': pathways['non_checked'][node]['inhibitors'] +
                            list(filter(lambda pathway: (pathway['consequent'] == node) and not pathway['activator'], solution_pathways))
                        } for node in network['nodes']}
                    # Update the node conditions for the outer loop
                    for network_node in network['nodes']:
                        new_node_pathways = pathways['non_checked'][network_node]['activators'] + \
                            pathways['non_checked'][network_node]['inhibitors']
                        node_conditions[network_node] = not new_node_pathways
                except NoSolutionException:
                    # Case in which a pathway could not be solved. There is no solution
                    return None
            # Update iteration
            iteration += 1
            # Assess break condition
            if all(node_conditions.values()):
                break
    else:
        raise KeyError('Invalid algorithm selected')
    # Assess condition of max iterations. Only the converging networks
    # are of interest
    if all(node_conditions.values()):
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
        # Return the result
        return network


def minterms2bnet(variables, minterms):
    """
    DESCRIPTION:
    A function to obtain the function expression from the minterms in boolnet format.
    It is supposed that all the nodes are in alphabetical order. Important, this 
    function only returns the SOP without the name of the function, what is at the 
    left side of the comma.
    :param nodes: [tuple] the function variables.
    :param minterms: [frozenset] the minterms to build the expression.
    :return: [str] the function expression in boolnet format.
    """
    # When the function is always false
    if not minterms:
        return '0'
    # Check when the minterms are the whole space
    if len(minterms) == 2 ** len(variables):
        return '1'
    # Generañ case
    qm = QuineMcCluskey(use_xor=False)
    simplified_minterms = qm.simplify([int(term, 2) for term in minterms])
    simplified_expression = [left_zfill(minterm.replace('-', '*'), len(variables)) for minterm in simplified_minterms]
    # Pass the expression to boolnet format
    n_variables = range(len(variables))
    bnet_expression = ' | '.join([
        '&'.join([variables[i] if int(minterm[i]) else '!' + variables[i] for i in n_variables if minterm[i] != '*'])
        for minterm in simplified_expression
    ])
    return bnet_expression


def network_formatter(network, min_attractors=2, max_attractors=4):
    """
    DESCRIPTION:
    A function to compute the attractors and format a given network to prepare it for 
    the storage.
    :param network: [dict] the boolean network whose attractors are going to be 
    computed.
    :param min/max_attractors: [int] a filtering parameters. The only networks of interest
    are those with a number of attractors within the limits.
    :return: [dict] the boolean network with the attractors stored. 
    """
    # Write every net function in BoolNet format
    network['network_terms'] = network['network'].copy()
    network['network'] = {
        node: f"{node}, {minterms2bnet(network['nodes'], network['network'][node])}" 
        for node in network['nodes']
    }
    # Obtain the state transition graph
    primes = bnet2primes('\n'.join(network['network'].values()))
    stg = primes2stg(primes, "synchronous")
    # Obtain the attractors
    steady, cyclic = compute_attractors_tarjan(stg)
    network['attractors'] = {'steady': steady, 'cyclic': cyclic}
    # Return only if the attractor condition is met
    if min_attractors <= len(steady + cyclic) and len(steady + cyclic) <= max_attractors:
        # Obtain the initial expression
        network['pre_network'] = {
            node: f"{node}, {minterms2bnet(network['nodes'], network['pre_network'][node])}" 
            for node in network['nodes']
        }
        # Simplify both initial and final expressions
        # Delete not needed variables
        del network['max_iterations']
        del network['graph_space']
        del network['nodes']
        del network['priority']
        del network['network_terms']
        # Return
        return network


def filter_boolean_networks(boolean_networks, attractors=None, n_attractors=None, partial=False):
    """
    DESCRIPTION:
    A function to filter the boolean networks according to different 
    criteria based on the attractors.
    :param boolean_networks: [list] the boolean networks to filter.
    :param attractors: [dict] the attractors that we want the networks
    to show. Format: {'steady': [...], 'cyclic': [...]}
    :param n_attractors: [int] the number of attractors we want the
    networks to show, steady + cyclic.
    :param partial: [bool] a flag to indicate if it is valid with showing
    just one of the attractors and a maximum of n_attractors, but less is
    permitted.
    :return: [list] the boolean networks that meet the criteria.
    """
    # Kernels
    def by_attractor(network):
        """
        DESCRIPTION:
        A function to filter by attractor. It is returned true only if
        the network has the attractors specified in attractors.
        :param network: [dict] the boolean network to filter.
        :return: [bool] value to indicate if the network meet the 
        criterium.
        """
        # Obtain all the attractors from the network
        network_attractors = network['attractors']['steady'] +\
            network['attractors']['cyclic']
        # If partial, look only for one attractor.
        if partial:
            condition = any([attractor in network_attractors 
                for attractor in attractors])
        else:
            condition = all([attractor in network_attractors 
                for attractor in attractors])
        return condition

    def by_n_attractor(network):
        """
        DESCRIPTION:
        A function to filter by the number of attractors specified.
        :param network: [dict] the boolean network to filter.
        :return: [bool] value to indicate if the network meet the 
        criterium.
        """
        # Obtain the number of attractors in the network
        n_network_attractors = len(network['attractors']['steady'] +\
            network['attractors']['cyclic'])
        # If partial, look only for networks that at most have the 
        # specified number of attractors
        condition = (n_network_attractors <= n_attractors 
            if partial else n_network_attractors == n_attractors)
        return condition

    # Filtering
    if attractors is not None:
        boolean_networks = list(filter(by_attractor, boolean_networks))
    if n_attractors is not None:
        boolean_networks = list(filter(by_n_attractor, boolean_networks))
    return boolean_networks


def prefilter_by_attractor(networks, attractors, partial=True):
    """
    DESCRIPTION:
    A function to filter boolean networks based on if they hold certain
    attractors or not.
    :param networks: [list] boolean networks (dict) to filter based on
    their attractors.
    :param attractors: [list] attractors (str) to filter the networks.
    :param partial: [bool]flag to indicate if the condition should be 
    relaxed. If partial == True, then it will be enough for the network
    to show just one attractor from the list to pass the filter. If False
    the network will need to show all the attractors in the list.
    """
    # Helper functions
    def check_node(node_index, nodes, attractor, network):
        node = nodes[node_index]
        # Check if the attractor value is 1 or 0
        if int(attractor[node_index]):
            result = attractor in network['network'][node]
        else:
            result = attractor not in network['network'][node]
        return result

    # Iterate over every network
    nodes = list(networks[0]['network'].keys())
    for network in networks:
        # Network might be a None resulting from unsuccessful conflict solving
        if network is not None:
            attractor_conditions = []
            # Check every attractor
            for attractor in attractors:
                node_conditions = [check_node(node_index, nodes, attractor, network) 
                    for node_index in range(len(nodes))]
                attractor_conditions.append(all(node_conditions))
            if all(attractor_conditions) or partial and any(attractor_conditions):
                yield network