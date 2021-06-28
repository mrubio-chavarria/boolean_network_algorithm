#!/home/mario/Projects/boolean_2/software/venv/bin/venv python

# Libraries
import itertools
from random import sample
from multiprocessing import Pool, cpu_count
from string import ascii_uppercase, digits
from sympy.logic import SOPform
from sympy import symbols
import tqdm
from uuid import uuid1
from ncbf_utils import ncbf_generator
from bn_utils import conflicts_solver, network_formatter, prefilter_by_attractor
from pytictoc import TicToc


# Classes
class Graph:

    # Methods
    def __init__(self, nodes, activators, inhibitors, n_simulations, multiprocess, n_free_cores, max_iterations, algorithm, attractors):
        """
        DESCRIPTION:
        The constructor of the Graph object. All the network inference
        is based on this object.
        :param nodes: [list] strings that represent the nodes that build the 
        graph.
        :param activators: [dict] dict in which the keys are the nodes, and the
        values their activators.
        :param inhibitors: [dict] the same but for the inhibitors.
        :param n_simulations: [int] the number of times that every boolean 
        network should be infered. It is exactly the number of priority matrices
        computed.
        :param multiprocess: [bool] parameter decide whether the computations 
        should be done with multiple processes or not.
        :param n_free_cores: [int] for the multiprocessing, to specify how much 
        computer capacity wants the user to let free.
        :param max_iterations: [int] maximum number of iterations around the 
        nodes to solve the conflicts in every given network.
        :param algorithm: [str] code to describe the algorithm used to solve the
        conflicts: new or original.
        :param attractors: [list] a list with the searched attractors in str 
        format.
        """
        # Always, the nodes are ordered alphabetically
        self.nodes = tuple(sorted(nodes))
        self.activators = {key: tuple(sorted(activators[key])) 
            for key in activators.keys()}
        self.inhibitors = {key: tuple(sorted(inhibitors[key])) 
            for key in inhibitors.keys()}
        self.min_attractors = 1
        self.max_attractors = 4
        self.n_simulations = n_simulations
        self.n_nodes = len(nodes)
        self.multiprocess = multiprocess
        self.used_cores = cpu_count() - n_free_cores
        self.max_iterations = max_iterations
        self.algorithm = algorithm
        self.attractors = attractors
        # Generate all the possible minterms in a space of len(nodes) variables.
        # IMPORTANT: the node position in every term is alphabetical: A:0, B:1...
        self.graph_space = frozenset('{:0{}b}'.format(i, self.n_nodes) 
            for i in range(2 ** self.n_nodes))
    
    def obtain_pathways_from_graph(self):
        """
        DESCRIPTION:
        The method to obtain all the possible groups of pathways from the graph 
        based on the pairs of canalising/canalised values described in the 
        article. The groups are added to the Graph object. Every pathway is a 
        tuple in which the fields have the following meaning:
        1. Antecedent: [str] the row key to obtain the priority in the 
        matrix.
        2. Consequent: [str] the column key to obtain the priority in the 
        matrix.
        3. Activator: [bool] a flag to indicate if the pathway is an activator 
        of the consequent (True) or not (False).
        4. Domain: [set] strings that represent the minterms of the expression
        present in the left side of the pathway. The right side is always
        exactly the letter shown in the consequent field, there is no need to 
        represent that function.
        """
        # Helper functions
        def pathway_serializer(antecedent, consequent, definition):
            """
            DESCRIPTION: 
            Helper function for the list comprehension of below. 
            :param antecedent: [str] value to write in the left side of the 
            pathway.
            :param consequent: [str] value to write in the right side of the 
            pathway.
            :param definition: [tuple] canalising/canalised value pair of the 
            pathway.
            :return: [dict] formatted pathway.
            """
            return {
                'id': ''.join(sample(ascii_uppercase + digits, 5)),
                'antecedent': antecedent,
                'consequent': consequent,
                'activator': bool(definition[1]),
                'domain': frozenset(filter(
                    lambda term: term[self.nodes.index(antecedent)] == str(definition[0]),
                    self.graph_space))}
        
        def pathway_manager(pathway_group):
            """
            DESCRIPTION:
            Helper function to organise the pathways first by node and after
            by activators and inhibitors. In the code, every activator will have
            a canalised value of 1 and vice versa.
            :param pathway_group: [tuple] the two lists of activators and 
            inhibitors.
            :return: [dict] the pathways grouped by node and activators/inhibitors.
            """
            pathways = [it for sb in pathway_group for it in sb]
            pathways = {node: {'activators': list(filter(lambda pathway: (pathway['consequent'] == node) and pathway['activator'], pathways)),
                                'inhibitors': list(filter(lambda pathway: (pathway['consequent'] == node) and not pathway['activator'], pathways))}
                for node in self.nodes}
            return pathways
        
        # Create all the pathways with both canalising/canalised pairs
        activator_pathways = [[None for activator in self.activators[node]] for node in self.nodes]
        inhibitor_pathways = [[None for inhibitor in self.inhibitors[node]] for node in self.nodes]
        i = 0
        for node in self.nodes:
            j = 0
            for activator in self.activators[node]:
                activator_pathways[i][j] = [
                    pathway_serializer(activator, node, (0, 0)),
                    pathway_serializer(activator, node, (1, 1))
                ]
                j += 1
            j = 0
            for inhibitor in self.inhibitors[node]:
                inhibitor_pathways[i][j] = [
                    pathway_serializer(inhibitor, node, (0, 1)),
                    pathway_serializer(inhibitor, node, (1, 0))
                ]
                j += 1
            i += 1
        activator_pathways = [it for sb in activator_pathways for it in sb]
        activator_pathways = itertools.product(*activator_pathways)
        inhibitor_pathways = [it for sb in inhibitor_pathways for it in sb]
        inhibitor_pathways = itertools.product(*inhibitor_pathways)
        pathways = itertools.product(activator_pathways, inhibitor_pathways)
        # Organise every group by node first and by activator/inhibitor second
        self.pathway_groups = [pathway_manager(group) for group in pathways]

    def generate_priority_matrices(self):
        """
        DESCRIPTION:
        A method that adds to the Graph object the piority matrices used in the 
        inference. There two priority matrices per simulation, one for activators
        and another for inhibitors.
        """
        # Obtain all the possible node combinations
        combinations = [itertools.combinations(self.nodes, i + 1)
            for i in range(self.n_nodes)]
        combinations = sorted([''.join(it) for sb in combinations for it in sb])
        # Generate all the possible priorities
        priorities = range(0, 1000)
        # Activators
        activator_matrices = []
        for _ in range(self.n_simulations):
            # Generate the random matrices
            activator_matrices.append({
                key: dict(zip(combinations, sample(priorities, len(combinations)))) 
                for key in combinations
                })
        # Inhibitors
        inhibitor_matrices = []
        for _ in range(self.n_simulations):
            # Generate the random matrices
            inhibitor_matrices.append({
                key: dict(zip(combinations, sample(priorities, len(combinations)))) 
                for key in combinations
                })
        # Store the matrices
        self.priority_matrices = zip(activator_matrices, inhibitor_matrices)
    
    def generate_NCBFs(self):
        """
        DESCRIPTION:
        A method to build all the groups of NCBF based on the pathway groups. 
        Multiple NCBF groups are built per pathway group. The information about the 
        canalising/canalised pairs is already stored in the pathways. The NCBF
        groups are added to the Graph object.
        """
        # Helper function
        def ncbf_formatter(ncbf_group, pathway_group_id):
            networks = [dict(zip(self.nodes, network)) for network in itertools.product(*ncbf_group)]
            return zip([pathway_group_id] * len(networks), networks)

        # Generate by node all the NCBFs
        total_ncbf = [[ncbf_generator(group[node]['activators'], group[node]['inhibitors'], self.graph_space) 
            for node in self.nodes] for group in self.pathway_groups]
        # Format all the NCBF groups conveniently: (pathway group position, NCBF
        # network in dict). 
        self.pre_networks = [it for sb in [ncbf_formatter(total_ncbf[i], i) for i in range(len(total_ncbf))] for it in sb]

    def generate_boolean_networks(self):
        """
        DESCRIPTION:
        A method to build all the boolean networks (BN), and store them in a 
        list within the Graph object. To process them independently, the BN are
        dicts formated with the following fields: 
        - priority_matrix: [dict] the random matrices that represent the priority
        between expressions.
        - pre_pathways: [tuple] the initial group of pathways. It is stored to 
        check the result at the end.
        - pre_network: [dict] frozensets with strings that represent the minterms
        of the node expressions. Every key is a different node.
        - network: [dict] the mutable version with sets of the pre_network. The
        network is where the inference takes place.
        - conflicts: [list] the conflicts, in order, handled during the inference.
        The conflicts are introduced in the list in order of appearance.
        - algorithm: [str] code to select the algorithm to use during the
        inference, new or original.   
        """
        # Helper functions
        def boolean_network_serializer(group):
            """
            DESCRIPTION:
            A function designed specifically to transle the output of the product
            from below into a dict with significant fields. The group is a tuple 
            of tuples. The first tuple is the priority matrices, the second are 
            the group of pathways to take and the network.
            """
            return {
                'priority': {'activators': group[0][0], 'inhibitors': group[0][1]},
                'pre_pathways': pathways[group[1][0]],
                'pre_network': group[1][1],
                'max_iterations': int(self.max_iterations),
                'graph_space': frozenset(self.graph_space),
                'nodes': tuple(self.nodes),
                # An empty field to store the conflicts during the inference
                'conflicts': [],
                'algorithm': self.algorithm
            }

        # Create all the groups of simulations, NCBFs and pathways
        pathways = self.pathway_groups
        boolean_networks = list(itertools.product(self.priority_matrices, self.pre_networks))
        # Format the networks and store
        self.boolean_networks = list(map(boolean_network_serializer, boolean_networks))
        print(f'Total networks generated: {len(list(self.boolean_networks))}')
    
    def solve_conflicts(self):
        """
        DESCRIPTION:
        A method to solve all the conflicts in the obtained boolean functions. No
        object is going to be created, just the previous networks are going to be 
        modified. Therefore, we create a new function that will be applied through
        a list comprehension to every boolean function.
        """
        # Measure the time taken in the inference
        # Solve the conflicts of every network
        if self.multiprocess:
            # Multiprocess inference
            with Pool(self.used_cores) as pool:
                # By default, the algorithm to iterate through the nodes is the original
                self.boolean_networks = pool.map(conflicts_solver, self.boolean_networks)
        else:
            # Single-process inference
            # By default, the algorithm to iterate through the nodes is the original
            self.boolean_networks = list(map(conflicts_solver, self.boolean_networks))

    def format_network(self):
        """
        DESCRIPTION:
        A method to compute all the information from the resulting networks. The 
        the networks are modified although they are referenced with the attribute converging
        networks.
        """
        # Measure the time taken to compute the information
        # Delete Nones
        self.boolean_networks = list(filter(lambda net: net is not None, self.boolean_networks))
        # Compute the network information with PyBoolNet
        print('Computing network information')
        if self.multiprocess:
            # Multiprocess
            with Pool(self.used_cores) as pool:
                self.boolean_networks = pool.map(network_formatter, self.boolean_networks)
        else:
            # Single-process
            self.boolean_networks = list(map(network_formatter, self.boolean_networks))

    def prefilter(self):
        """
        DESCRIPTION:
        A method to select only the networks among which the searched attractors are present.
        A cheap prefiltering before computing the attractors with the Tarjan algorithm and 
        PyBoolNet.
        """
        self.boolean_networks = list(filter(lambda network: network is not None, self.boolean_networks))
        self.boolean_networks = list(prefilter_by_attractor(self.boolean_networks, self.attractors))
        print(f'Total networks after prefiltering: {len(self.boolean_networks)}')
