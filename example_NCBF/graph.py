#!/home/mario/Projects/boolean_2/software/venv/bin/venv python

# Libraries
import itertools
from random import sample
from multiprocessing import Pool, cpu_count
from string import ascii_uppercase, digits
from ncbf_utils import ncbf_generator
from bn_utils import conflicts_solver, minterms2bnet, network_formatter, prefilter_by_attractor


# Classes
class Graph:

    # Methods
    def __init__(self, nodes, activators, inhibitors, n_simulations, multiprocess, n_free_cores, max_iterations, algorithm, attractors, partial, **kwargs):
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
        :param partial: [bool] an indication to relax the conditions imposed in 
        the prefiltering.
        :param mixed_pathways: [bool] a flag to indicate whether the pathways in 
        boolean network and pathways should be mixed or not. If not, both are 
        computed with the same pathways.
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
        self.attractors = None if attractors == "None" else attractors
        self.partial = partial
        self.mixed_pathways = kwargs.get("mixed_pathways", False)
        self.pathways_index = kwargs.get("pathways_index", None)
        # Generate all the possible minterms in a space of len(nodes) variables.
        # IMPORTANT: the node position in every term is alphabetical: A:0, B:1...
        self.graph_space = frozenset('{:0{}b}'.format(i, self.n_nodes) 
            for i in range(2 ** self.n_nodes))
        # Check for input nodes
        self.input_nodes = tuple([node for node in self.nodes 
            if not self.activators[node] and not self.inhibitors[node]])
    
    def obtain_pathways_from_graph(self):
        """
        DESCRIPTION:
        The method to obtain all the possible groups of pathways from the graph 
        based on the pairs of canalising/canalised values described in the 
        article. The groups are added to the Graph object. Every pathway is a 
        tuple in which the fields have the following meaning:
        1. Antecedent: [str] expression describing the cause of the effect. 
        This field is the row value in the priority matrix used later to solve 
        the conflicts between pathways.
        2. Consequent: [str] expression describing the nodes that suffers the 
        effect. This field is the column value in the priority matrix used
        later to solve the conflicts between pathways.
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
        
        def pathway_manager(pathway_group, input_pathways):
            """
            DESCRIPTION:
            Helper function to organise the pathways first by node and after
            by activators and inhibitors. In the code, every activator will have
            a canalised value of 1 and vice versa.
            :param pathway_group: [tuple] the two lists of activators and 
            inhibitors.
            :param input_pathways: [dict] pathways of the input nodes. They are 
            the same for every group.
            :return: [dict] the pathways grouped by node and activators/inhibitors.
            """
            # Format the pathways
            pathways = [it for sb in pathway_group for it in sb]
            pathways = {node: {'activators': list(filter(lambda pathway: (pathway['consequent'] == node) and pathway['activator'], pathways)),
                                'inhibitors': list(filter(lambda pathway: (pathway['consequent'] == node) and not pathway['activator'], pathways))}
                for node in self.nodes}
            # Introduce the input pathways
            [pathways.update({node: input_pathways[node]}) for node in input_pathways.keys()]
            return pathways
        
        # Create all the pathways with both canalising/canalised pairs
        activator_pathways = [[None] * len(self.activators[node]) for node in self.nodes]
        inhibitor_pathways = [[None] * len(self.inhibitors[node]) for node in self.nodes]
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
        # Obtain the pathways of the inputs (from and to themselves)
        input_pathways = {node: [] for node in self.input_nodes}
        for node in self.input_nodes:
            input_pathways[node] = {
                'activators': [pathway_serializer(node, node, (1, 1))],
                'inhibitors': [pathway_serializer(node, node, (0, 0))]
            }
        # Organise every group by node first and by activator/inhibitor second
        self.pathway_groups = [pathway_manager(group, input_pathways) for group in pathways]

    def generate_priority_matrices(self):
        """
        DESCRIPTION:
        A method that adds to the Graph object the priority matrices used in the 
        inference. There are two priority matrices per simulation, one for activators
        and another for inhibitors.
        """
        # Obtain all the possible node combinations
        combinations = [itertools.combinations(self.nodes, i + 1)
            for i in range(self.n_nodes)]
        combinations = sorted([''.join(it) for sb in combinations for it in sb])
        # Generate all the possible priorities
        # We set arbitrarily all the priorities between 0 and 1000, where 0 is the
        # lowest priority and 1000 the highest
        priorities = range(0, 1000)
        # Activators
        activator_matrices = []
        # We generate as many priority matrices as simulations that we will perform
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
        # Helper functions
        # There two depending on whether the pathways groups are mixed or not for
        # the ncbf and the conflicts strategy
        def ncbf_formatter_standard(ncbf_group, pathway_group_id):
            # This function assigns the ID of the pathways used to build the network
            networks = [dict(zip(self.nodes, network)) for network in itertools.product(*ncbf_group)]
            return zip([pathway_group_id] * len(networks), networks)

        def ncbf_formatter_mixed_pathways(ncbf_group):
            # This assigns the prefixed ID of pathways independently of the network
            networks = [dict(zip(self.nodes, network)) for network in itertools.product(*ncbf_group)]
            return zip([self.pathways_index] * len(networks), networks)

        # Generate by node all the NCBFs
        total_ncbf = [[ncbf_generator(group[node]['activators'], group[node]['inhibitors'], self.graph_space) 
            for node in self.nodes] for group in self.pathway_groups]
        # Format all the NCBF groups conveniently: (pathway group position, NCBF
        # network in dict). 
        if self.mixed_pathways:
            pre_networks = [it for sb in [ncbf_formatter_mixed_pathways(total_ncbf[i]) for i in range(len(total_ncbf))] for it in sb]
        else:
            pre_networks = [it for sb in [ncbf_formatter_standard(total_ncbf[i], i) for i in range(len(total_ncbf))] for it in sb]
        # Filter the equivalent networks
        codes = []
        final_pre_networks = []
        for network in pre_networks:
            # code = str(network[0]) + '$$' + '&'.join(['|'.join(sorted(net)) for _, net in sorted(network[1].items(), key=lambda x: x[0])])
            code = '&'.join(['|'.join(sorted(net)) for _, net in sorted(network[1].items(), key=lambda x: x[0])])
            if code not in codes:
                final_pre_networks.append(network)
                codes.append(code)
        # Store NCBF networks
        self.pre_networks = final_pre_networks
    
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
            # Calculate the expressions of the prenetwork (just for information)
            pre_expressions = dict(zip(
                group[1][1].keys(),
                [minterms2bnet(list(group[1][1].keys()), group[1][1][node]) for node in group[1][1].keys()]
            ))
            # Build and return the formatted network
            return {
                'priority': {'activators': group[0][0], 'inhibitors': group[0][1]},
                'pre_pathways': pathways[group[1][0]],
                'pre_network': group[1][1],
                'max_iterations': int(self.max_iterations),
                'graph_space': frozenset(self.graph_space),
                'nodes': tuple(self.nodes),
                'pre_expressions': pre_expressions,
                # An empty field to store the conflicts during the inference
                'conflicts': [],
                'algorithm': self.algorithm,
                'input_nodes': self.input_nodes
            }

        # Create all the groups of simulations, NCBFs and pathways
        pathways = self.pathway_groups
        boolean_networks = list(itertools.product(self.priority_matrices, self.pre_networks))
        # Format the networks and store
        self.boolean_networks = list(map(boolean_network_serializer, boolean_networks))
        print(f'Total networks generated: {len(self.boolean_networks)}')
    
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
        if self.boolean_networks:
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
        if (self.attractors is not None) and (self.attractors != []):
            self.boolean_networks = list(filter(lambda network: network is not None, self.boolean_networks))
            if self.boolean_networks:
                self.boolean_networks = list(prefilter_by_attractor(self.boolean_networks, self.attractors, self.partial))
            print(f'Total networks after prefiltering: {len(self.boolean_networks)}')
        else:
            print('No attractor-based prefiltering performed')
