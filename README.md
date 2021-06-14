# Algorithm for the inference of boolean networks

### Installation instructions
Install libraries
```
pip3 install requirements.txt
```
**Note**: PyBoolNet (version 2.31.0) should be installed manually, the instructions can be found [here](https://github.com/hklarner/PyBoolNet).

### How to use
The information to represent both the problem and software execution should be introduced in the file data.json. The file consists of the following fields:

- *Nodes* (nodes): name of the nodes that build the graph. They can only be letters in uppercase.
- *Activators* (```activators```): activators list for every node in the graph (sharp arrows). The *key* represents the activated node. 
- *Inhibitors* (```inhibitors```): inhibitors list for every node in the graph (plain arrows). The *key* represents the inhibited node. 
- *Attractors* (```attractors```): list of searched attractors. The attractors used to filter the networks.
- *Number of attractors* (```n_attractors```): the number of attractors that should have the desired networks.
- *Number of simulations* (```n_simulations```): the number of times that every inference is repeated, every time with a different priority matrix.
- *Multiprocess* (```multiprocess```): boolean flag to indicate if the code should use several cores instead of 1.
- *Number of free cores* (```n_free_cores```): when multiprocess=True, the code uses all the available cores but the number specified in this field.
- *Maximum number of iterations*  (```max_iterations```): the number of iterations over all the nodes before declaring the network as *not convergent*.
- *Algorithm* (```algorithm```): the algorithm use to manage the conflicts within the node. There are two alternatives: original and updated.
- *Load previous networks* (```load_session```): boolean flag to choose between infering and sotring new networks and loading previous networks.
- *Delete repeated* (```deleate_repeated```): when load_session=true, this flag deletes the repeated networks among the loaded networks. The resulting networks are stored in an independent file in the networks folder.
- *Partial* (```partial```): boolean flag to relax the filtering conditions. When partial=true, networks with at most n_attractors and at least one of the specified attractors can pass the filtering. 

### Files and folders

The code stores the obtained networks in the folder *networks* . Every execution creates a pickle file with all the networks in that execution. Besides, every time the code is launched, one example network is exported in a JSON file (example.json). 