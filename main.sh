#!/bin/bash

################################################################################
# DESCRIPTION
################################################################################
# A file to execute the script without debugging. Programing too many networks
# directly in the python script will overwhelm the memory. The shell scripts 
# repeats the process with a small amount of networks, which is equivalent.
################################################################################

# Load the virtual environment
source ./venv/bin/activate

# Adding working directory to python path
WORKING_DIRECTORY="/home/mario/Projects/boolean_2/software"
export PYTHONPATH="${PYTHONPATH}:$WORKING_DIRECTORY"

# Set the number of repeats
N_REPETITIONS=255

# Launch the script
echo "ANALYSIS STARTED"
echo "Launching instances"
for i in  $(seq 0 1 $N_REPETITIONS)
do
    echo "----------------------------------------------------------"
    echo "Instance $i"
    echo "----------------------------------------------------------"
    python3 main.py $i
done
echo "----------------------------------------------------------"
echo "ANALYSIS COMPLETED"