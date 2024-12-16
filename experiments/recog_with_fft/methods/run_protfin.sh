#!/bin/bash

# A script to execute the protfin.py command.

# Define variables for the command components
PYTHON_COMMAND="python3"
SCRIPT_NAME="protfin.py"
ACTION="create-db"
INPUT_FILE="protein.fa"

# Execute the command
$PYTHON_COMMAND $SCRIPT_NAME $ACTION $INPUT_FILE
