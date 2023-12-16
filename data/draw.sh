#!/bin/bash

# Request a GPU partition node and access to 1 GPU


# Ensures all allocated cores are on the same node
#SBATCH -N 1
#SBATCH --mem=64G
# Request 1 CPU core
#SBATCH -n 4
#SBATCH -t 10:30:00
#SBATCH -o with_pyhton.out
#SBATCH -e with_python.err


python3 handeling.py
