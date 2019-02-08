#!/bin/bash
# This file is named infinite_medium.sh

# Sbatch variables
partition=pre
time=1-00:00:00
ntasks=40
threads=4

sbatch_command="sbatch --partition=${partition} --time=${time} --ntasks=${ntasks} --cpus-per-task=${threads}"
${sbatch_command} mpirun -np ${ntasks} python2.7 infinite_medium.py --num_particles=1e8 --threads=${threads}
