#!/bin/bash
# This file is named dyson_sphere_mcnp.sh
#SBATCH --partition=pre
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

/home/lkersting/software/mcnp6.2/bin/mcnp6 i=dyson_sphere_mcnp.i o=dyson_sphere_mcnp.o tasks $SLURM_CPUS_PER_TASK
