#!/bin/sh
# This file is named broomstick.sh
#SBATCH --partition=pre
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=4

mpirun -np $SLURM_NTASKS python2.7 dyson_sphere.py --db_path=$DATABASE_PATH --sim_name="dyson_sphere" --num_particles=1e9 --threads=$SLURM_CPUS_PER_TASK
