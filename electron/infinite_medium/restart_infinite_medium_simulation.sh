#!/bin/bash
##---------------------------------------------------------------------------##
## -------------------------- FRENSIE test runner ---------------------------##
##---------------------------------------------------------------------------##
## Restart an infinite medium simulation
##---------------------------------------------------------------------------##

# Sbatch variables
partition=pre
time=1-00:00:00
ntasks=40
threads=4

##---------------------------------------------------------------------------##
## ------------------------- SIMULATION VARIABLES ---------------------------##
##---------------------------------------------------------------------------##

# Desired number of histories
num_particles=1e7
# Desired wall tim
wall_time=1350

# Set the rendezvous
RENDEZVOUS=$(realpath "$1")
# Set the file's current directory
current_directory=$PWD
# Set the original run directory
run_directory="$(dirname $(dirname $(dirname $(dirname "${RENDEZVOUS}"))))/"


sbatch_command="sbatch --partition=${partition} --time=${time} --ntasks=${ntasks} --cpus-per-task=${threads}"
run_date=$(date +'%Y-%m-%d')
mv_slurm_command="\nmv slurm-${SLURM_JOB_ID}.out ./results/${run_date}"

if ! type sbatch > /dev/null 2>&1; then
  sbatch_command=bash
  ntasks=1
fi

# Restart the simulation
echo "Restarting infinite medium test for ${num_particles} particles from rendezvous ${RENDEZVOUS}!"
python_command="mpirun -np ${ntasks} python2.7 -c \"import sys; sys.path.insert(1,\'${current_directory}\'); import infinite_medium_simulation; infinite_medium_simulation.restartInfiniteMediumSimulation( \'${RENDEZVOUS}\', ${num_particles}, ${threads}, ${wall_time} )\""

cd ${run_directory}

printf "#!/bin/bash\n${python_command}${mv_slurm_command}" > infinite_medium_temp.sh

${sbatch_command} infinite_medium_temp.sh
rm infinite_medium_temp.sh
