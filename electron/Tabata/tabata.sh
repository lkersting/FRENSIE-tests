#!/bin/sh
# This file is named tabata.sh
#SBATCH --partition=pre
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=4

##---------------------------------------------------------------------------##
## ---------------------------- FACEMC test runner --------------------------##
##---------------------------------------------------------------------------##
## FRENSIE benchmark test: Tabata charge deposition data.
## The charge deposition for several materials in 1-D is calculated by dividing
## the charge deposition by the subzone width (g/cm^2).
##---------------------------------------------------------------------------##
EXTRA_ARGS=$@

# Set the number of histories
HISTORIES=100000
# Set the max runtime (in minutes, 1 day = 1440 )
TIME=1300

# Run from the rendezvous
if [ "$#" -eq 1 ]; then
  # Set the rendezvous
  RENDEZVOUS="$1"

  # Restart the simulation
  echo "Restarting Facemc Tabata test for ${HISTORIES} particles with ${SLURM_NTASKS} MPI processes with ${SLURM_CPUS_PER_TASK} OpenMP threads each!"
  mpiexec -n ${SLURM_NTASKS} python -c "import tabata; tabata.restartSimulation(${SLURM_CPUS_PER_TASK}, ${HISTORIES}, ${TIME}, \"${RENDEZVOUS}\" )"

  directory="$(dirname "${RENDEZVOUS}")/"

# Run new simulation
else

  # Set the material ( Al, Be )
  MATERIAL="Al"

  # Set the data file type (ACE Native)
  FILE_TYPE=Native

  # Set the bivariate interpolation ( LOGLOGLOG LINLINLIN LINLINLOG )
  INTERP=LOGLOGLOG

  # Set the bivariate Grid Policy ( UNIT_BASE_CORRELATED CORRELATED UNIT_BASE )
  GRID_POLICY=UNIT_BASE_CORRELATED

  # Set the elastic distribution mode ( DECOUPLED COUPLED HYBRID )
  MODE=COUPLED

  # Set the elastic coupled sampling method ( ONE_D TWO_D MODIFIED_TWO_D )
  METHOD=MODIFIED_TWO_D

  ##---------------------------------------------------------------------------##
  ## ------------------------------- COMMANDS ---------------------------------##
  ##---------------------------------------------------------------------------##

  SUBZONE_WIDTH=0.0
  DENSITY=0.0
  if [ "${MATERIAL}" = "Al" ]; then
    # Set the subzone width (cm)
    SUBZONE_WIDTH=5.0
    # Set the material density (g/cm3)
    DENSITY=2.7

  elif [ "${MATERIAL}" = "Be" ]; then
    # Set the subzone width (cm)
    SUBZONE_WIDTH=10.0
    # Set the material density (g/cm3)
    DENSITY=1.85

  else
    echo "Material ${MATERIAL} is currently not supported!"
  fi


  # Create a unique python script and change the parameters
  python_script="tabata_${SLURM_JOB_ID}"
  cp tabata.py ${python_script}.py

  # Change the python_script parameters

  # Set the material
  command=s/material=.*/material=\"${MATERIAL}\"/
  sed -i "${command}" ${python_script}.py

  # Set the material density
  command=s/density=.*/density=${DENSITY}/
  sed -i "${command}" ${python_script}.py

  # Set the subzone width
  command=s/subzone_width=.*/subzone_width=${SUBZONE_WIDTH}/
  sed -i "${command}" ${python_script}.py

  # Set the file type
  command=s/file_type=Data.ElectroatomicDataProperties.*/file_type=Data.ElectroatomicDataProperties.${FILE_TYPE}_EPR_FILE/
  sed -i "${command}" ${python_script}.py

  # Set the interp
  command=s/interpolation=MonteCarlo.*/interpolation=MonteCarlo.${INTERP}_INTERPOLATION/
  sed -i "${command}" ${python_script}.py

  # Set 2D grid policy
  command=s/grid_policy=MonteCarlo.*/grid_policy=MonteCarlo.${GRID_POLICY}_GRID/
  sed -i "${command}" ${python_script}.py

  # Set the elastic distribution mode
  command=s/mode=MonteCarlo.*/mode=MonteCarlo.${MODE}_DISTRIBUTION/
  sed -i "${command}" ${python_script}.py

  # Set the elastic coupled sampling method
  command=s/method=MonteCarlo.*/method=MonteCarlo.${METHOD}_UNION/
  sed -i "${command}" ${python_script}.py

  # Create the results directory
  directory=$(python -c "import ${python_script}; ${python_script}.createResultsDirectory()" 2>&1)

  # Run the simulation
  echo "Running Facemc Tabata test with ${HISTORIES} particles with ${SLURM_NTASKS} MPI processes with ${SLURM_CPUS_PER_TASK} OpenMP threads each!"
  mpiexec -n ${SLURM_NTASKS} python -c "import ${python_script}; ${python_script}.runSimulation(${SLURM_CPUS_PER_TASK}, ${HISTORIES}, ${TIME})"

  # Remove the temperary python script
  rm ${python_script}.py*
fi

  # Move the slurm file to the results directory
  mv slurm-${SLURM_JOB_ID}.out ./${directory}