#!/bin/sh
# This file is named mclaughlin.sh
#SBATCH --partition=pre
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=4

##---------------------------------------------------------------------------##
## ---------------------------- FACEMC test runner --------------------------##
##---------------------------------------------------------------------------##
## FRENSIE benchmark test: McLaughlin dose depth data.
## The dose depth for a 1-D in several materials is calculated by dividing the
## energy deposition by the subzone width (g/cm^2).
##---------------------------------------------------------------------------##
EXTRA_ARGS=$@

# Set the number of histories
HISTORIES=1000000
# Set the max runtime (in minutes, 1 day = 1440 )
TIME=1300

# These parameters can be set if the cluster is not used
# SLURM_CPUS_PER_TASK=1
# SLURM_NTASKS=1

# Run from the rendezvous
if [ "$#" -eq 1 ]; then
  # Set the rendezvous
  RENDEZVOUS="$1"

  # Restart the simulation
  echo "Restarting Facemc McLaughlin test for ${HISTORIES} particles with ${SLURM_NTASKS} MPI processes with ${SLURM_CPUS_PER_TASK} OpenMP threads each!"
  mpiexec -n ${SLURM_NTASKS} python -c "import mclaughlin; mclaughlin.runSimulationFromRendezvous(${SLURM_CPUS_PER_TASK}, ${HISTORIES}, ${TIME}, \"${RENDEZVOUS}\" )"

  directory="$(dirname "${RENDEZVOUS}")/"

# Run new simulation
else

  # Set the material ( Al, polystyrene, polyethylene )
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


  ##--------------------------------------------------------------------------##
  ## ------------------------------- COMMANDS --------------------------------##
  ##--------------------------------------------------------------------------##

  ENERGY=0.0
  SUBZONE_WIDTH=0.0
  DENSITY=0.0
  if [ "${MATERIAL}" = "Al" ]; then
    # Set the source energy
    ENERGY=3.0
    # Set the subzone width (cm)
    SUBZONE_WIDTH=0.0148
    # Set the material density (g/cm3)
    DENSITY=2.7

  elif [ "${MATERIAL}" = "polystyrene" ]; then
    # Set the source energy
    ENERGY=0.1
    # Set the subzone width (cm)
    SUBZONE_WIDTH=0.0004
    # Set the material density (g/cm3)
    DENSITY=1.06

  elif [ "${MATERIAL}" = "polyethylene" ]; then
    # Set the source energy
    ENERGY=2.0
    # Set the subzone width (cm)
    SUBZONE_WIDTH=0.022
    # Set the material density (g/cm3)
    DENSITY=0.94

  else
    echo "Material ${MATERIAL} is currently not supported!"
  fi

  # Create a unique python script and change the parameters
  python_script="mclaughlin_${SLURM_JOB_ID}"
  cp mclaughlin.py ${python_script}.py

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

  # Set the energy
  command=s/energy=.*/energy=${ENERGY}/
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

  # Get the simulation name
  name=$(python -c "import ${python_script}; name = ${python_script}.getSimulationName(); print name" 2>&1)

  RENDEZVOUS="${name}_rendezvous_0.xml"

  # Run the simulation from the last rendezvous
  if [ -f ${RENDEZVOUS} ]; then

    # Get the last rendezvous
    i=0
    while [ -f "${name}_rendezvous_${i}.xml" ]; do
      RENDEZVOUS="${name}_rendezvous_${i}.xml"
      i=$[$i+1]
    done

    # Run the simulation from the last rendezvous
    echo "Running Facemc McLaughlin test with ${HISTORIES} particles with ${SLURM_NTASKS} MPI processes with ${SLURM_CPUS_PER_TASK} OpenMP threads each from the rendezvous '${RENDEZVOUS}'!"
    mpiexec -n ${SLURM_NTASKS} python -c "import ${python_script}; ${python_script}.runSimulationFromRendezvous(${SLURM_CPUS_PER_TASK}, ${HISTORIES}, ${TIME}, '${RENDEZVOUS}' )"
  else
    # Run the simulation from the start
    echo "Running Facemc McLaughlin test with ${HISTORIES} particles with ${SLURM_NTASKS} MPI processes with ${SLURM_CPUS_PER_TASK} OpenMP threads each!"
    mpiexec -n ${SLURM_NTASKS} python -c "import ${python_script}; ${python_script}.runSimulation(${SLURM_CPUS_PER_TASK}, ${HISTORIES}, ${TIME})"
  fi

  # Remove the temperary python script
  rm ${python_script}.py*
fi

  # Move the slurm file to the results directory
  mv slurm-${SLURM_JOB_ID}.out ./${directory}