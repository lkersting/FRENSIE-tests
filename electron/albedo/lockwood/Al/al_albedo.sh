#!/bin/sh
# This file is named albedo.sh
#SBATCH --partition=pre
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=4

##---------------------------------------------------------------------------##
## ---------------------------- FACEMC test runner --------------------------##
##---------------------------------------------------------------------------##
## FRENSIE benchmark test: Reflection Coeficient in semi-infinite slabs.
## The albedo for a particular material is calculated at several energies.
##---------------------------------------------------------------------------##
EXTRA_ARGS=$@

# Set the number of histories
HISTORIES=1e7
# Set the max runtime (in minutes, 1 day = 1440 )
TIME=1350

# These parameters can be set if the cluster is not used
# SLURM_CPUS_PER_TASK=4
# SLURM_NTASKS=1

# Run from the rendezvous
if [ "$#" -eq 1 ]; then
  # Set the rendezvous
  RENDEZVOUS="$1"

  # Restart the simulation
  echo "Restarting Facemc Albedo test for ${HISTORIES} particles with ${SLURM_NTASKS} MPI processes with ${SLURM_CPUS_PER_TASK} OpenMP threads each!"

  RENDEZVOUS="${PWD}/${RENDEZVOUS}"
  python_path="$(dirname "${PWD}")/"

  mpiexec -n ${SLURM_NTASKS} python2.7 -c "import sys; sys.path.insert(1,'${python_path}'); import albedo_simulation; albedo_simulation.runSimulationFromRendezvous(${SLURM_CPUS_PER_TASK}, ${HISTORIES}, ${TIME}, '${RENDEZVOUS}' )"

  directory="$(dirname "${RENDEZVOUS}")/"

# Run new simulation
else

  # Set the bivariate Grid Policy ( UNIT_BASE_CORRELATED CORRELATED UNIT_BASE )
  GRID_POLICY=UNIT_BASE

  # Set the elastic distribution mode ( DECOUPLED COUPLED HYBRID )
  MODE=COUPLED

  # Set the elastic coupled sampling method ( ONE_D TWO_D MODIFIED_TWO_D )
  METHOD=TWO_D

  # Set the electron cutoff energy
  CUTOFF=1e-4

  # Set the transport mode ( "forward", "adjoint" )
  TRANSPORT="forward"

  ## ------- FORWARD OPTIONS ------- ##
  # Set the bivariate interpolation ( LOGLOGLOG LINLINLIN LINLINLOG )
  INTERP=LOGLOGLOG

  # Set the test source energy
  ENERGY=1.033

  # Set the test source angle in degrees (0, 60)
  ANGLE=0.0

  # Set the data file type (ACE Native)
  FILE_TYPE=Native

  # Set if a refined grid should be used ( True, False )
  REFINED=True

  # Set if a spectrum source should be used ( True, False )
  SPECTRUM=True

  # Set if a isotrpoic source should be used ( True, False )
  ISOTROPIC=True

  ##---------------------------------------------------------------------------##
  ## ------------------------------- COMMANDS ---------------------------------##
  ##---------------------------------------------------------------------------##

  # Create a unique python script and change the parameters
  script_name="al_albedo_${SLURM_JOB_ID}"
  python_script="${script_name}.py"
  cp al_albedo.py ${python_script}

  # Change the python_script parameters

  # Set 2D grid policy
  command=s/grid_policy=MonteCarlo.*/grid_policy=MonteCarlo.${GRID_POLICY}_GRID/
  sed -i "${command}" ${python_script}

  # Set the elastic distribution mode
  command=s/mode=MonteCarlo.*/mode=MonteCarlo.${MODE}_DISTRIBUTION/
  sed -i "${command}" ${python_script}

  # Set the elastic coupled sampling method
  command=s/method=MonteCarlo.*/method=MonteCarlo.${METHOD}_UNION/
  sed -i "${command}" ${python_script}

  # Set the cutoff energy
  command=s/cutoff_energy=.*/cutoff_energy=${CUTOFF}/
  sed -i "${command}" ${python_script}

  if [ "${TRANSPORT}" = "forward" ]; then

    # Set if a spectrum source should be used
    command=s/spectrum_source=.*/spectrum_source=${SPECTRUM}/
    sed -i "${command}" ${python_script}

    # Set if a isotropic source should be used
    command=s/isotropic_source=.*/isotropic_source=${ISOTROPIC}/
    sed -i "${command}" ${python_script}

    # Set the source energy
    command=s/source_energy=.*/source_energy=${ENERGY}/
    sed -i "${command}" ${python_script}

    # Set the source angle
    command=s/source_angle=.*/source_angle=${ANGLE}/
    sed -i "${command}" ${python_script}

    # Set the interp
    command=s/interpolation=MonteCarlo.*/interpolation=MonteCarlo.${INTERP}_INTERPOLATION/
    sed -i "${command}" ${python_script}

    # Set the file type
    command=s/file_type=Data.ElectroatomicDataProperties.*/file_type=Data.ElectroatomicDataProperties.${FILE_TYPE}_EPR_FILE/
    sed -i "${command}" ${python_script}

    # Set if a refined grid should be used
    command=s/use_refined_grid=.*/use_refined_grid=${REFINED}/
    sed -i "${command}" ${python_script}

    # Get the simulation name
    name=$(python -c "import ${script_name}; ${script_name}.printSimulationName();" 2>&1)

  elif [ "${TRANSPORT}" = "adjoint" ]; then

    # Get the simulation name
    name=$(python -c "import ${script_name}; ${script_name}.printAdjointSimulationName();" 2>&1)

  fi

  RENDEZVOUS="${PWD}/${name}_rendezvous_0.xml"

  directory="$(dirname "${RENDEZVOUS}")/"

  # Run the simulation from the last rendezvous
  if [ -f ${RENDEZVOUS} ]; then

    # Get the last rendezvous
    i=0
    while [ -f "${PWD}/${name}_rendezvous_${i}.xml" ]; do
      RENDEZVOUS="${PWD}/${name}_rendezvous_${i}.xml"
      i=$[$i+1]
    done

    # Run the simulation from the last rendezvous
    echo "Running Facemc ${TRANSPORT} albedo test with ${HISTORIES} particles with ${SLURM_NTASKS} MPI processes with ${SLURM_CPUS_PER_TASK} OpenMP threads each from the rendezvous '${RENDEZVOUS}'!"
    mpiexec -n ${SLURM_NTASKS} python -c "from os import path; import sys; sys.path.insert(1, '../'); import albedo_simulation; albedo_simulation.runSimulationFromRendezvous(${SLURM_CPUS_PER_TASK}, ${HISTORIES}, ${TIME}, '${RENDEZVOUS}' )"

  else
    # Run the simulation from the start
    echo "Running Facemc ${TRANSPORT} albedo test with ${HISTORIES} particles with ${SLURM_NTASKS} MPI processes with ${SLURM_CPUS_PER_TASK} OpenMP threads each!"
    mpiexec -n ${SLURM_NTASKS} python ${python_script} --num_particles=${HISTORIES} --threads=$SLURM_CPUS_PER_TASK --time=${TIME} --transport=${TRANSPORT}
  fi

  # Remove the temperary python script
  rm ${python_script}*
fi

  # Move the slurm file to the results directory
  mv slurm-${SLURM_JOB_ID}.out ${directory}
