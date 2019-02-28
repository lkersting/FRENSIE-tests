#!/bin/bash
##---------------------------------------------------------------------------##
## ---------------------------- FACEMC test runner --------------------------##
##---------------------------------------------------------------------------##
## Validation runs comparing FRENSIE and MCNP.
## The electron angular distribution for a thin gold foil of .0009658 cm.
## The # of particles per steradian for scattering angle is found by dividing
## the surface current by 2pi * ( \mu_{i} - \mu_{i-1} ) where \mu_{0} is the
## lowest cosine bin (ie: -1). Surface current is needed so DagMC will be used.
## The #/steradians can be changed to #/square degree by multiplying by
## (pi/180)^2.
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## ---------------------------- TEST VARIABLES ------------------------------##
##---------------------------------------------------------------------------##
EXTRA_ARGS=$@

# Set the number of mpi processes and openMP threads
# NOTE: OpenMP threads should be a factor of 16 for univ and 20 for univ2
# NOTE: the max OpenMP threads should be <= 6
MPI_PROCESSES=40
OPEN_MP_THREADS=4

# Set the number of histories
HISTORIES=1e6
# Set the max runtime (in minutes, 1 day = 1440 )
TIME=1350

# Set the data file type (ACE Native)
file_types=( Native )

# Set if a refined grid should be used ( True False )
refined_grids=( True )

# Set the bivariate interpolation ( LOGLOGLOG LINLINLIN LINLINLOG )
interps=( LOGLOGLOG )

# Set the bivariate Grid Policy ( UNIT_BASE UNIT_BASE_CORRELATED CORRELATED )
grid_policys=( UNIT_BASE UNIT_BASE_CORRELATED CORRELATED )

# Set the elastic distribution mode ( DECOUPLED COUPLED HYBRID )
modes=( COUPLED )

# Set the elastic coupled sampling method
# ( ONE_D TWO_D MODIFIED_TWO_D )
methods=( TWO_D )

##---------------------------------------------------------------------------##
## ------------------------------- COMMANDS ---------------------------------##
##---------------------------------------------------------------------------##
script=hanson.sh

bold=$(tput bold)
normal=$(tput sgr0)

# Set the number of threads
command="s/\#SBATCH[[:space:]]--ntasks=.*/\#SBATCH --ntasks=${MPI_PROCESSES}/"
sed -i "${command}" ${script}
command="s/\#SBATCH[[:space:]]--cpus-per-task=.*/\#SBATCH --cpus-per-task=${OPEN_MP_THREADS}/"
sed -i "${command}" ${script}

# Set the wall time and number of histories
command=s/TIME=.*/TIME=${TIME}/
sed -i "${command}" ${script}
command=s/HISTORIES=.*/HISTORIES=${HISTORIES}/
sed -i "${command}" ${script}


for file_type in "${file_types[@]}"
do
  # Set the file type
  command=s/FILE_TYPE=.*/FILE_TYPE=${file_type}/
  sed -i "${command}" ${script}
  echo "Setting file type to ${bold}${file_type}${normal}"

  if [ "${file_type}" = "Native" ]; then

    for interp in "${interps[@]}"
    do
      # Set the interp
      command=s/INTERP=.*/INTERP=${interp}/
      sed -i "${command}" ${script}
      echo "  Setting interpolation to ${bold}${interp}${normal}"

      for grid_policy in "${grid_policys[@]}"
      do
        if [ "${interp}" == "LINLINLOG" ] && [ "${grid_policy}" == "CORRELATED" ]; then
          echo "    The interp (${bold}${interp}${normal}) and grid policy (${bold}${grid_policy}${normal}) combo will be skipped."
        else
          # Set bivariate grid policy
          command=s/GRID_POLICY=.*/GRID_POLICY=${grid_policy}/
          sed -i "${command}" ${script}
          echo "    Setting the bivariate grid policy to ${bold}${grid_policy}${normal}"

          # Set the refined grid mode on/off
          for refined_grid in "${refined_grids[@]}"
          do
            # Set if a refined grid should be used
            command=s/REFINED=.*/REFINED=${refined_grid}/
            sed -i "${command}" ${script}
            echo "      Setting refined grid mode to ${bold}${refined_grid}${normal}"

            for mode in "${modes[@]}"
            do
              # Set the elastic distribution mode
              command=s/MODE=.*/MODE=${mode}/
              sed -i "${command}" ${script}
              echo "        Setting elastic mode to ${bold}${mode}${normal}"

              if [ "${mode}" == "COUPLED" ]; then

                for method in "${methods[@]}"
                do
                  # Set the elastic coupled sampling method
                  command=s/METHOD=.*/METHOD=${method}/
                  sed -i "${command}" ${script}
                  echo "          Setting elastic coupled sampling method to ${bold}${method}${normal}"

                  sbatch ${script}

                done
              else
                  sbatch ${script}
              fi
            done
          done
        fi
      done
    done
  else
    sbatch ${script}
  fi

done
