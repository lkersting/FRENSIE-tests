#!/bin/bash
##---------------------------------------------------------------------------##
## ---------------------------- FACEMC test runner --------------------------##
##---------------------------------------------------------------------------##
## The electron albedo is found for a semi-infinite aluminum slab.
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
HISTORIES=1000000
# Set the max runtime (in minutes, 1 day = 1440 )
TIME=1350

# Set the bivariate Grid Policy ( UNIT_BASE_CORRELATED CORRELATED UNIT_BASE )
grid_policys=( UNIT_BASE_CORRELATED )

# Set the elastic distribution mode ( DECOUPLED COUPLED HYBRID )
modes=( COUPLED )

# Set the elastic coupled sampling method
# ( ONE_D TWO_D MODIFIED_TWO_D )
methods=( MODIFIED_TWO_D )

# Set the electron cutoff energy ( 1e-4 )
cutoffs=( 1e-4 )

# Set the transport mode ( "adjoint" "forward" )
transports=( "forward" )

## ------- FORWARD OPTIONS ------- ##

# Set the test energy (0.0002 0.0003 0.0004 0.0005 0.0006 0.0008 0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045 0.005 0.006 0.0093 0.01 0.011 0.0134 0.015 0.0173 0.02 0.0252 0.03 0.04 0.0415 0.05 0.06 0.0621 0.07 0.08 0.0818 0.1 0.102 0.121 0.146 0.172 0.196 0.2 0.238 0.256 )
# energies can be set to indivual energies e.g. ( 0.02 0.03 0.04 ) or "all"
energies="all"

# Set the data file type (ACE Native)
file_types=( Native )

# Set if a refined grid should be used ( "True" "False" )
refined_grids=( "False" )

# Set the bivariate interpolation ( LOGLOGLOG LINLINLIN LINLINLOG )
interps=( LOGLOGLOG )

## ------- ADJOINT OPTIONS ------- ##
# Set the nudge past max energy mode ( False True )
nudges=( True )

# Set the electro-ionization sampling mode ( OUTGOING_ENERGY KNOCK_ON )
ionizations=( KNOCK_ON )

##---------------------------------------------------------------------------##
## ------------------------------- COMMANDS ---------------------------------##
##---------------------------------------------------------------------------##

# Set the energies to all
if [ ${energies} == "all" ]; then
    energies=(0.0002 0.0003 0.0004 0.0005 0.0006 0.0008 0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045 0.005 0.006 0.0093 0.01 0.011 0.0134 0.015 0.0173 0.02 0.0252 0.03 0.04 0.0415 0.05 0.06 0.0621 0.07 0.08 0.0818 0.1 0.102 0.121 0.146 0.172 0.196 0.2 0.238 0.256 )
fi

# Material element
element="Al"; zaid=13000

# Move to element directory
cd ${element}

# Set the script name
script=al_albedo.sh

# Set the number of threads
command="s/\#SBATCH[[:space:]]--ntasks=.*/\#SBATCH --ntasks=${MPI_PROCESSES}/"
sed -i "${command}" ${script}
command="s/\#SBATCH[[:space:]]--cpus-per-task=.*/\#SBATCH --cpus-per-task=${OPEN_MP_THREADS}/"
sed -i "${command}" ${script}

command=s/TIME=.*/TIME=${TIME}/
sed -i "${command}" ${script}
command=s/HISTORIES=.*/HISTORIES=${HISTORIES}/
sed -i "${command}" ${script}

for transport in "${transports[@]}"
do
  # Set the transport mode
  command=s/TRANSPORT=.*/TRANSPORT=\"${transport}\"/
  sed -i "${command}" ${script}
  echo "Setting the transport mode to ${transport}"

  if [ "${transport}" = "forward" ]; then
    for file_type in "${file_types[@]}"
    do
      # Set the file type
      command=s/FILE_TYPE=.*/FILE_TYPE=${file_type}/
      sed -i "${command}" ${script}
      echo "Setting file type to ${file_type}"

      if [ "${file_type}" = "Native" ]; then

        for interp in "${interps[@]}"
        do
          # Set the interp
          command=s/INTERP=.*/INTERP=${interp}/
          sed -i "${command}" ${script}
          echo "  Setting interpolation to ${interp}"

          for grid_policy in "${grid_policys[@]}"
          do
            if [ "${interp}" == "LINLINLOG" ] && [ "${grid_policy}" == "CORRELATED" ]; then
              echo "    The interp (${interp}) and grid policy (${grid_policy}) combo will be skipped."
            else
              # Set 2D grid policy
              command=s/GRID_POLICY=.*/GRID_POLICY=${grid_policy}/
              sed -i "${command}" ${script}
              echo "    Setting grid policy to ${grid_policy}"

              # Set the refined grid mode on/off
              for refined_grid in "${refined_grids[@]}"
              do
                # Set if a refined grid should be used
                command=s/REFINED=.*/REFINED=${refined_grid}/
                sed -i "${command}" ${script}

                for mode in "${modes[@]}"
                do
                  # Set the elastic distribution mode
                  command=s/MODE=.*/MODE=${mode}/
                  sed -i "${command}" ${script}
                  echo "      Setting elastic mode to ${mode}"

                  if [ "${mode}" == "COUPLED" ]; then

                    for method in "${methods[@]}"
                    do
                      # Set the elastic coupled sampling method
                      command=s/METHOD=.*/METHOD=${method}/
                      sed -i "${command}" ${script}
                      echo "        Setting elastic coupled sampling method to ${method}"

                        # loop through test energies and run mpi script
                        for energy in "${energies[@]}"
                        do
                            # Set the energy
                            command=s/ENERGY=.*/ENERGY=${energy}/
                            sed -i "${command}" ${script}

                          echo -e "          Running Albedo at ${energy} MeV!\n"
                          sbatch ${script}
                        done
                    done
                  else
                    # loop through test numbers and run mpi script
                    for energy in "${energies[@]}"
                    do
                        # Set the energy
                        command=s/ENERGY=.*/ENERGY=${energy}/
                        sed -i "${command}" ${script}

                        echo -e "        Running Albedo at ${energy} MeV!\n"
                        sbatch ${script}
                    done
                  fi
                done
              done
            fi
          done
        done
      else
        # loop through test numbers and run mpi script
        for energy in "${energies[@]}"
        do
            # Set the energy
            command=s/ENERGY=.*/ENERGY=${energy}/
            sed -i "${command}" ${script}

            echo -e "  Running Albedo at ${energy} MeV!\n"
            sbatch ${script}
        done
      fi
    done
  elif [ "${transport}" = "adjoint" ]; then
    for grid_policy in "${grid_policys[@]}"
    do
      # Set 2D grid policy
      command=s/GRID_POLICY=.*/GRID_POLICY=${grid_policy}/
      sed -i "${command}" ${script}
      echo "  Setting grid policy to ${grid_policy}"

      # Set the nudge past max energy mode on/off
      for nudge in "${nudges[@]}"
      do
        # Set the nudge past max energy mode on/off
        command=s/NUDGE=.*/NUDGE=${nudge}/
        sed -i "${command}" ${script}
        echo "    Setting the nudge past max energy mode to ${nudge}"

        for ionization in "${ionizations[@]}"
        do
          # Set the electro-ionization sampling mode
          command=s/IONIZATION=.*/IONIZATION=${ionization}/
          sed -i "${command}" ${script}
          echo "      Setting electro-ionization sampling to ${ionization}"

          for mode in "${modes[@]}"
          do
            # Set the elastic distribution mode
            command=s/MODE=.*/MODE=${mode}/
            sed -i "${command}" ${script}
            echo "        Setting elastic mode to ${mode}"

            if [ "${mode}" == "COUPLED" ]; then

              for method in "${methods[@]}"
              do
                # Set the elastic coupled sampling method
                command=s/METHOD=.*/METHOD=${method}/
                sed -i "${command}" ${script}
                echo "            Setting elastic coupled sampling method to ${method}"

                sbatch ${script}
              done
            else
              sbatch ${script}
            fi
          done
        done
      done
    done
  fi
done
