#!/bin/bash
##---------------------------------------------------------------------------##
## ---------------------------- FACEMC test runner --------------------------##
##---------------------------------------------------------------------------##
## FRENSIE benchmark test: Mclaughlin dose depth data.
## The dose depth for a 1-D in several materials is calculated by dividing the
## energy deposition by the subzone width (g/cm^2).
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
TIME=1300

# Set the material ( Al polystyrene polyethylene )
materials=( Al )

# Set the data file type (ACE Native)
file_types=( Native )

# Set the bivariate interpolation ( LOGLOGLOG LINLINLIN LINLINLOG )
interps=( LOGLOGLOG )

# Set the bivariate Grid Policy ( UNIT_BASE_CORRELATED CORRELATED UNIT_BASE )
grid_policys=( UNIT_BASE_CORRELATED )

# Set the elastic distribution mode ( DECOUPLED COUPLED HYBRID )
modes=( COUPLED )

# Set the elastic coupled sampling method
# ( ONE_D TWO_D MODIFIED_TWO_D )
methods=( MODIFIED_TWO_D )

##---------------------------------------------------------------------------##
## ------------------------------- COMMANDS ---------------------------------##
##---------------------------------------------------------------------------##

# Set the number of threads
command="s/\#SBATCH[[:space:]]--ntasks=.*/\#SBATCH --ntasks=${MPI_PROCESSES}/"
sed -i "${command}" mclaughlin.sh
command="s/\#SBATCH[[:space:]]--cpus-per-task=.*/\#SBATCH --cpus-per-task=${OPEN_MP_THREADS}/"
sed -i "${command}" mclaughlin.sh

command=s/TIME=.*/TIME=${TIME}/
sed -i "${command}" mclaughlin.sh
command=s/HISTORIES=.*/HISTORIES=${HISTORIES}/
sed -i "${command}" mclaughlin.sh

for material in "${materials[@]}"
do
  # Set the material
  command="s/MATERIAL=.*/MATERIAL=\"${material}\"/"
  sed -i "${command}" mclaughlin.sh
  echo "Setting material to ${material}"

  for file_type in "${file_types[@]}"
  do
    # Set the file type
    command=s/FILE_TYPE=.*/FILE_TYPE=${file_type}/
    sed -i "${command}" mclaughlin.sh
    echo "  Setting file type to ${file_type}"

    if [ "${file_type}" = "Native" ]; then

      for interp in "${interps[@]}"
      do
        # Set the interp
        command=s/INTERP=.*/INTERP=${interp}/
        sed -i "${command}" mclaughlin.sh
        echo "    Setting interpolation to ${interp}"

        for grid_policy in "${grid_policys[@]}"
        do
          if [ "${interp}" == "LINLINLOG" ] && [ "${grid_policy}" == "CORRELATED" ]; then
            echo "      The interp (${interp}) and grid policy (${grid_policy}) combo will be skipped."
          else
            # Set 2D grid policy
            command=s/GRID_POLICY=.*/GRID_POLICY=${grid_policy}/
            sed -i "${command}" mclaughlin.sh
            echo "      Setting grid policy to ${grid_policy}"

            for mode in "${modes[@]}"
            do
              # Set the elastic distribution mode
              command=s/MODE=.*/MODE=${mode}/
              sed -i "${command}" mclaughlin.sh
              echo "        Setting elastic mode to ${mode}"

              if [ "${mode}" == "COUPLED" ]; then

                for method in "${methods[@]}"
                do
                  # Set the elastic coupled sampling method
                  command=s/METHOD=.*/METHOD=${method}/
                  sed -i "${command}" mclaughlin.sh
                  echo "          Setting elastic coupled sampling method to ${method}"

                  echo -e "            Running McLaughlin ${material}!\n"
                  sbatch mclaughlin.sh
                done
              else
                echo -e "          Running McLaughlin ${material}!\n"
                sbatch mclaughlin.sh
              fi
            done
          fi
        done
      done
    else
      echo -e "    Running McLaughlin ${material}!\n"
      sbatch mclaughlin.sh
    fi

  done
done