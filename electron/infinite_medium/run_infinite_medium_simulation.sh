#!/bin/bash
##---------------------------------------------------------------------------##
## -------------------------- FRENSIE test runner ---------------------------##
##---------------------------------------------------------------------------##
## Run an infinite medium simulation
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

# Materials ( "H" )
materials=( "H" )

# Sources ( "delta" "uniform" )
sources=( "delta" "uniform" )

# Energies ( 0.01 )
energies=( 0.01 )

# Reactions ( "all" "brem_only" "excitation_only" "ionization_only" )
reactions=( "all" "brem_only" "excitation_only" "ionization_only" )

# Desired transport ( "forward" "adjoint" )
transports=( "forward" "adjoint" )

# Desired elastic distribution modes ( "decoupled" "coupled" "hybrid" )
modes=( "coupled" )

# Desired elastic coupled sampling methods ( "modified 2D" "2D" "1D" )
methods=( "2D" )

# Desired bivariate grid policies ( "unit correlated" "unit base" "correlated" )
grid_policies=( "unit correlated" "unit base" "correlated" )

##---------------------------------------------------------------------------##
## ------------------------------- COMMANDS ---------------------------------##
##---------------------------------------------------------------------------##

sbatch_command="sbatch --partition=${partition} --time=${time} --ntasks=${ntasks} --cpus-per-task=${threads}"
run_date=$(date +'%Y-%m-%d')
mv_slurm_command="\nmv slurm-\'${SLURM_JOB_ID}.out\' ./results/${run_date}"

if ! type sbatch > /dev/null 2>&1; then
  sbatch_command=bash
  ntasks=1
fi

bold=$(tput bold)
normal=$(tput sgr0)

# Set the material mode
for material in "${materials[@]}"
do
  echo "Setting the material to ${bold}${material}${normal}"

  # Move to the material directory
  cd ${material}

  # Set the source type
  for source in "${sources[@]}"
  do
    echo "Setting the source energy type to ${bold}${source}${normal}"

    # Set the maximum problem energy
    for energy in "${energies[@]}"
    do
      echo "Setting the maximum problem energy to ${bold}${energy}${normal}"

      # Move to the source type and energy directory
      cd ${energy}_${source}

      # Set the reaction mode
      for reaction in "${reactions[@]}"
      do
        echo "Setting the reaction mode to ${bold}${reaction}${normal}"

        # Move to the reaction directory
        cd ${reaction}

        # Set the transport mode
        for transport in "${transports[@]}"
        do
          echo "  Setting transport mode to ${bold}${transport}${normal}"

          # Set the bivariate Grid Policy
          for grid_policy in "${grid_policies[@]}"
          do
            echo "    Setting bivariate grid policy to ${bold}${grid_policy}${normal}"

            # Set the elastic distribution mode
            for mode in "${modes[@]}"
            do
              echo "      Setting elastic mode to ${bold}${mode}${normal}"

              if [ "${mode}" == "coupled" ]; then

                # Set the elastic coupled sampling method
                for method in "${methods[@]}"
                do
                  echo "        Setting elastic coupled sampling method to ${bold}${method}${normal}"

                  python_command="mpirun -np ${ntasks} python2.7 ${transport}_infinite_medium.py --num_particles=${num_particles} --threads=${threads} --grid_policy=\'${grid_policy}\' --elastic_mode=\'${mode}\' --elastic_method=\'${method}\'"
                  printf "#!/bin/bash\n${python_command}${mv_slurm_command}" > infinite_medium_temp.sh

                  ${sbatch_command} infinite_medium_temp.sh
                  rm infinite_medium_temp.sh
                done
              else
                python_command="mpirun -np ${ntasks} python2.7 ${transport}_infinite_medium.py --num_particles=${num_particles} --threads=${threads} --grid_policy=\'${grid_policy}\' --elastic_mode=\'${mode}\'"
                printf "#!/bin/bash\n${python_command}${mv_slurm_command}" > infinite_medium_temp.sh

                ${sbatch_command} infinite_medium_temp.sh
                rm infinite_medium_temp.sh
              fi
            done
          done
        done
        # Move back to the sorce type and max energy directory
        cd ../
      done
      # Move back to the material directory
      cd ../
    done
  done
  # Move back to the start directory
  cd ../
done
