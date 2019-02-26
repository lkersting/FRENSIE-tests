#!/bin/bash
##---------------------------------------------------------------------------##
## ---------------------------- FACEMC test runner --------------------------##
##---------------------------------------------------------------------------##
## FRENSIE verification test: Self Adjoint test.
## The adjoint surface flux in source energy bins at the delta forward source
## energy
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
HISTORIES=1e8
# Set the max runtime (in minutes, 1 day = 1440 )
TIME=1350

# Set the scripts that will be edited ( 'forward.sh' 'adjoint.sh' )
scripts=( 'forward.sh' 'adjoint.sh' )

# Set the elastic distribution mode ( DECOUPLED COUPLED HYBRID )
modes=( DECOUPLED )

# Set the elastic coupled sampling method
# ( ONE_D TWO_D MODIFIED_TWO_D )
methods=( MODIFIED_TWO_D )

# Set the ionization sampling mode
# ( OUTGOING_ENERGY KNOCK_ON )
ionizations=( KNOCK_ON )

# Set the bivariate Grid Policy ( 'UNIT_BASE_CORRELATED' 'UNIT_BASE' )
grid_policies=( 'UNIT_BASE' 'UNIT_BASE_CORRELATED' )

# Set the nudge past max energy mode on/off ( 'on' 'off' )
nudge_modes=( 'on' )

# Turn individual physics options off ( ELASTIC EXCITATION BREM IONIZATION )
reactions_off=( )
# Turn inelastic physics options 'off'
inelastic_reactions=''

##---------------------------------------------------------------------------##
## ------------------------------- COMMANDS ---------------------------------##
##---------------------------------------------------------------------------##

bold=$(tput bold)
normal=$(tput sgr0)

# Set the number of threads
command1="s/\#SBATCH[[:space:]]--ntasks=.*/\#SBATCH --ntasks=${MPI_PROCESSES}/"
command2="s/\#SBATCH[[:space:]]--cpus-per-task=.*/\#SBATCH --cpus-per-task=${OPEN_MP_THREADS}/"

# Set the wall time and number of histories
command3=s/TIME=.*/TIME=${TIME}/
command4=s/HISTORIES=.*/HISTORIES=${HISTORIES}/

for script in ${scripts[@]}
do
  sed -i "${command1}" ${script}
  sed -i "${command2}" ${script}
  sed -i "${command3}" ${script}
  sed -i "${command4}" ${script}
done

for mode in "${modes[@]}"
do
  # Set the elastic distribution mode
  command=s/MODE=.*/MODE=${mode}/
  for script in ${scripts[@]}; do
    sed -i "${command}" ${script}
  done
  echo "Setting elastic mode to ${bold}${mode}${normal}"

  if [ "${mode}" == "COUPLED" ]; then

    for method in "${methods[@]}"
    do
      # Set the elastic coupled sampling method
      command=s/METHOD=.*/METHOD=${method}/
      for script in ${scripts[@]}; do
        sed -i "${command}" ${script}
      done
      echo "  Setting elastic coupled sampling method to ${bold}${method}${normal}"

      for ionization in "${ionizations[@]}"
      do
        # Set the ionization sampling mode
        command=s/IONIZATION_MODE=.*/IONIZATION_MODE=${ionization}/
        for script in ${scripts[@]}; do
          sed -i "${command}" ${script}
        done
        echo "    Setting electro-ionization sampling to ${bold}${ionization}${normal}"

        for grid_policy in "${grid_policies[@]}"
        do
          # Set the bivariate Grid Policy
          command=s/GRID_POLICY=.*/GRID_POLICY=${grid_policy}/
          for script in ${scripts[@]}; do
            sed -i "${command}" ${script}
          done
          echo "      Setting the bivariate grid policy to ${bold}${grid_policy}${normal}"

          for nudge_mode in "${nudge_modes[@]}"
          do
            # Set the nudge mode
            command=s/NUDGE_PAST_MAX=.*/NUDGE_PAST_MAX=\'${nudge_mode}\'/
            sed -i "${command}" adjoint.sh
            echo "        Setting the nudge past max energy mode to ${bold}${nudge_mode}${normal}"

            for reaction in "${reactions_off[@]}"
            do
              # Set the reaction off
              command=s/${reaction}=.*/${reaction}=\'off\'/
              for script in "${scripts[@]}"; do
                sed -i "${command}" ${script}
              done
              echo "          Turning the ${bold}${reaction}${normal} reaction off"

              for script in "${scripts[@]}"; do
                sbatch ${script}
              done

              # Set the reaction on
              command=s/${reaction}=.*/${reaction}=\'\'/
              for script in ${scripts[@]}; do
                sed -i "${command}" ${script}
              done

            done

          done

          if [ "${reactions_off}" == "" ]; then
              for script in "${scripts[@]}"; do
                sbatch ${script}
              done
          fi

          if [ "${inelastic_reactions}" == "off" ]; then
            command_off_1=s/EXCITATION=.*/EXCITATION=\'off\'/
            command_off_2=s/BREM=.*/BREM=\'off\'/
            command_off_3=s/IONIZATION=.*/IONIZATION=\'off\'/

            command_on_1=s/EXCITATION=.*/EXCITATION=\'\'/
            command_on_2=s/BREM=.*/BREM=\'\'/
            command_on_3=s/IONIZATION=.*/IONIZATION=\'\'/

            for script in ${scripts[@]}; do
              sed -i "${command_off_1}" ${script}
              sed -i "${command_off_2}" ${script}
              sed -i "${command_off_3}" ${script}
              sbatch ${script}

              sed -i "${command_on_1}" ${script}
              sbatch ${script}
              sed -i "${command_off_1}" ${script}

              sed -i "${command_on_2}" ${script}
              sbatch ${script}
              sed -i "${command_off_2}" ${script}

              sed -i "${command_on_3}" ${script}
              sbatch ${script}
              sed -i "${command_on_1}" ${script}
              sed -i "${command_on_2}" ${script}

            done
          fi

        done
      done
    done
  else

    for ionization in "${ionizations[@]}"
    do
      # Set the ionization sampling mode
      command=s/IONIZATION_MODE=.*/IONIZATION_MODE=${ionization}/
      for script in ${scripts[@]}; do
        sed -i "${command}" ${script}
      done
      echo "  Setting electro-ionization sampling to ${bold}${ionization}${normal}"

      for grid_policy in "${grid_policies[@]}"
      do
        # Set the bivariate Grid Policy
        command=s/GRID_POLICY=.*/GRID_POLICY=${grid_policy}/
        for script in ${scripts[@]}; do
          sed -i "${command}" ${script}
        done
        echo "    Setting the bivariate grid policy to ${bold}${grid_policy}${normal}"


        for nudge_mode in "${nudge_modes[@]}"
        do
          # Set the nudge mode
          command=s/NUDGE_PAST_MAX=.*/NUDGE_PAST_MAX=\'${nudge_mode}\'/
          sed -i "${command}" adjoint.sh
          echo "      Setting the nudge past max energy mode to ${bold}${nudge_mode}${normal}"

          for reaction in "${reactions_off[@]}"
          do
            # Set the reaction off
            command=s/${reaction}=.*/${reaction}=\'off\'/
            for script in "${scripts[@]}"; do
              sed -i "${command}" ${script}
            done
            echo "        Turning the ${bold}${reaction}${normal} reaction off"

            for script in "${scripts[@]}"; do
              sbatch ${script}
            done

            # Set the reaction on
            command=s/${reaction}=.*/${reaction}=\'\'/
            for script in ${scripts[@]}; do
              sed -i "${command}" ${script}
            done

          done

          if [ "${reactions_off}" == "" ]; then
              for script in "${scripts[@]}"; do
                sbatch ${script}
              done
          fi

          if [ "${inelastic_reactions}" == "off" ]; then
            command_off_1=s/EXCITATION=.*/EXCITATION=\'off\'/
            command_off_2=s/BREM=.*/BREM=\'off\'/
            command_off_3=s/IONIZATION=.*/IONIZATION=\'off\'/

            command_on_1=s/EXCITATION=.*/EXCITATION=\'\'/
            command_on_2=s/BREM=.*/BREM=\'\'/
            command_on_3=s/IONIZATION=.*/IONIZATION=\'\'/

            for script in ${scripts[@]}; do
              sed -i "${command_off_1}" ${script}
              sed -i "${command_off_2}" ${script}
              sed -i "${command_off_3}" ${script}
              sbatch ${script}

              sed -i "${command_on_1}" ${script}
              sbatch ${script}
              sed -i "${command_off_1}" ${script}

              sed -i "${command_on_2}" ${script}
              sbatch ${script}
              sed -i "${command_off_2}" ${script}

              sed -i "${command_on_3}" ${script}
              sbatch ${script}
              sed -i "${command_on_1}" ${script}
              sed -i "${command_on_2}" ${script}

            done
          fi

        done

      done
    done

  fi
done
