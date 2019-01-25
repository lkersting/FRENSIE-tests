#!/bin/bash
##---------------------------------------------------------------------------##
## -------------------------- Adjoint File Updater --------------------------##
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## ------------------------------- VARIABLES --------------------------------##
##---------------------------------------------------------------------------##
EXTRA_ARGS=$@

# Set the ionization sampling mode ( "Knock-On" "Outgoing Energy" )
ionizations=( "Knock-On" "Outgoing Energy" )
# Set the bivariate Grid Policy ( "UnitBaseCorrelated" "UnitBase" )
grid_policies=( "UnitBaseCorrelated" "UnitBase" )
# Set the nudge past max energy mode on/off ( 'on' 'off' )
nudge_modes=( 'off' 'on')

##---------------------------------------------------------------------------##
## ------------------------------- COMMANDS ---------------------------------##
##---------------------------------------------------------------------------##

version=0

for ionization in "${ionizations[@]}"
do
  # Set the ionization sampling mode
  command=s/ionization=.*/ionization=\'${ionization}\'/
  sed -i "${command}" update_adjoint_test_files.sh
  echo "Setting ionization sampling mode to ${ionization}"

  for grid_policy in "${grid_policies[@]}"
  do
    # Set the bivariate Grid Policy
    command=s/grid_policy=.*/grid_policy=\'${grid_policy}\'/
    sed -i "${command}" update_adjoint_test_files.sh
    echo "  Setting bivariate Grid Policy to ${grid_policy}"

    for nudge_mode in "${nudge_modes[@]}"
    do
      # Set the nudge mode
      command=s/scatter_above_max_mode=.*/scatter_above_max_mode=\'${nudge_mode}\'/
      sed -i "${command}" update_adjoint_test_files.sh
      echo "    Setting nudge past max energy to ${nudge_mode}"

      # Set the version
      command=s/version=.*/version=${version}/
      sed -i "${command}" update_adjoint_test_files.sh
      echo "     Setting version number to ${version}"

      sbatch update_adjoint_test_files.sh

      # Update the version number
      version=$((version + 1))

    done
  done
done
