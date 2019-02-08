#!/bin/bash
# This file is named update_adjoint_test_files.sh

##---------------------------------------------------------------------------##
## ------------------------- adjoint test file updater ----------------------##
##---------------------------------------------------------------------------##
## This scripts generates the the adjoint tests files and adds them to the
## database
##---------------------------------------------------------------------------##
EXTRA_ARGS=$@

# Set the ionization sampling mode ( "Knock-On" "Outgoing Energy" )
ionizations=( "Knock-On" "Outgoing Energy" )
# Set the bivariate Grid Policy ( "UnitBaseCorrelated" "UnitBase" )
grid_policies=( "UnitBaseCorrelated" "UnitBase" )
# Set the nudge past max energy mode on/off ( 'on' 'off' )
nudge_modes=( 'on' 'off' )

# Sbatch variables
partition=pre
time=1-00:00:00
ntasks=1
cpus=5

##---------------------------------------------------------------------------##
## ------------------------------- COMMANDS ---------------------------------##
##---------------------------------------------------------------------------##

database='/home/lkersting/software/mcnp6.2/MCNP_DATA/database.xml'
if [ ! -f "${database}" ]; then
  database='/home/software/mcnpdata/database.xml'
fi

sbatch_command="sbatch --partition=${partition} --time=${time} --ntasks=${ntasks} --cpus-per-task=${cpus}"
# sbatch_command=bash
if ! type "sbatch" > /dev/null; then
  sbatch_command=bash
fi

for ionization in "${ionizations[@]}"
do
  # Set the ionization sampling mode
  echo "Setting ionization sampling mode to ${ionization}"

  for grid_policy in "${grid_policies[@]}"
  do
    # Set the bivariate Grid Policy
    echo "  Setting bivariate Grid Policy to ${grid_policy}"

    # Set the tolerances
    if [ "${grid_policy}" = "UnitBaseCorrelated" ]; then
      xs_convergence_tol=1e-4
      brem_convergence_tol=1e-4
      ion_convergence_tol=1e-3
      brem_eval_tol=1e-6
      version=0

      if [ "${ionization}" = "Outgoing Energy" ]; then
        ion_eval_tol=1e-5
        version=$((version + 4))
      else
        ion_eval_tol=1e-6
      fi

    elif [ "${grid_policy}" = "UnitBase" ]; then
      xs_convergence_tol=1e-4
      brem_convergence_tol=1e-4
      ion_convergence_tol=1e-3
      brem_eval_tol=1e-7
      version=2

      if [ "${ionization}" = "Outgoing Energy" ]; then
        ion_eval_tol=1e-6
        version=$((version + 4))
      else
        ion_eval_tol=1e-7
      fi

    elif [ "${grid_policy}" = "Correlated" ]; then
      xs_convergence_tol=5e-3
      brem_convergence_tol=5e-3
      ion_convergence_tol=5e-3
      brem_eval_tol=1e-7
      version=8

      if [ "${ionization}" = "Outgoing Energy" ]; then
        ion_eval_tol=1e-4
        version=$((version + 4))
      else
        ion_eval_tol=1e-5
      fi

    else
      echo "The grid policy ${grid_policy} is currently not supported!"
    fi

    for nudge_mode in "${nudge_modes[@]}"
    do
      # Set the nudge mode
      echo "    Setting nudge past max energy to ${nudge_mode}"

      convergence_tol="--ion_grid_convergence=${ion_convergence_tol} --brem_grid_convergence=${brem_convergence_tol} --xs_grid_convergence=${xs_convergence_tol}"
      eval_tol="--ion_eval_tol=${ion_eval_tol} --brem_eval_tol=${brem_eval_tol}"

      # Update the test file
      if [ "${nudge_mode}" = "on" ]; then
        # Set the version
        echo "     Setting version number to ${version}"

        python_command="python ../../update_adjoint_test_files.py -d ${database} -z 13000 -g ${grid_policy} -i "${ionization}" -v ${version} ${convergence_tol} ${eval_tol}"
        printf "#!/bin/bash\n${python_command}" > temp.sh
        ${sbatch_command} temp.sh
        rm temp.sh
      else
        version=$((version + 1))
        # Set the version
        echo "     Setting version number to ${version}"

        python_command="python ../../update_adjoint_test_files.py -d ${database} -z 13000 -g ${grid_policy} -i "${ionization}" -v ${version} ${convergence_tol} ${eval_tol} --scatter_above_max_mode_off"
        printf "#!/bin/bash\n${python_command}" > temp.sh
        ${sbatch_command} temp.sh
        rm temp.sh
      fi
    done
  done
done
