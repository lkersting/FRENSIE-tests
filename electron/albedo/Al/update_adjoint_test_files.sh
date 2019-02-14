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
ionizations=( "Knock-On" )
# Set the bivariate Grid Policy ( "UnitBase" "UnitBaseCorrelated" "Correlated" )
grid_policies=( "UnitBase" "UnitBaseCorrelated" "Correlated" )
# Set the nudge past max energy mode on/off ( 'on' 'off' )
nudge_modes=( 'on' )

# Sbatch variables
partition=pre
time=1-00:00:00
ntasks=1
cpus=5

##---------------------------------------------------------------------------##
## ------------------------------- COMMANDS ---------------------------------##
##---------------------------------------------------------------------------##

sbatch_command="sbatch --partition=${partition} --time=${time} --ntasks=${ntasks} --cpus-per-task=${cpus}"
# sbatch_command=bash
if ! type sbatch > /dev/null 2>&1; then
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
    if [ "${grid_policy}" = "UnitBase" ]; then
      xs_convergence_tol=1e-4
      brem_convergence_tol=1e-4
      ion_convergence_tol=1e-3
      brem_eval_tol=1e-7
      version=0

      if [ "${ionization}" = "Outgoing Energy" ]; then
        ion_eval_tol=1e-6
        version=$((version + 6))
      else
        ion_eval_tol=1e-7
      fi

    elif [ "${grid_policy}" = "UnitBaseCorrelated" ]; then
      xs_convergence_tol=1e-4
      brem_convergence_tol=1e-4
      ion_convergence_tol=1e-3
      brem_eval_tol=1e-6
      version=1

      if [ "${ionization}" = "Outgoing Energy" ]; then
        ion_eval_tol=1e-5
        version=$((version + 6))
      else
        ion_eval_tol=1e-6
      fi

    elif [ "${grid_policy}" = "Correlated" ]; then
      xs_convergence_tol=5e-3
      brem_convergence_tol=5e-3
      ion_convergence_tol=5e-3
      brem_eval_tol=1e-7
      version=2

      if [ "${ionization}" = "Outgoing Energy" ]; then
        ion_eval_tol=1e-4
        version=$((version + 6))
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
        echo "      Setting version number to ${version}"

        python_command="python ../../update_adjoint_test_files.py -d ${DATABASE_PATH} -z 13000 -e 0.256 -g ${grid_policy} -i \"${ionization}\" -v ${version} ${convergence_tol} ${eval_tol}"
        printf "#!/bin/bash\n${python_command}" > update_Al_adjoint_temp.sh
        ${sbatch_command} update_Al_adjoint_temp.sh
        rm update_Al_adjoint_temp.sh
      else
        version=$((version + 3))
        # Set the version
        echo "      Setting version number to ${version}"

        python_command="python ../../update_adjoint_test_files.py -d ${DATABASE_PATH} -z 13000 -e 0.256 -g ${grid_policy} -i \"${ionization}\" -v ${version} ${convergence_tol} ${eval_tol} --scatter_above_max_mode_off"
        printf "#!/bin/bash\n${python_command}" > update_Al_adjoint_temp.sh
        ${sbatch_command} update_Al_adjoint_temp.sh
        rm update_Al_adjoint_temp.sh
      fi
    done
  done
done
