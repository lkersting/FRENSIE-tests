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
grid_policies=( "UnitBase" )
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

bold=$(tput bold)
normal=$(tput sgr0)

for ionization in "${ionizations[@]}"
do
  # Set the ionization sampling mode
  echo "Setting electro-ionization sampling to ${bold}${ionization}${normal}"

  for grid_policy in "${grid_policies[@]}"
  do
    # Set the bivariate Grid Policy
    echo "  Setting the bivariate grid policy to ${bold}${grid_policy}${normal}"

    # Set the tolerances
    if [ "${grid_policy}" = "UnitBase" ]; then
      tabular_eval_tol=1e-7
      xs_convergence_tol=1e-4
      brem_convergence_tol=1e-4
      ion_convergence_tol=1e-4
      brem_eval_tol=1e-7
      version=0

      if [ "${ionization}" = "Outgoing Energy" ]; then
        ion_eval_tol=1e-6
        version=$((version + 6))
      else
        ion_eval_tol=1e-7
      fi

    elif [ "${grid_policy}" = "UnitBaseCorrelated" ]; then
      tabular_eval_tol=1e-7
      xs_convergence_tol=1e-4
      brem_convergence_tol=1e-4
      ion_convergence_tol=1e-4
      brem_eval_tol=1e-6
      version=1

      if [ "${ionization}" = "Outgoing Energy" ]; then
        ion_eval_tol=1e-5
        version=$((version + 6))
      else
        ion_eval_tol=1e-6
      fi

    elif [ "${grid_policy}" = "Correlated" ]; then
      tabular_eval_tol=1e-7
      xs_convergence_tol=1e-4
      brem_convergence_tol=1e-4
      ion_convergence_tol=1e-4
      brem_eval_tol=1e-7
      version=2

      if [ "${ionization}" = "Outgoing Energy" ]; then
        ion_eval_tol=1e-4
        version=$((version + 6))
      else
        ion_eval_tol=1e-5
      fi

    else
      echo "The grid policy ${bold}${grid_policy}${normal} is currently not supported!"
    fi

    for nudge_mode in "${nudge_modes[@]}"
    do
      # Set the nudge mode
      echo "    Setting the nudge past max energy mode to ${bold}${nudge_mode}${normal}"

      convergence_tol="--ion_grid_convergence=${ion_convergence_tol} --brem_grid_convergence=${brem_convergence_tol} --xs_grid_convergence=${xs_convergence_tol}"
      eval_tol="--ion_eval_tol=${ion_eval_tol} --brem_eval_tol=${brem_eval_tol} --tabular_evaluation_tol=${tabular_eval_tol}"

      # Update the test file
      if [ "${nudge_mode}" = "on" ]; then
        # Set the version
        echo "      Setting version number to ${bold}${version}${normal}"

        python_command="python ../update_adjoint_test_files.py -d ${DATABASE_PATH} -z 1000 -e 0.01 -g ${grid_policy} -i \"${ionization}\" -v ${version} ${convergence_tol} ${eval_tol}"
        printf "#!/bin/bash\n${python_command}" > update_H_adjoint_temp.sh
        ${sbatch_command} update_H_adjoint_temp.sh
        rm update_H_adjoint_temp.sh
      else
        version=$((version + 3))
        # Set the version
        echo "      Setting version number to ${bold}${version}${normal}"

        python_command="python ../update_adjoint_test_files.py -d ${DATABASE_PATH} -z 1000 -e 0.01 -g ${grid_policy} -i \"${ionization}\" -v ${version} ${convergence_tol} ${eval_tol} --scatter_above_max_mode_off"
        printf "#!/bin/bash\n${python_command}" > update_H_adjoint_temp.sh
        ${sbatch_command} update_H_adjoint_temp.sh
        rm update_H_adjoint_temp.sh
      fi
    done
  done
done
