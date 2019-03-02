#!/bin/bash
# This file is named update_adjoint_test_files.sh

##---------------------------------------------------------------------------##
## ------------------------- adjoint test file updater ----------------------##
##---------------------------------------------------------------------------##
## This scripts generates the the adjoint tests files and adds them to the
## database
##---------------------------------------------------------------------------##
EXTRA_ARGS=$@

# Set the bivariate Grid Policy ( "UnitBase" "UnitBaseCorrelated" "Correlated" )
grid_policies=( "UnitBase" )

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

for grid_policy in "${grid_policies[@]}"
do
  # Set the bivariate Grid Policy
  echo "Setting the bivariate grid policy to ${bold}${grid_policy}${normal}"

  # Set the tolerances
  if [ "${grid_policy}" = "UnitBase" ]; then
    tabular_eval_tol=1e-7
    xs_convergence_tol=1e-4
    brem_convergence_tol=1e-4
    ion_convergence_tol=1e-3
    brem_eval_tol=1e-7
    version=0

  elif [ "${grid_policy}" = "UnitBaseCorrelated" ]; then
    tabular_eval_tol=5e-6
    xs_convergence_tol=1e-4
    brem_convergence_tol=1e-4
    ion_convergence_tol=1e-3
    brem_eval_tol=1e-6
    version=1

  elif [ "${grid_policy}" = "Correlated" ]; then
    tabular_eval_tol=5e-6
    xs_convergence_tol=5e-3
    brem_convergence_tol=5e-3
    ion_convergence_tol=5e-3
    brem_eval_tol=1e-7
    version=2

  else
    echo "The grid policy ${bold}${grid_policy}${normal} is currently not supported!"
  fi

  convergence_tol="--ion_grid_convergence=${ion_convergence_tol} --brem_grid_convergence=${brem_convergence_tol} --xs_grid_convergence=${xs_convergence_tol}"
  eval_tol="--ion_eval_tol=${ion_eval_tol} --brem_eval_tol=${brem_eval_tol} --tabular_evaluation_tol=${tabular_eval_tol}"


  # Set the version
  echo "  Setting version number to ${bold}${version}${normal}"

  python_command="python ../../update_adjoint_test_files.py -d ${DATABASE_PATH} -z 13000 -e 1.033 -g ${grid_policy} -v ${version} ${convergence_tol} ${eval_tol}"
  printf "#!/bin/bash\n${python_command}" > update_Al_adjoint_temp${version}.sh
  ${sbatch_command} update_Al_adjoint_temp${version}.sh
  rm update_Al_adjoint_temp${version}.sh

done
