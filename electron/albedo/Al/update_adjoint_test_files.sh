#!/bin/sh
# This file is named update_adjoint_test_files.sh
#SBATCH --partition=pre
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

##---------------------------------------------------------------------------##
## ------------------------- adjoint test file updater ----------------------##
##---------------------------------------------------------------------------##
## This scripts generates the the adjoint tests files and adds them to the
## database
##---------------------------------------------------------------------------##
EXTRA_ARGS=$@

##---------------------------------------------------------------------------##
## ------------------------------- COMMANDS ---------------------------------##
##---------------------------------------------------------------------------##

database='/home/lkersting/software/mcnp6.2/MCNP_DATA/database.xml'
if [ ! -f "${database}" ]; then
  database='/home/software/mcnpdata/database.xml'
fi

# Set the version number
version=7

# Set the grid policy ('UnitBaseCorrelated', 'UnitBase')
grid_policy='UnitBase'

# Set the ionization sampling mode ('Knock-On', 'Outgoing Energy')
ionization='Outgoing Energy'

# Turn scatter above max on/off
scatter_above_max_mode='off'

if [ "${grid_policy}" = "UnitBaseCorrelated" ]; then
  xs_convergence_tol=1e-4
  brem_convergence_tol=1e-4
  ion_convergence_tol=1e-3
  brem_eval_tol=1e-6

  if [ "${ionization}" = "Outgoing Energy" ]; then
    ion_eval_tol=1e-5
  else
    ion_eval_tol=1e-6
  fi

elif [ "${grid_policy}" = "UnitBase" ]; then
  xs_convergence_tol=1e-4
  brem_convergence_tol=1e-4
  ion_convergence_tol=1e-3
  brem_eval_tol=1e-7

  if [ "${ionization}" = "Outgoing Energy" ]; then
    ion_eval_tol=1e-6
  else
    ion_eval_tol=1e-7
  fi

elif [ "${grid_policy}" = "Correlated" ]; then
  xs_convergence_tol=5e-3
  brem_convergence_tol=5e-3
  ion_convergence_tol=5e-3
  brem_eval_tol=1e-7

  if [ "${ionization}" = "Outgoing Energy" ]; then
    ion_eval_tol=1e-4
  else
    ion_eval_tol=1e-5
  fi

else
  echo "The grid policy ${grid_policy} is currently not supported!"
fi

convergence_tol="--ion_grid_convergence=${ion_convergence_tol} --brem_grid_convergence=${brem_convergence_tol} --xs_grid_convergence=${xs_convergence_tol}"
eval_tol="--ion_eval_tol=${ion_eval_tol} --brem_eval_tol=${brem_eval_tol}"

# Update the test file
if [ "${scatter_above_max_mode}" = "on" ]; then
  python ../../update_adjoint_test_files.py -d ${database} -z 13000 -g ${grid_policy} -i "${ionization}" -v ${version} ${convergence_tol} ${eval_tol}
else
  python ../../update_adjoint_test_files.py -d ${database} -z 13000 -g ${grid_policy} -i "${ionization}" -v ${version}  ${convergence_tol} ${eval_tol} --scatter_above_max_mode_off
fi
