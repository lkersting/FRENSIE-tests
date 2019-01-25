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
# database='/home/software/mcnpdata/database.xml'

# Set the version number
version=7

# Set the grid policy ('UnitBaseCorrelated', 'UnitBase')
grid_policy='UnitBase'

# Set the ionization sampling mode ('Knock-On', 'Outgoing Energy')
ionization='Outgoing Energy'

# Turn scatter above max on/off
scatter_above_max_mode='on'

# Update the test file
if [ "${scatter_above_max_mode}" = "on" ]; then
  python ./update_adjoint_test_files.py -d ${database} -g ${grid_policy} -i ${ionization} -v ${version}
else
  python ./update_adjoint_test_files.py -d ${database} -g ${grid_policy} -i ${ionization} -v ${version} --scatter_above_max_mode_off
fi
