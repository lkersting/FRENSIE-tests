#!/bin/sh
# This file is named update_adjoint_test_files.sh
#SBATCH --partition=pre
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

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

# Update the test files
python ./update_adjoint_test_files.py -d ${database} -g "UnitBaseCorrelated" -v 0
python ./update_adjoint_test_files.py -d ${database} -g "UnitBaseCorrelated" -v 1 --scatter_above_max_mode_off
python ./update_adjoint_test_files.py -d ${database} -g "UnitBase" -v 2
python ./update_adjoint_test_files.py -d ${database} -g "UnitBase" -v 3 --scatter_above_max_mode_off
python ./update_adjoint_test_files.py -d ${database} -g "UnitBase" -v 4 -i "Outgoing Energy"
python ./update_adjoint_test_files.py -d ${database} -g "UnitBase" -v 5 -i "Outgoing Energy Ratio"