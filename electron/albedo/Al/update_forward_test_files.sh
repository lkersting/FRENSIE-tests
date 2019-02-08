#!/bin/bash
##---------------------------------------------------------------------------##
## Al Albedo test data updater
##---------------------------------------------------------------------------##

# Sbatch variables
partition=pre
time=1-00:00:00
ntasks=1
cpus=1

##---------------------------------------------------------------------------##
## ------------------------------- COMMANDS ---------------------------------##
##---------------------------------------------------------------------------##

if ! type "sbatch" > /dev/null; then
  echo "no sbatch"
fi
if type "bash" > /dev/null; then
  echo "bash available"
fi

database='/home/lkersting/software/mcnp6.2/MCNP_DATA/database.xml'
# sbatch_command=''
sbatch_command="sbatch --partition=${partition} --time=${time} --ntasks=${ntasks} --cpus-per-task=${cpus}"

if [ ! -f "${database}" ]; then
  database='/home/software/mcnpdata/database.xml'
  sbatch_command=''
fi

if [ -f "${database}" ]; then

    # Update Al data version 0
    printf "Updating the Al native version 0 test data...\n"
    ${sbatch_command} python ../../update_forward_test_files.py --db_name="${database}" -z 13 -g "UnitBaseCorrelatedGrid" -v 0
    if [ $? -eq 0 ]; then
        printf "Al native data version 0 updated successfully!\n\n"
    else
        printf "Al native data version 0 FAILED to update!\n"
        exit 1
    fi

    # Update Al data version 1
    printf "Updating the Al native version 1 test data...\n"
    ${sbatch_command} python ../../update_forward_test_files.py --db_name="${database}" -z 13 -g "UnitBaseGrid" -v 1 --refine_electron_secondary_grids
    if [ $? -eq 0 ]; then
        printf "Al native data version 1 updated successfully!\n\n"
    else
        printf "Al native data version 1 FAILED to update!\n"
        exit 1
    fi

    # Update Al data version 2
    printf "Updating the Al native version 2 test data...\n"
    ${sbatch_command} python ../../update_forward_test_files.py --db_name="${database}" -z 13 -g "UnitBaseCorrelatedGrid" -v 2 --refine_electron_secondary_grids
    if [ $? -eq 0 ]; then
        printf "Al native data version 2 updated successfully!\n\n"
    else
        printf "Al native data version 2 FAILED to update!\n"
        exit 1
    fi

    # Update Al data version 3
    printf "Updating the Al native version 3 test data...\n"
    ${sbatch_command} python ../../update_forward_test_files.py --db_name="${database}" -z 13 -g "CorrelatedGrid" -v 3 --refine_electron_secondary_grids
    if [ $? -eq 0 ]; then
        printf "Al native data version 3 updated successfully!\n\n"
    else
        printf "Al native data version 3 FAILED to update!\n"
        exit 1
    fi

else
    printf "\nERROR: Invalid database file: ${database}\n"
fi