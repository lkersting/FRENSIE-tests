#!/bin/bash
##---------------------------------------------------------------------------##
## Adjoint test forward data updater
##---------------------------------------------------------------------------##

# Sbatch variables
partition=pre
time=1-00:00:00
ntasks=1
cpus=1

##---------------------------------------------------------------------------##
## ------------------------------- COMMANDS ---------------------------------##
##---------------------------------------------------------------------------##

database='/home/lkersting/software/mcnp6.2/MCNP_DATA/database.xml'
# sbatch_command=''
sbatch_command="sbatch --partition=${partition} --time=${time} --ntasks=${ntasks} --cpus-per-task=${cpus}"

if [ ! -f "${database}" ]; then
  database='/home/software/mcnpdata/database.xml'
  sbatch_command=''
fi

if [ -f "${database}" ]; then

    # Update H data version 0
    printf "Updating the H native version 0 test data...\n"
    ${sbatch_command} python ../update_forward_test_files.py --db_name="${database}" -z 1 -g "UnitBaseCorrelatedGrid" -v 0
    if [ $? -eq 0 ]; then
        printf "H native data version 0 updated successfully!\n\n"
    else
        printf "H native data version 0 FAILED to update!\n"
        exit 1
    fi

    # Update H data version 1
    printf "Updating the H native version 1 test data...\n"
    ${sbatch_command} python ../update_forward_test_files.py --db_name="${database}" -z 1 -g "UnitBaseGrid" -v 1 --refine_electron_secondary_grids
    if [ $? -eq 0 ]; then
        printf "H native data version 1 updated successfully!\n\n"
    else
        printf "H native data version 1 FAILED to update!\n"
        exit 1
    fi

    # Update H data version 2
    printf "Updating the H native version 2 test data...\n"
    ${sbatch_command} python ../update_forward_test_files.py --db_name="${database}" -z 1 -g "UnitBaseCorrelatedGrid" -v 2 --refine_electron_secondary_grids
    if [ $? -eq 0 ]; then
        printf "H native data version 2 updated successfully!\n\n"
    else
        printf "H native data version 2 FAILED to update!\n"
        exit 1
    fi

    # Update H data version 3
    printf "Updating the H native version 3 test data...\n"
    ${sbatch_command} python ../update_forward_test_files.py --db_name="${database}" -z 1 -g "CorrelatedGrid" -v 3 --refine_electron_secondary_grids
    if [ $? -eq 0 ]; then
        printf "H native data version 3 updated successfully!\n\n"
    else
        printf "H native data version 3 FAILED to update!\n"
        exit 1
    fi

else
    printf "\nERROR: Invalid database file: ${database}\n"
fi