#!/bin/bash
##---------------------------------------------------------------------------##
## Adjoint test forward data updater
##---------------------------------------------------------------------------##

# Set the database directory path.
while getopts d: option
do case "${option}"
   in
       d) database_directory=${OPTARG};;
   esac
done

if [ -d "$database_directory" ]; then

    # Update H data version 0
    printf "Updating the H native version 0 test data...\n"
    python ../update_forward_test_files.py --db_name="$database_directory/database.xml" -z 1 -g "UnitBaseCorrelatedGrid" -v 0
    if [ $? -eq 0 ]; then
        printf "H native data version 0 updated successfully!\n\n"
    else
        printf "H native data version 0 FAILED to update!\n"
        exit 1
    fi

    # Update H data version 1
    printf "Updating the H native version 1 test data...\n"
    python ../update_forward_test_files.py --db_name="$database_directory/database.xml" -z 1 -g "UnitBaseGrid" -v 1 --refine_electron_secondary_grids
    if [ $? -eq 0 ]; then
        printf "H native data version 1 updated successfully!\n\n"
    else
        printf "H native data version 1 FAILED to update!\n"
        exit 1
    fi

else
    printf "\nERROR: Invalid cross section directory!\n"
    printf "  update_test_files.sh -d cross_sectin_directory\n\n"
fi