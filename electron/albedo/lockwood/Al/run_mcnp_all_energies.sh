#!/bin/bash
##---------------------------------------------------------------------------##
## MCNP6.2 multiple mpi tests runner
##---------------------------------------------------------------------------##
## Validation runs comparing FRENSIE and MCNP.
## This script will run the run_al_mpi.sh script at several different energies:
## .0093, .011, .0134, .0173, .0252, .0415, .0621, .0818 and .102 MeV
##---------------------------------------------------------------------------##


# Set cross_section.xml directory path.
EXTRA_ARGS=$@
MCNP6_1=/home/software/mcnp6.1.1/bin/mcnp611_linux_x86_64_omp
MCNP6_2=/home/software/mcnp6.2/bin/mcnp6
MCNP=$MCNP6_2
TODAY=$(date +%Y-%m-%d)
OUTPUT_DIR="./results/mcnp/${TODAY}/"

THREADS="8"
if [ "$#" -eq 1 ]; then
    # Set the number of threads used
    THREADS="$1"
fi

# Lockwood source energies
energies=( 0.032 0.058 0.084 0.109 0.314 0.521 1.033 )

# Set the input file name
NAME="mcnp.in"
mkdir -p $OUTPUT_DIR

# loop through energies and run mpi script
for i in "${energies[@]}"
do
    pattern="ERG=${i} POS=-20 0 0 DIR=1 VEC=1 0 0 PAR=e"
    sed -i 's,ERG=.*,'"$pattern"',' "mcnp.in"

    OUTPUT="mcnp_${i}"

    echo "Running MCNP:"
    echo "${MCNP} i=${NAME} n=${OUTPUT}. tasks ${THREADS}"
    ${MCNP} i=${NAME} n="${OUTPUT}." tasks ${THREADS}

    NEW_NAME=${OUTPUT_DIR}${OUTPUT}

    # Move output files to test directory
    mv ${OUTPUT}.o ${NEW_NAME}.o
    mv ${OUTPUT}.r ${NEW_NAME}.r
    mv ${OUTPUT}.m ${NEW_NAME}.m

    # Move to output directory
    cd ${OUTPUT_DIR}

    echo "Processing the results:"
    python ../../../../mcnp_data_processor.py -f ${OUTPUT} -e ${i}
    echo "The processed data is located at: ${OUTPUT_DIR}"
    cd ../../../

done
