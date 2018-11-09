#!/bin/bash
##---------------------------------------------------------------------------##
## MCNP6 test runner
##---------------------------------------------------------------------------##

# Set cross_section.xml directory path.
EXTRA_ARGS=$@
MCNP6=/home/software/mcnp6.1.1/bin/mcnp611_linux_x86_64_omp
CROSS_SECTION_XML_PATH=/home/software/mcnpdata/
TODAY=$(date +%Y-%m-%d)
OUTPUT_DIR="./results/${TODAY}/"

THREADS="8"
if [ "$#" -eq 1 ]; then
    # Set the number of threads used
    THREADS="$1"
fi

# energy in MeV (.005 .0093 .01 .011 .0134 .015 .0173 .02 .0252 .03 .04 .0415 .05 .06 .0621 .0818 .102)
energy=.005

# Set the input file name
NAME="mcnp.in"
mkdir -p $OUTPUT_DIR

pattern="ERG=${energy} POS=-20 0 0 DIR=1 VEC=1 0 0 PAR=e"
sed -i 's,ERG=.*,'"${pattern}"',' "mcnp.in"

OUTPUT="mcnp_${energy}"

echo "Running MCNP 6.2:"
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
python ../../../../mcnp_data_processor.py -f ${OUTPUT} -e ${energy}
echo "The processed data is located at: ${OUTPUT_DIR}"
