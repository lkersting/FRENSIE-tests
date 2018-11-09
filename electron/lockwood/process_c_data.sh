#!/bin/bash
##----------------------------------------------------------------------------##
## ---------------------------- FACEMC test runner ---------------------------##
##----------------------------------------------------------------------------##
## FRENSIE benchmark test: Lockwood dose depth data.
## The dose depth for several ranges and material is calculated by dividing the
## energy deposition by the calorimeter thickness (g/cm^2).
##----------------------------------------------------------------------------##
EXTRA_ARGS=$@

# Run from the rendezvous
if [ "$#" -eq 1 ]; then
  # Set the rendezvous
  RENDEZVOUS="$1"

  echo "Enter the test number:"
  read test_number

  # ranges for 1.0 MeV source (g/cm2)
  ranges=( 0.007805 0.060883 0.112662 0.164441 0.215082 0.268568 0.323192 0.382937 0.444389 0.502427 )
  calorimeter_thickness=1.561e-02
  range=${ranges[$test_number]}

  # Restart the simulation
  echo "Processing ${RENDEZVOUS}!"
  cd ./C/C_1/dagmc
  python -c "from os import path; import sys; sys.path.insert(1,'../../../'); import lockwood; lockwood.processDataFromRendezvous(\"${RENDEZVOUS}\", ${range}, ${calorimeter_thickness} )"

# Run new simulation
else
  echo "Error, bad RENDEZVOUS, ${RENDEZVOUS}"
fi
