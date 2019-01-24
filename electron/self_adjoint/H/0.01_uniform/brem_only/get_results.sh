#!/bin/bash
##---------------------------------------------------------------------------##
## Cluster results retriever
##---------------------------------------------------------------------------##

delete=""
while getopts "d" opt; do
  case $opt in
    d)
      echo "Results will be deleted from the cluster."
      delete="true"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

INSTALL="/home/lkersting/frensie0.4_debug/tests/electron"
# INSTALL="/home/lkersting/frensie0.4_release/tests/electron"

# Get brem_only self-adjoint results
cd ./results
echo -e "\nGet brem_only self-adjoint results:"
  # Copy results to this location
  scp -r aci2:${INSTALL}/self_adjoint/H/0.01_uniform/brem_only/results/* ./

  # Erase files from cluster
  if [ "$delete" = "true" ]; then
    ssh aci2 "rm -rf ${INSTALL}/self_adjoint/H/0.01_uniform/brem_only/results/*"
  fi
cd ../