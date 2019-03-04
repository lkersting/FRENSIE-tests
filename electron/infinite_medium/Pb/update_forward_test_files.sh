#!/bin/bash
##---------------------------------------------------------------------------##
## Pb test forward data updater
##---------------------------------------------------------------------------##

# Sbatch variables
partition=pre
time=1-00:00:00
ntasks=1
cpus=1

# Set the grid policies ( UnitBase UnitBaseCorrelated Correlated )
grid_policies=( UnitBase )

##---------------------------------------------------------------------------##
## ------------------------------- COMMANDS ---------------------------------##
##---------------------------------------------------------------------------##

sbatch_command="sbatch --partition=${partition} --time=${time} --ntasks=${ntasks} --cpus-per-task=${cpus}"
# sbatch_command=bash
if ! type sbatch > /dev/null 2>&1; then
  sbatch_command=bash
fi

echo "Setting Grid Refinement Off"
# Update Pb data version 0
python_command="python ../../update_forward_test_files.py --db_name="${DATABASE_PATH}" -z 82 -g 'UnitBaseGrid' -v 0"
printf "#!/bin/bash\n${python_command}" > update_Pb_0_temp.sh
${sbatch_command} update_Pb_0_temp.sh
if [ ! $? -eq 0 ]; then
    printf "\nPb native data version 0 FAILED to update!\n"
    rm update_Pb_0_temp.sh
    exit 1
fi
rm update_Pb_0_temp.sh

echo "Setting Grid Refinement On"

for grid_policy in "${grid_policies[@]}"
do
  # Set the grid policy
  echo "  Setting the bivariate grid policy to ${bold}${grid_policy}${normal}"

  # Set the version
  if [ "${grid_policy}" = "UnitBase" ]; then
    version=1
  elif [ "${grid_policy}" = "UnitBaseCorrelated" ]; then
    version=2
  elif [ "${grid_policy}" = "Correlated" ]; then
    version=3

  # Update Pb data version
  python_command="python ../../update_forward_test_files.py --db_name="${DATABASE_PATH}" -z 82 -g '${grid_policy}Grid' -v ${version} --refine_electron_secondary_grids"
  printf "#!/bin/bash\n${python_command}" > update_Pb_${version}_temp.sh
  ${sbatch_command} update_Pb_${version}_temp.sh
  if [ ! $? -eq 0 ]; then
      printf "\nPb native data version ${version} FAILED to update!\n"
      rm update_Pb_${version}_temp.sh
      exit 1
  fi
  rm update_Pb_${version}_temp.sh
done
