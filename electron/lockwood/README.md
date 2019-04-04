## Lockwood experiments ##

# Experimental
Electron energy deposition in a film is measured for various incident energies.

# Setup
Weight and energy pulse height estimator is used to measure the energy deposition in the calorimeter.
Thin films are added in front of the calorimeter to change the depth.
Source is normal incident and mono-energetic.
Each depth must be a separate run.

## Run the simulation with FRENSIE

Select the desired material (e.g. 'al')
Set the desired physics option at the top of run script (e.g. run_al_lockwood.sh ).
run `run_al_lockwood.sh` on the UW-Cluster.
Use scp to copy the rendezvous and albedo files from the results directory on
the UW-Cluster to a local computer.

## Run the simulation with MCNP6.2

Move to the directory for the desired material (e.g. `cd ./Al`.)
Set the desired physics option at the top of run script (e.g. run_mcnp_all_ranges.sh ).
run `run_mcnp_all_ranges.sh`.

## Plotting results

Move to the directory for the desired material and energy (e.g. `cd ./Al/Al_0.314`.)
Run `plot_results.py -o output_name /path/to/processed/results/file/a /path/to/processed/results/file/b /path/to/processed/results/file/c ...`
The `-e` flag will plots the designated FRENSIE results against experimental results.
The `-m` flag will plots the designated FRENSIE results against mcnp results.
The `-o` flag can be use to designate the name of the output file.
Multiple FRENSIE spectrums can be designated when plotting.
