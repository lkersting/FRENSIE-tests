## Albedo experiments ##

# Experimental
The reflection coefficient in a semi-infinite slab is measured at various
energies and angles.

# Forward Setup
Surface current estimator with cosines bin (-1.0, 0.0, 1.0).
A semi-infinite slab.
Each incident angle must be a separate run.
Each source energy must be a separate run or source energy bins must be used.

# Adjoint Setup
Surface current estimator with cosines bins around the incident angle and energy bins.
The same geometry is used as the forward.
There can be one adjoint run for all source energies and angles.

## Run the simulation

# Running the simulation for FRENSIE - normal incident angle

Set the desired physics option at the top of run_albedo.sh.
run `run_albedo.sh` on the UW-Cluster.
Use scp to copy the rendezvous and albedo files from the results directory on
the UW-Cluster to a local computer.

# Running the simulation for FRENSIE - multiply incident angles

Set the desired physics option at the top of run_lockwood_albedo.sh.
run `run_lockwood_albedo.sh` on the UW-Cluster.
Use scp to copy the rendezvous and albedo files from the results directory on
the UW-Cluster to a local computer.

# Running the simulation for MCNP6.2

Move to the directory for the desired material (e.g. `cd ./Al`.)
Set the path to mcnp6.2 in the run_mcnp.sh script.
run `run_mcnp.sh N` where N is the desired number of cores.

## Plotting results

# Plotting Forward Discrete Energy Data
Combine the albedo results (FRENSIE or MCNP) for the different source energies into a single txt document.
Run `data_combiner.py /path/to/first/albedo.txt /path/to/second/albedo.txt /path/to/third/albedo.txt ...`
Move to the directory for the desired material (e.g. `cd ./Al`.)
Run `albedo-plot.py -o output_name /path/to/combined/file/a /path/to/combined/file/b /path/to/combined/file/c ...`


# Plotting Forward Spectrum Data
Move to the directory for the desired material (e.g. `cd ./Al`.)
Run `albedo-plot.py -a incident_angle -o output_name /path/to/rendezvous/file/a /path/to/rendezvous/file/b /path/to/rendezvous/file/c ...`


# Plotting Forward vs Adjoint Data
Move to the directory for the desired material (e.g. `cd ./Al`.)
Run `adjoint_albedo-plot.py -f /path/to/forward/rendezvous/file -a /path/to/adjoint/rendezvous/file -o output_name`
