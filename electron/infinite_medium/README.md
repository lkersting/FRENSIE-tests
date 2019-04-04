## Infinite medium experiments ##

# Experimental
The surface flux on concentric spheres around a point source in an infinite medium.

# Forward Setup
Surface flux on concentric spheres around an isotropic point source.
Energy bins are used.

# Adjoint Setup
Surface flux on concentric spheres around an isotropic point source.
Source energy bins are used.

## Run the simulation

Set the desired physics option at the top of run_infinite_medium_simulation.sh.
run `run_infinite_medium_simulation.sh` on the UW-Cluster.
Use scp to copy the rendezvous and albedo files from the results directory on
the UW-Cluster to a local computer.


## Plotting results

Move to the directory for the desired physics (e.g. `cd ./H/0.01_uniform/all`.)
Run `all-plot.py -o output_name -f /path/to/forward/rendezvous/file.xml -a /path/to/adjoint/rendezvous/file.xml`
The python script name (e.g. all-plot.py) should be replaced with the appropriate plot script for the selected physics
