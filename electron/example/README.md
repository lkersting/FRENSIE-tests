## Example H Sphere Test ##

# Setup
H sphere of density 8.988E-5 g/cm3 and radius 20.0, 0.5, or 0.01 cm surrounded by a void.
0.1, 0.01, or 0.001 MeV isotropic delta source at the center of the sphere.
Measure flux and current on the surface of the sphere and the track flux in the sphere. Set energy bins for all estimators.

# Trelis geometry commands
To construct the geometry run 'construct_geometry.sh' and enter the desired energy.

## Run the simulation

# Running the simulation for FRENSIE

Set the desired physics option at the top of run_example.sh.
run `run_example.sh` on the UW-Cluster.
Use scp to copy the rendezvous and albedo files from the results directory on
the UW-Cluster to a local computer.

# Running the simulation for MCNP6.2

Set the path to mcnp6.2 in the run_mcnp.sh script.
run `run_mcnp.sh N` where N is the desired number of cores.

## Plotting results
Run `plot_data.py -m mcnp/processed/data/file -a ace/processed/data/file -e endl/processed/data/file -i moments/processed/data/file -o output_name`
