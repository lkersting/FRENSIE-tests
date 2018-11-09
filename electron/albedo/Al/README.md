## Aluminum albedo experiment ##

# Experimental
The reflection coefficient in a semi-infinite slab Al is measured at various
energies.

# Setup
Surface current estimator with cosines bin (-1.0, -0.99, 0.0, 1.0).
For Al make the infinite slab 30 cm long

## Trelis geometry commands
brick x 60.0 y 60.0 z 30.0
move volume 1 x 0.0 y 0.0 z 15.0

brick x 61.0 y 61.0 z 31.0
brick x 62.0 y 62.0 z 32.0
subtract volume 2 from volume 3
move volume 4 x 0.0 y 0.0 z 15.0

imprint body all
merge tol 5e-7
merge all
group "termination.cell" add vol 4
group "material_1_density_-2.6989" add vol 1
group "estimator_1.surface.current.e" add surface 1, 2
group "estimator_2.cell.tl.flux.e" add vol 1
group "reflecting.surface" add surface 3 to 6
export dagmc "path-to-albedo/Al/geom.h5m" faceting_tolerance 1.e-5 make_watertight

# Running the simulation if FRENSIE

Set the desired physics option at the top of run_al_albedo.sh.
run `run_al_albedo.sh` on the UW-Cluster.
Use scp to copy the rendezvous and albedo files from the results directory on
the UW-Cluster to a local computer.

# Running the simulation if MCNP6.2

Set the path to mcnp6.2 in the run_mcnp.sh script.
run `run_mcnp.sh N` where N is the desired number of cores.

# Plotting results