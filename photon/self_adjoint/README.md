## Self Adjoint H Sphere Test ##

# Forward Setup
H sphere of density 8.988E-5 g/cm3 and radius 20.0, 0.5, or 0.01 cm surrounded by a void.
0.1, 0.01, or 0.001 MeV isotropic delta source at the center of the sphere.
Measure flux and current on the surface of the sphere and the track flux in the sphere. Set energy bins for all estimators.

# Adjoint Setup
H sphere of density 8.988E-5 g/cm3 and radius 20.0, 0.5, or 0.01 cm surrounded by a void.
cutoff energy to 0.1, 0.01, or 0.001 MeV isotropic isotropic uniform source at the center of the sphere.
Measure adjoint flux and current on the surface of the sphere and the adjoint track flux in the sphere. Set source energy bins for all estimators. Set a delta distribution at the forward source energy as a response function.

# Trelis geometry commands
To construct the geometry run 'construct_geometry.sh' and enter the desired energy.

radius=20.0
length1=45.0
length2=50.0
name=geom_100keV.h5m
tol="1e-3"
elif [ ${energy} = "0.01" ]; then
  Constructing geometry for source energy 0.01 MeV!"
  radius="200.0"
  length1="450.5"
  length2="500.0"
  name="geom_10keV.h5m"
  tol="1e-5"
elif [ ${energy} = "0.001" ]; then
  Constructing geometry for source energy 0.001 MeV!"
  radius="0.01"
  length1="0.05"
  length2="0.1"
  name="geom_1keV.h5m"
  tol="1e-5"

sphere r ${radius}
sphere r 450.5
sphere r ${radius}

# Create termination cell
brick x 450.5 y 450.5 z 450.5
brick x ${length2} y ${length2} z ${length2}
subtract volume 2 from volume 3

# Imprint and merge
imprint body all
merge tol 5e-7
merge all

# Set groups
group 'termination.cell' add vol 4
group 'material_1_density_-0.00008988' add vol 1
group 'estimator_1.surface.flux.p' add surface 1
group 'estimator_2.surface.flux.p*' add surface 1

# export .h5m file
<!-- export dagmc 'geom_100keV.h5m' faceting_tolerance 1e-3 make_watertight -->
export dagmc 'geom_10keV.h5m' faceting_tolerance 1e-4 make_watertight
<!-- export dagmc 'geom_1keV.h5m' faceting_tolerance 1e-4 make_watertight -->

# comment out this line to not automatically exit Trelis
exit

# Running the simulation if FRENSIE

# Running the simulation if MCNP6.2

# Plotting results