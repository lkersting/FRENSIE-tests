## Sanford experiment ##

# Experimental
0.75 MeV electrons are incident on a carbon convertor (1.76 g/cm3).
The convertor is a cylinder 0.48 cm thick, and 1 cm in radius.
Electrons are normally incident in the direction of the cylinder axis, centered,
and with a radius of 0.1 cm.
A CaF2 detector is down-stream from the electron source by 30.22 cm, sandwiched
by 0.22 cm thick aluminum equilibration layers on each side.
There is air between the convertor and detector, at scaled density for
Albuquerque.

CaF2 = 48.668% F + 51.332% Ca, density = 3.180 g/cm3
Air = 0.0124% C + 75.5268% N + 23.1781% O + 1.2827% Ar, density = 0.0009988815 g/cm3

# Setup
Energy Deposition Tally in the CaF2 in 2.5 cm radial regions out to 60 cm.

# Trelis geometry commands
## C Converter
cylinder z 0.48 radius 1.0
move volume 1 x 0.0 y 0.0 z 0.24
## CaF2 Detector
cylinder z 0.89 radius 60.0
cylinder z 0.89 radius 57.5
cylinder z 0.89 radius 55.0
cylinder z 0.89 radius 52.5
cylinder z 0.89 radius 50.0
cylinder z 0.89 radius 47.5
cylinder z 0.89 radius 45.0
cylinder z 0.89 radius 42.5
cylinder z 0.89 radius 40.0
cylinder z 0.89 radius 37.5
cylinder z 0.89 radius 35.0
cylinder z 0.89 radius 32.5
cylinder z 0.89 radius 30.0
cylinder z 0.89 radius 27.5
cylinder z 0.89 radius 25.0
cylinder z 0.89 radius 22.5
cylinder z 0.89 radius 20.0
cylinder z 0.89 radius 17.5
cylinder z 0.89 radius 15.0
cylinder z 0.89 radius 12.5
cylinder z 0.89 radius 10.0
cylinder z 0.89 radius 7.5
cylinder z 0.89 radius 5.0
cylinder z 0.89 radius 2.5
subtract volume 3 from volume 2 keep
subtract volume 4 from volume 3 keep
subtract volume 5 from volume 4 keep
subtract volume 6 from volume 5 keep
subtract volume 7 from volume 6 keep
subtract volume 8 from volume 7 keep
subtract volume 9 from volume 8 keep
subtract volume 10 from volume 9 keep
subtract volume 11 from volume 10 keep
subtract volume 12 from volume 11 keep
subtract volume 13 from volume 12 keep
subtract volume 14 from volume 13 keep
subtract volume 15 from volume 14 keep
subtract volume 16 from volume 15 keep
subtract volume 17 from volume 16 keep
subtract volume 18 from volume 17 keep
subtract volume 19 from volume 18 keep
subtract volume 20 from volume 19 keep
subtract volume 21 from volume 20 keep
subtract volume 22 from volume 21 keep
subtract volume 23 from volume 22 keep
subtract volume 24 from volume 23 keep
subtract volume 25 from volume 24 keep
delete volume 2 to 24
move volume 25 to 48 x 0.0 y 0.0 z 30.665
## Aluminum Equilibration Layers
cylinder z 0.22 radius 60.0
Volume 49 copy move z 1.11
move volume 49 to 50 x 0.0 y 0.0 z 30.11
## Air
cylinder z 1.33 radius 60.0
move volume 51 x 0.0 y 0.0 z 30.665
cylinder z 31.335 radius 60.0
move volume 52 x 0.0 y 0.0 z 15.6625
subtract volume 1 from volume 52 keep
delete volume 52
subtract volume 51 from volume 53
## Termination Cell
cylinder z 32.0 radius 60.5
move volume 55 x 0.0 y 0.0 z 15.5
cylinder z 33.0 radius 61.0
move volume 56 x 0.0 y 0.0 z 15.5
subtract volume 55 from volume 56

imprint body all
merge tol 5e-7
merge all
group "termination.cell" add vol 57
group "material_1_density_-1.76" add vol 1
group "material_2_density_-3.18" add vol 25 to 48
group "material_3_density_-2.7" add vol 49, 50
group "material_3_density_-0.0009988815" add vol 54
group "estimator_1.cell.pulse.height.e" add vol 25 to 48
group "estimator_1.cell.pulse.height.p" add vol 25 to 48
export dagmc "path-to-frensie-tests/Sanford/geom.h5m" faceting_tolerance 1.e-5 make_watertight

# Running the simulation if FRENSIE

# Running the simulation if MCNP6.2

# Plotting results