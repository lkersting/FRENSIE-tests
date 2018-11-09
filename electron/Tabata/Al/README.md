## Tabata aluminum experiment ##

# Experimental
14.9 MeV electrons on aluminum for a 1D charge deposition profile.

# Setup
Tally in 40 uniform subzones of 5 cm.

# Trelis geometry commands
brick x 100.0 y 100.0 z 5.0
move volume 1 x 0.0 y 0.0 z 2.5
Volume 1 copy move z 5.0 repeat 39

brick x 105.0 y 105.0 z 210.0
move volume 41 x 0.0 y 0.0 z 100.0
brick x 110.0 y 110.0 z 220.0
move volume 42 x 0.0 y 0.0 z 100.0
subtract body 41 from body 42

imprint body all
merge tol 5e-7
merge all
group "termination.cell" add vol 43
group "material_1_density_-2.7" add vol 1 to 40
group "estimator_1.cell.pulse.height.e" add vol 1 to 40
group "reflecting.surface" add surface 3 to 6, 9 to 12, 15 to 18, 21 to 24, 27 to 30, 33 to 36, 39 to 42, 45 to 48, 51 to 54, 57 to 60, 63 to 66, 69 to 72, 75 to 78, 81 to 84, 87 to 90, 93 to 96, 99 to 102, 105 to 108, 111 to 114, 117 to 120, 123 to 126, 129 to 132, 135 to 138, 141 to 144, 147 to 150, 153 to 156, 159 to 162, 165 to 168, 171 to 174, 177 to 180, 183 to 186, 189 to 192, 195 to 198, 201 to 204, 207 to 210, 213 to 216, 219 to 222, 225 to 228, 231 to 234, 237 to 240
export dagmc "path-to-Tabata/Al/geom.h5m" faceting_tolerance 1.e-5 make_watertight

# Running the simulation if FRENSIE

# Running the simulation if MCNP6.2

# Plotting results