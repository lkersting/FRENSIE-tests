## McLaughlin polystyrene experiment ##

# Experimental
100 keV electrons on polystyrene for a 1D dose profile
The polystyrene consists of H (7.7418%) and C (92.2582%) and has a density
of 1.06 g/cm3.

# Setup
Tally in 50 uniform subzones over 0.02 cm (0.0004 cm each).

# Trelis geometry commands
brick x 0.05 y 0.05 z 0.0004
move volume 1 x 0.0 y 0.0 z 0.0002
Volume 1 copy move z 0.0004 repeat 49

brick x 0.1 y 0.1 z 0.05
move volume 51 x 0.0 y 0.0 z 0.01
brick x 0.2 y 0.2 z 0.1
move volume 52 x 0.0 y 0.0 z 0.01
subtract body 51 from body 52

merge tol 5e-6
imprint tolerant body all
merge tol 5e-7
merge all
group "termination.cell" add vol 53
group "material_1_density_-1.06" add vol 1 to 50
group "estimator_1.cell.pulse.height.e" add vol 1 to 50
group "reflecting.surface" add surface 3 to 6, 9 to 12, 15 to 18, 21 to 24, 27 to 30, 33 to 36, 39 to 42, 45 to 48, 51 to 54, 57 to 60, 63 to 66, 69 to 72, 75 to 78, 81 to 84, 87 to 90, 93 to 96, 99 to 102, 105 to 108, 111 to 114, 117 to 120, 123 to 126, 129 to 132, 135 to 138, 141 to 144, 147 to 150, 153 to 156, 159 to 162, 165 to 168, 171 to 174, 177 to 180, 183 to 186, 189 to 192, 195 to 198, 201 to 204, 207 to 210, 213 to 216, 219 to 222, 225 to 228, 231 to 234, 237 to 240, 243 to 246, 249 to 252, 255 to 258, 261 to 264, 267 to 270, 273 to 276, 279 to 282, 285 to 288, 291 to 294, 297 to 300
export dagmc "path-to-McLaughlin/polystyrene/geom.h5m" faceting_tolerance 1.e-5 make_watertight

# Running the simulation if FRENSIE

# Running the simulation if MCNP6.2

# Plotting results