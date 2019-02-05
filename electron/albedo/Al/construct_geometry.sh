#!/bin/sh
# This file is named contruct_geometry.sh
##---------------------------------------------------------------------------##
## ---------------------- Al Albedo Geometry Constructor --------------------##
##---------------------------------------------------------------------------##
## This script can be used to construct the Trelis geometry.
## To construct a geometry run 'geom_constructor.sh'.
##---------------------------------------------------------------------------##

temp_file=$(mktemp)
name="geom.h5m"
tol="1e-4"


# Create semi-inifite Al slab
echo "brick x 60.0 y 60.0 z 30.0" >> temp_file # Vol 1
echo "move volume 1 x 0.0 y 0.0 z 15.0" >> temp_file

# Create termination cell
echo "brick x 61.0 y 61.0 z 31.0" >> temp_file # Vol 2
echo "brick x 62.0 y 62.0 z 32.0" >> temp_file # Vol 3
echo "subtract volume 2 from volume 3" >> temp_file # Vol 4
echo "move volume 4 x 0.0 y 0.0 z 15.0" >> temp_file

# Imprint and merge
echo "imprint body all" >> temp_file
echo "merge tol 5e-7" >> temp_file
echo "merge all" >> temp_file

# Set groups
echo "group 'termination.cell' add vol 4" >> temp_file
echo "group 'material_1_density_-2.6989' add vol 1" >> temp_file
echo "group 'estimator_1.surface.flux.e' add surface 1 2" >> temp_file
echo "group 'estimator_2.surface.flux.e*' add surface 1 2" >> temp_file
echo "group 'reflecting.surface' add surface 3 to 6" >> temp_file


# export .h5m file
echo "export dagmc '${name}' faceting_tolerance ${tol} make_watertight" >> temp_file

# comment out this line to not automatically exit Trelis
echo "exit" >> temp_file

trelis temp_file

rm temp_file
rm *.jou