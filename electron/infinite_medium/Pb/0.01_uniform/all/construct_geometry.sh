#!/bin/sh
# This file is named contruct_geometry.sh
##---------------------------------------------------------------------------##
## -------------------- Self Adjoint Geometry Constructor -------------------##
##---------------------------------------------------------------------------##
## This script can be used to construct the Trelis geometry.
## To construct a geometry run 'contruct_geometry.sh'.
##---------------------------------------------------------------------------##

temp_file=$(mktemp)


radius1="1.0"
radius2="2.0"
radius3="5.0"
radius4="10.0"
radius5="20.0"
radius6="40.0"
length1="200.0"
length2="210.0"
length3="220.0"
name="geom.h5m"
tol="5e-4"

# Create inner spheres
echo "sphere r ${radius1}" >> temp_file # Vol 1
echo "sphere r ${radius2}" >> temp_file # Vol 2
echo "sphere r ${radius3}" >> temp_file # Vol 3
echo "sphere r ${radius4}" >> temp_file # Vol 4
echo "sphere r ${radius5}" >> temp_file # Vol 5
echo "sphere r ${radius6}" >> temp_file # Vol 6

# Create inifite medium
echo "brick x ${length1} y ${length1} z ${length1}" >> temp_file # Vol 7
echo "subtract volume 6 from volume 7 keep" >> temp_file # Vol 8
echo "delete volume 7" >> temp_file
echo "subtract volume 5 from volume 6 keep" >> temp_file # Vol 9
echo "delete volume 6" >> temp_file
echo "subtract volume 4 from volume 5 keep" >> temp_file # Vol 10
echo "delete volume 5" >> temp_file
echo "subtract volume 3 from volume 4 keep" >> temp_file # Vol 11
echo "delete volume 4" >> temp_file
echo "subtract volume 2 from volume 3 keep" >> temp_file # Vol 12
echo "delete volume 3" >> temp_file
echo "subtract volume 1 from volume 2 keep" >> temp_file # Vol 13
echo "delete volume 2" >> temp_file

# Create termination cell
echo "brick x ${length2} y ${length2} z ${length2}" >> temp_file # Vol 14
echo "brick x ${length3} y ${length3} z ${length3}" >> temp_file # Vol 15
echo "subtract volume 14 from volume 15" >> temp_file # Vol 16

# Imprint and merge
echo "imprint body all" >> temp_file
echo "merge tol 5e-7" >> temp_file
echo "merge all" >> temp_file

# Set groups
echo "group 'termination.cell' add vol 16" >> temp_file
echo "group 'material_1_density_-0.000000000000000000000000000001064' add vol 1 8 9 10 11 12 13" >> temp_file
echo "group 'estimator_1.surface.flux.e' add surface 1 19 21 23 25 27" >> temp_file
echo "group 'estimator_2.surface.flux.e*' add surface 1 19 21 23 25 27" >> temp_file


# export .h5m file
echo "export dagmc '${name}' faceting_tolerance ${tol} make_watertight" >> temp_file

# comment out this line to not automatically exit Trelis
echo "exit" >> temp_file

trelis temp_file

rm temp_file
rm *.jou