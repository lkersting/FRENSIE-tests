reset
clear

## ---------- Plot Flux ---------- ##
flux='flux'

set terminal svg enhanced font 'Verdana,10' size 1000,500
set xlabel 'Energy (MeV)'
set ylabel 'Flux'

set title 'Adjoint Flux on surface of Hydrogen sphere of radius 0.002 cm'
set output flux.'.svg'
plot flux.'.txt' title '' with errorbars
clear

## ---------- Plot Current ---------- ##
current='current'
set ylabel 'Current'

set title 'Adjoint Current on surface of Hydrogen sphere of radius 0.002 cm'
set output current.'.svg'
plot current.'.txt' title '' with errorbars
clear

## ---------- Plot Track Length Flux ---------- ##
track_flux='track_flux'
set ylabel 'Flux' 

set title 'Adjoint Flux in Hydrogen sphere of radius 0.002 cm' 
set output track_flux.'.svg'
plot track_flux.'.txt' title '' with errorbars
clear

