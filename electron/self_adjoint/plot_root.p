reset
clear

track_flux='track_flux'

set terminal pngcairo nocrop enhanced font 'Verdana,10' size 1000,500
set xlabel 'Energy (MeV)' 
set ylabel 'Flux' 

set output track_flux.'.png'
set title 'Adjoint Flux in Hydrogen sphere of radius 0.002 cm' 
plot track_flux.'.txt' title '' with errorbars
clear
