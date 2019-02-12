#! /usr/bin/env python
from os import path
from optparse import *
import sys

# Add the parent directory to the path
sys.path.insert(1,path.dirname(path.dirname(path.abspath(__file__))))
from albedo_simulation_plot import plotAlbedoSimulationSpectrum

dir=path.dirname(path.abspath(__file__))

if __name__ == "__main__":

    # Parse the command line arguments
    parser = OptionParser()
    parser.add_option("--forward_file", type="string", dest="forward_file",
                      help="the rendezvous file to load")
    parser.add_option("--adjoint_file", type="string", dest="adjoint_file",
                      help="the rendezvous file to load")
    options,args = parser.parse_args()

    forward_filename = dir + "/" + options.forward_file
    adjoint_filename = dir + "/" + options.adjoint_file

    top_ylims = [0.0, 0.8]
    bottom_ylims = [0.0, 3.0]
    legend_pos = (0.95,0.95)

    # Plot the spectrum
    plotAlbedoSimulationSpectrum( forward_filename,
                                  adjoint_filename,
                                  True,
                                  top_ylims,
                                  bottom_ylims,
                                  None,
                                  legend_pos )
