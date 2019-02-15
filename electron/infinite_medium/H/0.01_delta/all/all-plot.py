#! /usr/bin/env python
from os import path
import argparse as ap
import sys

# Add the parent directory to the path
sys.path.insert(1,path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))
from infinite_medium_simulation_plot import plotInfiniteMediumSimulationSurfaceFlux

dir=path.dirname(path.abspath(__file__))

if __name__ == "__main__":

    # Set up the argument parser
    description = "This script asks for processed forward and adjoint surface "\
                  "flux result files to plots against each other."

    parser = ap.ArgumentParser(description=description)

    parser.add_argument("-f", "--forward_file", dest="forward_file",
                      help="the processed forward results file", required=True)
    parser.add_argument("-a", "--adjoint_file", dest="adjoint_file",
                      help="he processed adjoint results file", required=True)
    parser.add_argument("-o", "--output_name", dest="output_name",
                      help="the plot output name", required=False)

    # Parse the user's arguments
    user_args = parser.parse_args()

    top_ylims = [0.0, 4e3]
    bottom_ylims = [0.7, 1.3]
    xlims = [0.0085,0.01]
    legend_pos = (0.95,0.95)

    # top_ylims = None
    # bottom_ylims = None
    # legend_pos = None

    # Plot the results
    plotInfiniteMediumSimulationSurfaceFlux( user_args.forward_file,
                                             user_args.adjoint_file,
                                             user_args.output_name,
                                             top_ylims,
                                             bottom_ylims,
                                             xlims,
                                             legend_pos )
