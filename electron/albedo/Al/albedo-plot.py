#! /usr/bin/env python
from os import path
import argparse as ap
import sys

# Add the parent directory to the path
sys.path.insert(1,path.dirname(path.dirname(path.abspath(__file__))))
from albedo_simulation_plot import plotAlbedoSimulationForwardSpectrum, plotAlbedoSimulationForwardCombined

dir=path.dirname(path.abspath(__file__))

if __name__ == "__main__":

    # Set up the argument parser
    description = "This script asks for albedo rendezvous and data files "\
                  "which it then plots against experimental data."

    parser = ap.ArgumentParser(description=description)

    parser.add_argument("-a", "--angle", dest="source_angle",
                      help="the problem source angle", required=False, default=0.0)
    parser.add_argument("-o", "--output_name", dest="output_name",
                      help="the plot output name", required=False)
    parser.add_argument("forward_files", nargs='*',
                        help="combined forward files or forward rendezvous spectrum files")

    # Parse the user's arguments
    user_args = parser.parse_args()

    top_ylims = [0.1, 0.28]
    bottom_ylims = [0.9, 1.1]
    xlims = [1, 256]
    legend_pos = 1

    # Plot the spectrum
    if '.xml' in user_args.forward_files[0]:
      plotAlbedoSimulationForwardSpectrum( user_args.forward_files,
                                           None,
                                           user_args.source_angle,
                                           True,
                                           user_args.output_name,
                                           top_ylims,
                                           bottom_ylims,
                                           xlims,
                                           legend_pos )

    else:
      plotAlbedoSimulationForwardCombined( user_args.forward_files,
                                           user_args.source_angle,
                                           True,
                                           user_args.output_name,
                                           top_ylims,
                                           xlims,
                                           legend_pos )
