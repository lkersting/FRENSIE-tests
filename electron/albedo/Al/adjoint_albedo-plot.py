#! /usr/bin/env python
from os import path
import argparse as ap
import sys

# Add the parent directory to the path
sys.path.insert(1,path.dirname(path.dirname(path.abspath(__file__))))
from albedo_simulation_plot import plotAlbedoSimulationSpectrum

dir=path.dirname(path.abspath(__file__))

if __name__ == "__main__":

    # Set up the argument parser
    description = "This script asks for albedo rendezvous and data files "\
                  "which it then plots against experimental data."

    parser = ap.ArgumentParser(description=description)

    parser.add_argument("-f", "--forward_file", dest="forward_file",
                        help="the forward rendezvous file to load", required=True)
    parser.add_argument("-a", "--adjoint_file", dest="adjoint_file",
                      help="the rendezvous file to load", required=True)
    parser.add_argument("-o", "--output_name", dest="output_name",
                      help="the plot output name", required=False)
    parser.add_argument("combined_forward_files", nargs='*',
                        help="the combined forward discrete file to load")

    # Parse the user's arguments
    user_args = parser.parse_args()

    forward_filename = dir + "/" + user_args.forward_file
    adjoint_filename = dir + "/" + user_args.adjoint_file

    if user_args.combined_forward_files:
      combined_forward_files = dir + "/" + user_args.combined_forward_files
    else:
      combined_forward_files = None

    top_ylims = [0.0, 0.6]
    bottom_ylims = [0.5, 2.0]
    legend_pos = (0.95,0.95)

    exp_files = ['assad', 'bienlein','bishop', 'bongeler', 'bronshtein', 'cosslett', 'drescher', 'el_gomati', 'heinrich', 'kanter', 'kulenkampff', 'lockwood', 'neubert', 'reimer', 'shimizu', 'soum', 'trump', 'wittry' ]

    for i in range(len(exp_files)):
      exp_files[i] = dir + "/experimental_results/" + exp_files[i] + ".tsv"
      print exp_files[i]


    # Plot the spectrum
    plotAlbedoSimulationSpectrum( forward_filename,
                                  adjoint_filename,
                                  exp_files,
                                  combined_forward_files,
                                  0.0,
                                  "Al",
                                  user_args.output_name,
                                  top_ylims,
                                  bottom_ylims,
                                  None,
                                  legend_pos )
