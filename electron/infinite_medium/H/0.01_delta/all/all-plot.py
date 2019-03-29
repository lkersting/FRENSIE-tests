#! /usr/bin/env python
from os import path, environ
import argparse as ap
import sys

# Add the parent directory to the path
sys.path.insert(1,path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))
from infinite_medium_simulation_plot import plotAllInfiniteMediumSimulationSurfaceFlux
import infinite_medium_simulation as simulation
import simulation_setup as setup

import PyFrensie.Utility as Utility
import PyFrensie.MonteCarlo as MonteCarlo
import PyFrensie.MonteCarlo.Collision as Collision
import PyFrensie.MonteCarlo.Manager as Manager

dir=path.dirname(path.abspath(__file__))

if __name__ == "__main__":

    # Set up the argument parser
    description = "This script asks for a forward and adjoint rendezvous "\
                  "file to plot surface flux data against each other."

    parser = ap.ArgumentParser(description=description)

    parser.add_argument("-f", "--forward_file", dest="forward_rendezvous_file",
                      help="the forward rendezvous file", required=True)
    parser.add_argument("-a", "--adjoint_file", dest="adjoint_rendezvous_file",
                      help="the adjoint rendezvous file", required=True)
    parser.add_argument("-o", "--output_name", dest="output_name",
                      help="the plot output name", required=False)

    # Parse the user's arguments
    user_args = parser.parse_args()

    # Activate just-in-time initialization to prevent automatic loading of the
    # geometry and data tables
    Utility.activateJustInTimeInitialization()

    # Set the database path
    Collision.FilledGeometryModel.setDefaultDatabasePath( environ['DATABASE_PATH'] )

    entity_ids = [1, 27, 25, 23, 21, 19]
    radii = [1, 2, 5, 10, 20, 40]

    if "no_excitation" in user_args.forward_rendezvous_file:
      top_ylims = [ [0.0, 4e3], [0.0, 750], [0.0, 44], [0.0, 6], [0.0, 0.79], [0.0, 0.13] ]
      bottom_ylims = [ [0.81, 1.29], [0.81, 1.29], [0.81, 1.29], [0.61, 1.29], [0.75, 1.25], [0.85, 1.15] ]
      xlims = [ [0.0093,0.01], [0.008,0.01] ]
      legend_pos = [ 2, 2, 2, 2, 2, 2 ]
    else:
      top_ylims = [ [0.0, 3e3], [0.0, 440], [0.0, 34], [0.0, 6], [0.0, 0.69], [0.0, 0.115] ]
      bottom_ylims = [ [0.3, 1.4], [0.3, 1.4], [0.3, 1.3], [0.0, 1.9], [0.0, 1.9], [0.0, 1.9] ]
      xlims = [ [0.0093,0.01], [0.007,0.01] ]
      legend_pos = [ 2, 2, 2, 2, 2, 2 ]

    output = None
    if not user_args.output_name is None:
      output = user_args.output_name
    else:
      output = user_args.forward_rendezvous_file.split("forward_H_")[0] + 'H_0.01_delta_all'

    for j in range(2):

      forward_data = [None]*3
      adjoint_data = [None]*3

      for i in range(0, 3):
        # Load forward data from file
        manager = Manager.ParticleSimulationManagerFactory( user_args.forward_rendezvous_file ).getManager()
        event_handler = manager.getEventHandler()
        estimator = event_handler.getEstimator( 1 )
        forward_data[i] = estimator.getEntityBinProcessedData( entity_ids[j*3+i] )

        # delete manager
        manager = []

        # Load adjoint data from file
        manager = Manager.ParticleSimulationManagerFactory( user_args.adjoint_rendezvous_file ).getManager()
        event_handler = manager.getEventHandler()
        estimator = event_handler.getEstimator( 2 )
        adjoint_data[i] = estimator.getEntityBinProcessedData( entity_ids[j*3+i] )

        energy_bins = list(estimator.getSourceEnergyDiscretization())

        # delete manager
        manager = []

      output_data_name = output + "_" + str(j)

      # Plot the results
      plotAllInfiniteMediumSimulationSurfaceFlux( forward_data,
                                                  adjoint_data,
                                                  energy_bins,
                                                  output_data_name,
                                                  radii[j*3:j*3+3],
                                                  top_ylims[j*3:j*3+3],
                                                  bottom_ylims[j*3:j*3+3],
                                                  xlims[j],
                                                  legend_pos[j*3:j*3+3] )
