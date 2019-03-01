#! /usr/bin/env python
from os import path, environ
import argparse as ap
import sys

# Add the parent directory to the path
sys.path.insert(1,path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))
from infinite_medium_simulation_plot import plotInfiniteMediumSimulationSurfaceFlux
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

    if "unit_base" in user_args.forward_rendezvous_file:
      top_ylims = [ [0.0, 10], [0.0, 3], [0.0, 0.4], [0.0, 0.1], [0.0, .03], [0.0, 0.006] ]
      bottom_ylims = [ [0.9, 1.1], [0.85, 1.15], [0.85, 1.15], [0.85, 1.15], [0.8, 1.2], [0.0, 2.0] ]
      xlims = [ [0.0,0.01], [0.0,0.01], [0.0,0.01], [0.0,0.01], [0.0,0.01], [0.0,0.01] ]
      legend_pos = [ (0.95,0.95), (0.95,0.95), (0.95,0.95), (0.95,0.95), (0.95,0.95), (0.95,0.95) ]
    elif "unit_correlated" in user_args.forward_rendezvous_file:
      top_ylims = [ [0.0, 250], [0.0, 70], [0.0, 8], [0.0, 0.4], [0.0, .005], [0.0, 0.0001] ]
      bottom_ylims = [ [0.5, 1.5], [0.4, 1.5], [0.2, 1.5], [0.0, 2.0], [0.0, 2.0], [0.7, 1.3] ]
      xlims = [ [0.0,0.01], [0.0,0.01], [0.0,0.01], [0.0,0.01], [0.0,0.01], [0.0,0.01] ]
      legend_pos = [ (0.95,0.95), (0.95,0.95), (0.95,0.95), (0.95,0.95), (0.95,0.95), (0.95,0.95) ]
    elif "correlated" in user_args.forward_rendezvous_file:
      top_ylims = [ [0.0, 10], [0.0, 3], [0.0, 0.4], [0.0, 0.1], [0.0, .03], [0.0, 0.006] ]
      bottom_ylims = [ [0.9, 1.1], [0.85, 1.15], [0.85, 1.15], [0.85, 1.15], [0.8, 1.2], [0.0, 2.0] ]
      xlims = [ [0.0,0.01], [0.0,0.01], [0.0,0.01], [0.0,0.01], [0.0,0.01], [0.0,0.01] ]
      legend_pos = [ (0.95,0.95), (0.95,0.95), (0.95,0.95), (0.95,0.95), (0.95,0.95), (0.95,0.95) ]

    # output = None
    # if not user_args.output_name is None:
    #   output = user_args.output_name
    # else:
    #   output = user_args.forward_rendezvous_file.split("_rendezvous_")[0]
    #   output = output.split("forward_0.01_loglog_")[-1]

    for i in range(0, len(entity_ids)):
      # Load forward data from file
      manager = Manager.ParticleSimulationManagerFactory( user_args.forward_rendezvous_file ).getManager()
      event_handler = manager.getEventHandler()
      estimator = event_handler.getEstimator( 1 )
      forward_data = estimator.getEntityBinProcessedData( entity_ids[i] )

      # delete manager
      manager = []

      # Load adjoint data from file
      manager = Manager.ParticleSimulationManagerFactory( user_args.adjoint_rendezvous_file ).getManager()
      event_handler = manager.getEventHandler()
      estimator = event_handler.getEstimator( 2 )
      adjoint_data = estimator.getEntityBinProcessedData( entity_ids[i] )

      energy_bins = list(estimator.getSourceEnergyDiscretization())

      # delete manager
      manager = []

      # output_data_name = output + "_" + str(radii[i])

      # Plot the results
      plotInfiniteMediumSimulationSurfaceFlux( forward_data,
                                               adjoint_data,
                                               energy_bins,
                                               None,
                                               radii[i],
                                               top_ylims[i],
                                               bottom_ylims[i],
                                               xlims[i],
                                               legend_pos[i] )
