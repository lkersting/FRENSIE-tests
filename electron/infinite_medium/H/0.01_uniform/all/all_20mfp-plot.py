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

    entity_id = 19
    radius = 40

    if "no_excitation" in user_args.forward_rendezvous_file:
      top_ylims = [0.0, 0.0075]
      bottom_ylims = [0.93, 1.05]
      xlims = [0.0,0.01]
      legend_pos = 1
    else:
      top_ylims = [0.0, 0.0067]
      bottom_ylims = [0.93, 1.035]
      xlims = [0.0,0.01]
      legend_pos = 1

    output = None
    if not user_args.output_name is None:
      output = user_args.output_name
    else:
      output = user_args.forward_rendezvous_file.split("forward_H_")[0] + 'H_0.01_uniform_all_20mfps'

    # Load forward data from file
    manager = Manager.ParticleSimulationManagerFactory( user_args.forward_rendezvous_file ).getManager()
    event_handler = manager.getEventHandler()
    estimator = event_handler.getEstimator( 1 )
    forward_data = estimator.getEntityBinProcessedData( entity_id )

    # delete manager
    manager = []

    # Load adjoint data from file
    manager = Manager.ParticleSimulationManagerFactory( user_args.adjoint_rendezvous_file ).getManager()
    event_handler = manager.getEventHandler()
    estimator = event_handler.getEstimator( 2 )
    adjoint_data = estimator.getEntityBinProcessedData( entity_id )

    energy_bins = list(estimator.getSourceEnergyDiscretization())

    # delete manager
    manager = []

    output_data_name = output

    # Plot the results
    plotInfiniteMediumSimulationSurfaceFlux( forward_data,
                                             adjoint_data,
                                             energy_bins,
                                             output_data_name,
                                             radius,
                                             "H",
                                             top_ylims,
                                             bottom_ylims,
                                             xlims,
                                             legend_pos )
