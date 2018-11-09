import numpy
import math as m
import matplotlib.pyplot as plt
import os
import sys
import PyFrensie.Geometry.DagMC as DagMC
import PyFrensie.Utility as Utility
import PyFrensie.MonteCarlo as MonteCarlo
import PyFrensie.MonteCarlo.Collision as Collision
import PyFrensie.MonteCarlo.Event as Event
import PyFrensie.MonteCarlo.Manager as Manager
from spectrum_plot_tools import plotSpectralDataWithErrors

##---------------------------------------------------------------------------##
## Plot the frensie and mcnp spectral results
def plotDysonSphereSimulationSpectrum( rendezvous_file,
                                       estimator_id,
                                       entity_id,
                                       mcnp_file,
                                       mcnp_file_start,
                                       mcnp_file_end,
                                       is_a_current,
                                       top_ylims = None,
                                       bottom_ylims = None,
                                       xlims = None,
                                       legend_pos = None ):

    # Set the database path
    Collision.FilledGeometryModel.setDefaultDatabasePath( os.environ['DATABASE_PATH'] )
    
    # Reload the simulation
    manager = Manager.ParticleSimulationManagerFactory( rendezvous_file ).getManager()
    
    # Extract the estimator of interest
    estimator = manager.getEventHandler().getEstimator( estimator_id )

    full_entity_bin_data = estimator.getEntityBinProcessedData( entity_id )

    # start_index = 0
    # end_index = estimator.getNumberOfBins( Event.OBSERVER_ENERGY_DIMENSION )
    start_index = estimator.getNumberOfBins( Event.OBSERVER_ENERGY_DIMENSION )
    end_index = 2*start_index

    entity_bin_data = {"mean": [], "re": [], "e_bins": []}

    for i in range(start_index, end_index):
        entity_bin_data["mean"].append( full_entity_bin_data["mean"][i] )
        entity_bin_data["re"].append( full_entity_bin_data["re"][i] )

    entity_bin_data["e_bins"] = list(estimator.getEnergyDiscretization())

    # Extract the mcnp data from the output file
    mcnp_file = open( mcnp_file, "r" )
    mcnp_file_lines = mcnp_file.readlines()
    
    mcnp_bin_data = {"e_up": [], "mean": [], "re": []}

    mcnp_first_nonzero_index = 0
    mcnp_first_nonzero_index_found = False
    
    for i in range(mcnp_file_start,mcnp_file_end+1):
        split_line = mcnp_file_lines[i-1].split()

        mean_value = float(split_line[1])

        if mean_value == 0.0 and not mcnp_first_nonzero_index_found:
            mcnp_first_nonzero_index += 1
        else:
            mcnp_bin_data["e_up"].append( float(split_line[0]) )
            mcnp_bin_data["mean"].append( mean_value )
            mcnp_bin_data["re"].append( float(split_line[2]) )
            mcnp_first_nonzero_index_found = True

    # Filter out zero values
    del entity_bin_data["e_bins"][0:mcnp_first_nonzero_index]
    del entity_bin_data["mean"][0:mcnp_first_nonzero_index]
    del entity_bin_data["re"][0:mcnp_first_nonzero_index]

    for i in range(0,len(entity_bin_data["e_bins"])-1):
        print entity_bin_data["e_bins"][i+1], entity_bin_data["mean"][i], entity_bin_data["re"][i], " ", mcnp_bin_data["e_up"][i], mcnp_bin_data["mean"][i], mcnp_bin_data["re"][i]

    output_file_name = "h_broomstick_"
    output_file_names = []

    if is_a_current:
        output_file_names.append( output_file_name + "current.eps" )
        output_file_names.append( output_file_name + "current.png" )
        data_type = "Current"
    else:
        output_file_names.append( output_file_name + "flux.eps" )
        output_file_names.append( output_file_name + "flux.png" )
        data_type = "Flux"
        
    # Plot the data
    plotSpectralDataWithErrors( "FRENSIE",
                                entity_bin_data,
                                "MCNP6",
                                mcnp_bin_data,
                                data_type,
                                log_spacing = False,
                                per_lethargy = False,
                                top_ylims = top_ylims,
                                bottom_ylims = bottom_ylims,
                                xlims = xlims,
                                legend_pos = legend_pos,
                                output_plot_names = output_file_names )

##---------------------------------------------------------------------------##
## Extract the mcnp estimator data
def extractMCNPEstimatorData( mcnp_file,
                              mcnp_file_start,
                              mcnp_file_end,
                              output_file_name ):

    # Extract the mcnp data from the output file
    mcnp_file = open( mcnp_file, "r" )
    mcnp_file_lines = mcnp_file.readlines()
    
    mcnp_bin_data = {"e_bins": [], "mean": [], "re": []}

    for i in range(mcnp_file_start,mcnp_file_end+1):
        split_line = mcnp_file_lines[i-1].split()

        mcnp_bin_data["e_bins"].append( float(split_line[0]) )
        
        if i != mcnp_file_start:
            mcnp_bin_data["mean"].append( float(split_line[1]) )
            mcnp_bin_data["re"].append( float(split_line[2]) )

    # Save the data to the desired file
    output_file = open( output_file_name, "w" )

    for i in range(0,len(mcnp_bin_data["mean"])):
        output_file.write(str(mcnp_bin_data["e_bins"][i])+" "+str(mcnp_bin_data["e_bins"][i+1])+" "+str(mcnp_bin_data["mean"][i])+" "+str(mcnp_bin_data["re"][i])+"\n")

##---------------------------------------------------------------------------##
## Extract the frensie estimator data
def extractFrensieEstimatorData( rendezvous_file,
                                 estimator_id,
                                 entity_id,
                                 output_file_name ):
  
    # Set the database path
    Collision.FilledGeometryModel.setDefaultDatabasePath( os.environ['DATABASE_PATH'] )
    
    # Reload the simulation
    manager = Manager.ParticleSimulationManagerFactory( rendezvous_file ).getManager()
    
    # Extract the estimator of interest
    estimator = manager.getEventHandler().getEstimator( estimator_id )
    
    full_entity_bin_data = estimator.getEntityBinProcessedData( entity_id )
    
    start_index = estimator.getNumberOfBins( Event.OBSERVER_ENERGY_DIMENSION )
    end_index = 2*start_index

    entity_bin_data = {"mean": [], "re": [], "e_bins": []}

    for i in range(start_index, end_index):
        entity_bin_data["mean"].append( full_entity_bin_data["mean"][i] )
        entity_bin_data["re"].append( full_entity_bin_data["re"][i] )

    entity_bin_data["e_bins"] = list(estimator.getEnergyDiscretization())

    # Save the data to the desired file
    output_file = open( output_file_name, "w" )

    for i in range(0,len(entity_bin_data["mean"])):
        output_file.write(str(entity_bin_data["e_bins"][i])+" "+str(entity_bin_data["e_bins"][i+1])+" "+str(entity_bin_data["mean"][i])+" "+str(entity_bin_data["re"][i])+"\n")
