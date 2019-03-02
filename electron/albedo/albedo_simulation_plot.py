# import numpy
import math
import matplotlib.pyplot as plt
import os
import PyFrensie.Utility as Utility
import PyFrensie.Geometry.DagMC as DagMC
import PyFrensie.MonteCarlo.Collision as Collision
import PyFrensie.MonteCarlo.Event as Event
import PyFrensie.MonteCarlo.Manager as Manager
from spectrum_plot_tools import plotSpectralDataWithErrors

def plotAlbedoSimulationSpectrum( forward_rendezvous_file,
                                  adjoint_rendezvous_file,
                                  combined_forward_files,
                                  source_angle,
                                  include_experimental = False,
                                  output_plot_name = None,
                                  top_ylims = None,
                                  bottom_ylims = None,
                                  xlims = None,
                                  legend_pos = None ):

    # Activate just-in-time initialization to prevent automatic loading of the
    # geometry and data tables
    Utility.activateJustInTimeInitialization()

    # Set the database path
    Collision.FilledGeometryModel.setDefaultDatabasePath( os.environ['DATABASE_PATH'] )

    # Reload the adjoint simulation
    manager = Manager.ParticleSimulationManagerFactory( adjoint_rendezvous_file ).getManager()

    # Extract the estimator of interest
    estimator = manager.getEventHandler().getEstimator( 2 )
    full_entity_bin_data = estimator.getEntityBinProcessedData( 2 )

    num_bins = estimator.getNumberOfBins( Event.OBSERVER_ENERGY_DIMENSION )

    adjoint_data = {"mean": [None]*num_bins, "re": [None]*num_bins, "e_bins": []}

    # Get the adjoint energy bins
    adjoint_data["e_bins"] = list(estimator.getEnergyDiscretization())

    angles = estimator.getCosineDiscretization()
    ADJOINT_NORM=1.0

    if float(source_angle) == 0.0:
      mu_start=-1
    elif float(source_angle) == 15.0:
      mu_start=-3
    elif float(source_angle) == 30.0:
      mu_start=-5
    elif float(source_angle) == 45.0:
      mu_start=-7
    elif float(source_angle) == 60.0:
      mu_start=-9
    elif float(source_angle) == 75.0:
      mu_start=-11

    start_index = (len(angles)-1-mu_start)*num_bins

    adjoint_respone=(adjoint_data["e_bins"][-1]-adjoint_data["e_bins"][0])*0.5/(angles[mu_start]-angles[mu_start-1])*ADJOINT_NORM
    for i in range(0, num_bins):
        j = start_index + i
        # print j, full_entity_bin_data["mean"][j], full_entity_bin_data["re"][j]
        adjoint_data["mean"][i] = full_entity_bin_data["mean"][j]*adjoint_respone
        adjoint_data["re"][i] = full_entity_bin_data["re"][j]

    manager = None
    # Reload the forward simulation
    manager = Manager.ParticleSimulationManagerFactory( forward_rendezvous_file ).getManager()

    # Extract the estimator of interest
    estimator = manager.getEventHandler().getEstimator( 1 )
    full_entity_bin_data = estimator.getEntityBinProcessedData( 2 )

    num_bins = estimator.getNumberOfBins( Event.OBSERVER_SOURCE_ENERGY_DIMENSION )

    forward_data = {"mean": [None]*num_bins, "re": [None]*num_bins, "e_bins": []}

    # Get the forward energy bins
    forward_data["e_bins"] = list(estimator.getSourceEnergyDiscretization())
    for i in range(0, len(forward_data["e_bins"])):
        if not forward_data["e_bins"][i] == adjoint_data["e_bins"][i]:
          raise ValueError( "The forward and adjoint energy bins must match!")

    FORWARD_NORM=1.0
    forward_respone=(forward_data["e_bins"][-1]-forward_data["e_bins"][0])*FORWARD_NORM
    start_index = (len(estimator.getCosineDiscretization())-2)*num_bins
    for i in range(0, num_bins):
        j = start_index + i
        # print j, full_entity_bin_data["mean"][j], full_entity_bin_data["re"][j]
        forward_data["mean"][i] = full_entity_bin_data["mean"][j]*forward_respone
        forward_data["re"][i] = full_entity_bin_data["re"][j]

    if output_plot_name is None:
      output_plot_name = "al_albedo"
    output_plot_names = []

    output_plot_names.append( output_plot_name + ".eps" )
    output_plot_names.append( output_plot_name + ".png" )

    data_type = "Current"
    forward_data_name = "Forward"
    adjoint_data_name = "Adjoint"

    linestyles = [(0, ()), (0, (5, 5)), (0, (3, 5, 1, 5)), (0, (1, 1)), (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (1, 5)), (0, (5, 10)), (0, (3, 10, 1, 10)), (0, (3, 10, 1, 10, 1, 10))]

    markers = ["o","*","v","^","<",">","+","x","1","2","3","4","p","s","h","D","d","H","8","o","*"]
    exp_names = ['assad', 'bienlein','bishop', 'bongeler', 'bronshtein', 'cosslett', 'drescher', 'el_gomati', 'heinrich', 'kanter', 'kulenkampff', 'lockwood', 'neubert', 'reimer', 'shimizu', 'soum', 'trump', 'wittry' ]

    # Compute the bin norm constants and convert the mean values to mean per energy
    forward_normalized_mean = [None]*num_bins
    forward_error = [None]*num_bins
    adjoint_normalized_mean = [None]*num_bins
    adjoint_error = [None]*num_bins

    for i in range(0, num_bins):
        bin_norm_const = forward_data["e_bins"][i+1] - forward_data["e_bins"][i]
        forward_normalized_mean[i] = forward_data["mean"][i]/bin_norm_const
        forward_error[i] = forward_data["re"][i]*forward_normalized_mean[i]
        # bin_norm_const = adjoint_data["e_bins"][i+1] - adjoint_data["e_bins"][i]
        adjoint_normalized_mean[i] = adjoint_data["mean"][i]/bin_norm_const
        adjoint_error[i] = adjoint_data["re"][i]*adjoint_normalized_mean[i]

    # Compute the F/T values and uncertainties
    f_over_t = []
    f_over_t_unc = []

    for i in range(0, len(forward_normalized_mean)):
        f_over_t.append( forward_normalized_mean[i]/adjoint_normalized_mean[i] )

        sigma_f = forward_normalized_mean[i]*forward_data["re"][i]
        sigma_t = adjoint_normalized_mean[i]*adjoint_data["re"][i]

        f_squared = forward_normalized_mean[i]*forward_normalized_mean[i]
        t_squared = adjoint_normalized_mean[i]*adjoint_normalized_mean[i]

        f_over_t_unc.append( math.sqrt( sigma_f*sigma_f + (f_squared/t_squared)*sigma_t*sigma_t )/adjoint_normalized_mean[i] )

    edge_thickness = 1.1

    title='Electron Albedos for an infinite slab of Al'
    if float(source_angle) == 60.0:
      title+=' 60 Degrees Incident Source'

    # Initialize the plot
    fig, ax = plt.subplots(2, 1, sharex=True)
    plt.subplots_adjust( top=0.95, bottom=0.1, hspace=0.0 )
    ax[0].set_title(title, size=14)

    # Set up the top subplot

    if include_experimental:
      # Plot experimental data
      directory = os.path.dirname(os.path.abspath(__file__))

      for i in range(len(exp_names)):
        filename = directory + "/Al/experimental_results/" + exp_names[i] +".tsv"
        with open(filename) as input:
            name = input.readline().strip()
            input.readline()
            data = zip(*(line.strip().split('\t') for line in input))
            x = [None] * len(data[0][:])
            x = [0 for k in range(len(data[0][:]))]
            y = [0 for k in range(len(data[1][:]))]
            for j in range(len(x)):
              x[j] = float(data[0][j])*1e-3
              y[j] = float(data[1][j])

        if i == 1:
          ax[0].scatter(x, y, label="Experimental", marker=markers[1], s=50, facecolors='none', edgecolors='b' )
        else:
          ax[0].scatter(x, y, marker=markers[1], s=50, facecolors='none', edgecolors='b' )

    if not combined_forward_files == None:
      for i in range(len(combined_forward_files)):
        filename = combined_forward_files[i]
        with open(filename) as input:
            name = input.readline().strip()
            input.readline()
            data = zip(*(line.strip().split('\t') for line in input))
            x = [None] * len(data[0][:])
            x = [0 for k in range(len(data[0][:]))]
            y = [0 for k in range(len(data[1][:]))]
            for j in range(len(x)):
              x[j] = float(data[0][j])
              y[j] = float(data[1][j])

        ax[0].scatter(x, y, label=name, marker=markers[i+2], s=50, facecolors='none', edgecolors='g' )

    # Plot forward histogram of results
    label = forward_data_name
    if not FORWARD_NORM == 1.0:
      label += "*" + str(FORWARD_NORM)
    m_bin, bins, plt1 = ax[0].hist(forward_data["e_bins"][:-1], bins=forward_data["e_bins"], weights=forward_normalized_mean, histtype='step', label=label, color='g', linestyle=linestyles[0], linewidth=1.8 )

    # Plot error bars
    forward_mid = 0.5*(bins[1:] + bins[:-1])
    # plt2 = ax[0].errorbar(forward_mid, m_bin, yerr=forward_error, ecolor='b', fmt=None)

    # Plot adjoint histogram of results
    label = adjoint_data_name
    if not ADJOINT_NORM == 1.0:
      label += "*" + str(ADJOINT_NORM)
    m_bin, bins, plt3 = ax[0].hist(adjoint_data["e_bins"][:-1], bins=adjoint_data["e_bins"], weights=adjoint_normalized_mean, histtype='step', label=label, color='r', linestyle=linestyles[1], linewidth=1.8 )

    # Plot error bars
    adjoint_mid = 0.5*(bins[1:] + bins[:-1])
    # plt4 = ax[0].errorbar(adjoint_mid, m_bin, yerr=adjoint_error, ecolor='g', fmt=None)

    ax[0].set_xscale("log")

    y_label = data_type + " Spectrum"

    ax[0].set_ylabel( y_label )

    if not legend_pos is None:
        ax[0].legend(frameon=True, bbox_to_anchor=legend_pos)
    else:
        ax[0].legend(frameon=True)

    # Turn on the grid
    ax[0].grid(True, linestyle=':', linewidth=1)

    # Set the x limits
    if not xlims is None:
        ax[0].set_xlim( xlims[0], xlims[-1] )
    else:
        ax[0].set_xlim( forward_data["e_bins"][0], forward_data["e_bins"][-1] )

    if not top_ylims is None:
        ax[0].set_ylim( top_ylims[0], top_ylims[1] )

    # Set the y tic labels
    yticklabels = ax[0].yaxis.get_ticklabels()
    # yticklabels[0].set_visible(False)
    yticklabels[-1].set_visible(False)

    # Set the tic properties
    ax[0].yaxis.set_ticks_position("both")
    ax[0].xaxis.set_ticks_position("both")
    ax[0].tick_params(direction="in", width=edge_thickness)
    ax[0].tick_params(which="minor", direction="in", width=edge_thickness)

    for axis in ['top','bottom','left','right']:
        ax[0].spines[axis].set_linewidth(edge_thickness)

    # Set up the bottom subplot
    ax[1].set_xscale("log")

    ax[1].errorbar( forward_mid, f_over_t, yerr=f_over_t_unc, capsize=1.5, fmt='o', ecolor="black", color="black", linewidth=0.5, markersize=1.9 )
    ax[1].set_ylabel( forward_data_name + "/" + adjoint_data_name )
    ax[1].set_xlabel( "Energy (MeV)" )

    # Turn on the grid
    ax[1].grid(True, linestyle=':', linewidth=1)

    # Set the x limits
    if not xlims is None:
        ax[1].set_xlim( xlims[0], xlims[-1] )
    else:
        ax[1].set_xlim( forward_data["e_bins"][0], forward_data["e_bins"][-1] )

    if not bottom_ylims is None:
        ax[1].set_ylim( bottom_ylims[0], bottom_ylims[1] )

    # Set the y tic labels
    yticklabels = ax[1].yaxis.get_ticklabels()
    yticklabels[0].set_visible(False)
    yticklabels[-1].set_visible(False)

    # Set the tic properties
    ax[1].yaxis.set_ticks_position("both")
    ax[1].xaxis.set_ticks_position("both")
    ax[1].tick_params(direction="in", width=edge_thickness)
    ax[1].tick_params(which="minor", direction="in", width=edge_thickness)

    for axis in ['top','bottom','left','right']:
        ax[1].spines[axis].set_linewidth(edge_thickness)

    # Save the figure
    for i in range(0,len(output_plot_names)):
      fig.savefig( output_plot_names[i] )

    plt.show()


def plotAlbedoSimulationForwardSpectrum( forward_rendezvous_files,
                                         source_angle,
                                         include_experimental = False,
                                         output_plot_name = None,
                                         top_ylims = None,
                                         bottom_ylims = None,
                                         xlims = None,
                                         legend_pos = None ):

    # Activate just-in-time initialization to prevent automatic loading of the
    # geometry and data tables
    Utility.activateJustInTimeInitialization()

    # Set the database path
    Collision.FilledGeometryModel.setDefaultDatabasePath( os.environ['DATABASE_PATH'] )

    if output_plot_name is None:
      output_plot_name = "al_albedo"
    output_plot_names = []

    output_plot_names.append( output_plot_name + ".eps" )
    output_plot_names.append( output_plot_name + ".png" )

    data_type = "Current"
    forward_spectrum_data_name = "Forward"

    top_ylims = [0.1, 0.5]

    linestyles = [(0, ()), (0, (5, 5)), (0, (3, 5, 1, 5)), (0, (1, 1)), (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (1, 5)), (0, (5, 10)), (0, (3, 10, 1, 10)), (0, (3, 10, 1, 10, 1, 10))]

    markers = ["o","*","v","^","<",">","+","x","1","2","3","4","p","s","h","D","d","H","8","o","*"]
    exp_names = ['assad', 'bienlein','bishop', 'bongeler', 'bronshtein', 'cosslett', 'drescher', 'el_gomati', 'heinrich', 'kanter', 'kulenkampff', 'lockwood', 'neubert', 'reimer', 'shimizu', 'soum', 'trump', 'wittry' ]
    marker_color = ['g', 'r', 'c', 'm', 'y', 'k', 'g', 'r', 'c', 'm', 'y', 'k' ]

    edge_thickness = 1.1

    # Initialize the plot
    fig, ax = plt.subplots(2, 1, sharex=True)
    plt.subplots_adjust( top=0.95, bottom=0.1, hspace=0.0 )

    title='Electron Albedos for an infinite slab of Al'
    if float(source_angle) == 60.0:
      title+=' 60 Degrees Incident Source'
      include_experimental=False

    ax[0].set_title(title, size=14)

    num_spectrum_files = len(forward_rendezvous_files)

    # Compute the bin norm constants and convert the mean values to mean per energy
    forward_energy_bins = [[]]*num_spectrum_files
    forward_normalized_mean = [[]]*num_spectrum_files
    forward_error = [[]]*num_spectrum_files

    for k in range(0,num_spectrum_files):
      # Reload the forward simulation
      manager = []
      manager = Manager.ParticleSimulationManagerFactory( forward_rendezvous_files[k] ).getManager()

      # Extract the estimator of interest
      estimator = manager.getEventHandler().getEstimator( 1 )
      full_entity_bin_data = estimator.getEntityBinProcessedData( 2 )

      num_bins = estimator.getNumberOfBins( Event.OBSERVER_SOURCE_ENERGY_DIMENSION )

      # Get the forward energy bins
      forward_energy_bins[k] = list(estimator.getSourceEnergyDiscretization())
      if( k > 0 ):
        for i in range(0, num_bins):
          if not forward_energy_bins[k][i] == forward_energy_bins[0][i]:
            raise ValueError( "The forward spectrum energy bins must match!")

      FORWARD_NORM=(forward_energy_bins[k][-1]-forward_energy_bins[k][0])
      start_index = (len(estimator.getCosineDiscretization())-2)*num_bins
      norm_mean = []
      err = []
      for i in range(0, num_bins):
          j = start_index + i
          bin_norm_const = forward_energy_bins[k][i+1] - forward_energy_bins[k][i]
          norm_mean.append( full_entity_bin_data["mean"][j]*FORWARD_NORM/bin_norm_const )
          err.append( full_entity_bin_data["re"][j]*norm_mean[i] )
      forward_normalized_mean[k] = norm_mean
      forward_error[k] = err

    # Compute the F/T values and uncertainties
    f_over_t = [[]]*(num_spectrum_files -1)
    f_over_t_unc = [[]]*(num_spectrum_files -1)

    for k in range(1,num_spectrum_files):
      ratio = []
      unc = []
      for i in range(0, len(forward_normalized_mean[k])):
          ratio.append( forward_normalized_mean[k][i]/forward_normalized_mean[0][i] )

          sigma_f = forward_error[k][i]
          sigma_t = forward_error[0][i]

          f_squared = forward_normalized_mean[k][i]*forward_normalized_mean[k][i]
          t_squared = forward_normalized_mean[0][i]*forward_normalized_mean[0][i]

          unc.append( math.sqrt( sigma_f*sigma_f + (f_squared/t_squared)*sigma_t*sigma_t )/forward_normalized_mean[0][i] )
      f_over_t[k-1] = ratio
      f_over_t_unc[k-1] = unc

    # Set up the top subplot

    # Plot experimental data
    directory = os.path.dirname(os.path.abspath(__file__))

    for i in range(len(exp_names)):
      filename = directory + "/Al/experimental_results/" + exp_names[i] +".tsv"
      with open(filename) as input:
          name = input.readline().strip()
          input.readline()
          data = zip(*(line.strip().split('\t') for line in input))
          x = [None] * len(data[0][:])
          x = [0 for k in range(len(data[0][:]))]
          y = [0 for k in range(len(data[1][:]))]
          for j in range(len(x)):
            x[j] = float(data[0][j])*1e-3
            y[j] = float(data[1][j])

      if i == 1:
        ax[0].scatter(x, y, label="Experimental", marker=markers[1], s=50, facecolors='none', edgecolors='b' )
      else:
        ax[0].scatter(x, y, marker=markers[1], s=50, facecolors='none', edgecolors='b' )

    for k in range(0,num_spectrum_files):
      # Plot forward histogram of results
      label = forward_spectrum_data_name + "_" + str(k)
      if not FORWARD_NORM == 1.0:
        label += "*" + str(FORWARD_NORM)

      m_bin, bins, plt1 = ax[0].hist(forward_energy_bins[k][:-1], bins=forward_energy_bins[k], weights=forward_normalized_mean[k], histtype='step', label=label, color=marker_color[k], linestyle=linestyles[k], linewidth=1.8 )

      # Plot error bars
      forward_mid = 0.5*(bins[1:] + bins[:-1])
      # plt2 = ax[0].errorbar(forward_mid, m_bin, yerr=forward_error, ecolor='b', fmt=None)

      ax[0].set_xscale("log")

      y_label = data_type + " Spectrum"

      ax[0].set_ylabel( y_label )

      if not legend_pos is None:
          ax[0].legend(frameon=True, bbox_to_anchor=legend_pos)
      else:
          ax[0].legend(frameon=True)

      # Turn on the grid
      ax[0].grid(True, linestyle=':', linewidth=1)

      # Set the x limits
      if not xlims is None:
          ax[0].set_xlim( xlims[0], xlims[-1] )
      else:
          ax[0].set_xlim( forward_energy_bins[k][0], forward_energy_bins[k][-1] )

      if not top_ylims is None:
          ax[0].set_ylim( top_ylims[0], top_ylims[1] )

      # Set the y tic labels
      yticklabels = ax[0].yaxis.get_ticklabels()
      # yticklabels[0].set_visible(False)
      yticklabels[-1].set_visible(False)

      # Set the tic properties
      ax[0].yaxis.set_ticks_position("both")
      ax[0].xaxis.set_ticks_position("both")
      ax[0].tick_params(direction="in", width=edge_thickness)
      ax[0].tick_params(which="minor", direction="in", width=edge_thickness)

      for axis in ['top','bottom','left','right']:
          ax[0].spines[axis].set_linewidth(edge_thickness)

    for k in range(1,num_spectrum_files):
      # Set up the bottom subplot
      ax[1].set_xscale("log")

      label = "Forward_" + str(k) + "/Forward_0"
      ax[1].errorbar( forward_mid, f_over_t[k-1], yerr=f_over_t_unc[k-1], capsize=1.5, fmt='o', ecolor=marker_color[k], color=marker_color[k], linewidth=0.5, markersize=1.9, label=label )

      ax[1].set_ylabel( "F/R" )
      ax[1].set_xlabel( "Energy (MeV)" )

      if not legend_pos is None:
          ax[1].legend(frameon=True, bbox_to_anchor=legend_pos)
      else:
          ax[1].legend(frameon=True)

      # Turn on the grid
      ax[1].grid(True, linestyle=':', linewidth=1)

      # Set the x limits
      if not xlims is None:
          ax[1].set_xlim( xlims[0], xlims[-1] )
      else:
          ax[1].set_xlim( forward_energy_bins[k][0], forward_energy_bins[k][-1] )

      if not bottom_ylims is None:
          ax[1].set_ylim( bottom_ylims[0], bottom_ylims[1] )

      # Set the y tic labels
      yticklabels = ax[1].yaxis.get_ticklabels()
      yticklabels[0].set_visible(False)
      yticklabels[-1].set_visible(False)

      # Set the tic properties
      ax[1].yaxis.set_ticks_position("both")
      ax[1].xaxis.set_ticks_position("both")
      ax[1].tick_params(direction="in", width=edge_thickness)
      ax[1].tick_params(which="minor", direction="in", width=edge_thickness)

      for axis in ['top','bottom','left','right']:
          ax[1].spines[axis].set_linewidth(edge_thickness)

      # Save the figure
      for i in range(0,len(output_plot_names)):
        fig.savefig( output_plot_names[i] )

    plt.show()