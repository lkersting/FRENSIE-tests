import os
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})

class bcolors:
    HEADER = '\033[95m'
    SIGMA2 = '\033[94m'
    SIGMA1 = '\033[92m'
    SIGMA3 = '\033[93m'
    NO_SIGMA = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def plotInfiniteMediumSimulationSurfaceFlux( forward_data,
                                             adjoint_data,
                                             energy_bins,
                                             output_name,
                                             radius,
                                             element,
                                             top_ylims = None,
                                             bottom_ylims = None,
                                             xlims = None,
                                             legend_pos = None ):

  # Make sure that the forward data is a dictionary
  if not isinstance(forward_data, dict):
      print "The forward data must be a dictionary with keys \"e_bins\", \"mean\" and \"re\""
      sys.exit(1)

  # Make sure that the tadjoint data is a dictionary
  if not isinstance(adjoint_data, dict):
      print "The adjoint data must be a dictionary with keys \"mean\" and \"re\""
      sys.exit(1)

  # Get the forward energy bin boundaries
  energy_bins = np.array(energy_bins)
  # Get the x bin widths
  bin_widths = (energy_bins[1:] - energy_bins[:-1])

  # Get forward data
  # Average the flux to the bin width
  forward_y = np.array(forward_data['mean'])/bin_widths
  # Calculate the error for the bin averaged surface flux
  forward_error = np.array(forward_data['re'])*forward_y

  # Get Adjoint Data
  NORM = 1.0

  # Average the flux to the bin width
  adjoint_y = np.array(adjoint_data['mean'])/bin_widths*NORM
  # Calculate the error for the bin averaged surface flux
  adjoint_error = np.array(adjoint_data['re'])*adjoint_y

  # Plot
  fig = plt.figure(figsize=(10,6))

  # set height ratios for sublots
  gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])

  # the first subplot
  ax0 = plt.subplot(gs[0])

  if element == "Pb":
    atom_name = "Lead"
  elif element == "H":
    atom_name = "Hydrogen"

  mfps = ('%f' % (radius/2.0)).rstrip('0').rstrip('.')

  plot_title = r'\textbf{Surface Flux in an Infinite Medium of ' + element + ' at ' + str(mfps) +' mfps}'
  x_label = r'\textbf{Energy (MeV)}'
  plt.xlabel(x_label)
  plt.ylabel(r'\textbf{Surface Flux}')
  plt.title( plot_title, size=18)
  ax=plt.gca()

  if not top_ylims is None:
    plt.ylim(top_ylims[0],top_ylims[1])

  markers = ["--v","-.o",":^","--<","-.>",":+","--x","-.1",":2","--3","-.4",":8","--p","-.P",":*","--h","-.H",":X","--D","-.d"]
  markerssizes = [6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6]
  marker_color = ['g', 'r', 'c', 'm', 'y', 'k', 'w', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

  linestyles = [(0, ()), (0, (5, 5)), (0, (3, 5, 1, 5)), (0, (1, 1)), (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (1, 5)), (0, (5, 10)), (0, (3, 10, 1, 10)), (0, (3, 10, 1, 10, 1, 10))]

  plots = []
  labels = []

  # Plot Adjoint Data

  # Plot histogram of results
  # label = "Adjoint~\,- No Atomic Excitation"
  label = "Adjoint"
  if not NORM == 1.0:
    label += "*" + str(NORM)
  m, bins, plt1 = plt.hist(energy_bins[:-1], bins=energy_bins, weights=adjoint_y, histtype='step', label=r"\textbf{" + label + '}', color='b', linestyle=linestyles[0], linewidth=1.8 )

  # Plot error bars
  mid = 0.5*(bins[1:] + bins[:-1])
  # plt2 = plt.errorbar(mid, m, yerr=adjoint_error, ecolor='b', fmt=None)

  # plt.errorbar(mid, adjoint_y, yerr=adjoint_error, label="adjoint", fmt="--s", markersize=6, color='b' )

  handle1 = Line2D([], [], c='b', linestyle='--', dashes=linestyles[0][1], linewidth=1.8)
  plots.append( handle1 )
  labels.append(label)

  # Plot Forward Data

  # Plot histogram of results
  # label = "Forward - No Atomic Excitation"
  label = "Forward"
  m, bins, plt1 = plt.hist(energy_bins[:-1], bins=energy_bins, weights=forward_y, histtype='step', label=r"\textbf{" + label + '}', color='g', linestyle=linestyles[1], linewidth=1.8 )

  # Plot error bars
  mid = 0.5*(bins[1:] + bins[:-1])
  # plt2 = plt.errorbar(mid, m, yerr=forward_error, ecolor='g', fmt=None)

  # plt.errorbar(mid, forward_y, yerr=forward_error, label="forward", fmt="--s", markersize=6, color='b' )

  handle1 = Line2D([], [], c='g', linestyle='--', dashes=linestyles[1][1], linewidth=1.8)
  plots.append( handle1 )
  labels.append(label)


  plt.legend(loc=legend_pos, fontsize=14)
  # ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

  markers = ["v","o","^","<",">","+","x","1","2","3","4","8","p","P","*","h","H","X","D","d"]

  # The C/R subplot (with shared x-axis)
  ax1 = plt.subplot(gs[1], sharex = ax0)
  plt.xlabel(x_label)
  plt.ylabel(r'\textbf{Adjoint/Forward}')

  yerr = np.sqrt( ((1.0/forward_y)**2)*(adjoint_error)**2 + ((energy_bins[:-1]/forward_y**2)**2)*(forward_error)**2 )
  y = adjoint_y/forward_y

  if output_name is None:
    output_name = element + "_infinite_medium_" + str(radius)

  output_data_name = output_name + "_3_sigma.txt"

  f = open(output_data_name, 'w')
  f.write( "\n#Energy\tRatio\tUncertainty\n" )

  # calculate % of C/R values within 1,2,3 sigma
  num_in_one_sigma = 0
  num_in_two_sigma = 0
  num_in_three_sigma = 0
  num_below = 0
  num_above = 0

  N=0
  max_k = 0
  length = len(y)-N
  for i in range(N, len(y)):
    # Print C/R results
    # print energy_bins[i+1], ": ", (1.0-y[i])*100, u"\u00B1", yerr[i]*100, "%"
    # print energy_bins[i+1], ": ", y[i], "\t",forward_y[i]
    if not np.isfinite( y[i] ):
      # print energy_binsx[i+1], ": ", y[i], "\t",forward_y[i], "\t",adjoint_y[i]
      if forward_y[i] == adjoint_y[i]:
        y[i] = 1.0
        yerr[i] = 0.0
      else:
        y[i] = 0
        yerr[i] = 1
    else:
      max_k = i

    # Calculate number above and below reference
    if y[i] < 1.0:
      num_below += 1
    else:
      num_above += 1

    diff = abs( 1.0 - y[i] )
    # message = '%.4e' % energy_bins[i+1] + ": " + '%.6f' % (y[i]) + u"\u00B1" + '%.6f' % (yerr[i]) + "%"
    message = '%.4e' % energy_bins[i+1] + "\t" + '%.6f' % (y[i]) +"\t"+ '%.6f' % (yerr[i])

    sigma = bcolors.NO_SIGMA
    if diff <= 3*yerr[i]:
        num_in_three_sigma += 1
        sigma = bcolors.SIGMA3
    if diff <= 2*yerr[i]:
        num_in_two_sigma += 1
        sigma = bcolors.SIGMA2
    if diff <= yerr[i]:
        num_in_one_sigma += 1
        sigma = bcolors.SIGMA1

    f.write( message +"\n")
    message = sigma + message + bcolors.ENDC
    # print message

  message = "----------------------------------------------------------------"
  print message
  f.write( message +"\n")
  message = '%.3f' % (float(num_above)/length*100) + "% above reference"
  print "  ", message
  f.write( message +"\n")
  message = '%.3f' % (float(num_below)/length*100) + "% below reference"
  print "  ", message
  f.write( message +"\n")
  message = "----------------------------------------------------------------"
  print message
  f.write( message +"\n")
  message = '%.3f' % (float(num_in_one_sigma)/length*100) + "% C/R within 1 sigma"
  print "  ", bcolors.SIGMA1, message, bcolors.ENDC
  f.write( message +"\n")
  message = "----------------------------------------------------------------"
  print message
  f.write( message +"\n")
  message = '%.3f' % (float(num_in_two_sigma)/length*100) + "% C/R within 2 sigma"
  print "  ", bcolors.SIGMA2, message, bcolors.ENDC
  f.write( message +"\n")
  message = "----------------------------------------------------------------"
  print message
  f.write( message +"\n")
  message = '%.3f' % (float(num_in_three_sigma)/length*100) + "% C/R within 3 sigma"
  print "  ", bcolors.SIGMA3, message, bcolors.ENDC
  f.write( message +"\n")
  f.close()

  max_k += 1

  # Plot histogram of results
  m, bins, _ = ax1.hist(energy_bins[:max_k], bins=energy_bins[:max_k+1], weights=y[:max_k], histtype='step', label=r"\textbf{ratio}", color='b', linestyle=linestyles[0], linewidth=1.8 )
  # Plot error bars
  mid = 0.5*(bins[1:] + bins[:-1])
  ax1.errorbar(mid, m, yerr=yerr[:max_k], ecolor='b', fmt=None)

  # make x ticks for first suplot invisible
  plt.setp(ax0.get_xticklabels(), visible=False)

  # remove first tick label for the first subplot
  yticks = ax0.yaxis.get_major_ticks()
  yticks[0].label1.set_visible(False)
  ax0.grid(linestyle=':')
  ax1.grid(linestyle=':')

  output_plot_names = []
  output_plot_names.append( output_name + ".eps" )
  output_plot_names.append( output_name + ".png" )
  output_plot_names.append( output_name + ".pdf" )

  if not xlims is None:
    plt.xlim(xlims[0],xlims[1])
  if not bottom_ylims is None:
    plt.ylim(bottom_ylims[0],bottom_ylims[1])

  # remove vertical gap between subplots
  plt.subplots_adjust(hspace=.0)

  print "Plot outputted to: ",output_name
  # Save the figure
  for i in range(0,len(output_plot_names)):
    fig.savefig( output_plot_names[i], bbox_inches='tight', dpi=600)
  plt.draw()
  plt.show()


def plotAllInfiniteMediumSimulationSurfaceFlux( forward_data,
                                                adjoint_data,
                                                energy_bins,
                                                output_name,
                                                radius,
                                                top_ylims = None,
                                                bottom_ylims = None,
                                                xlim = None,
                                                legend_pos = None ):

  linestyles = [(0, ()), (0, (5, 5)), (0, (3, 5, 1, 5)), (0, (1, 1)), (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (1, 5)), (0, (5, 10)), (0, (3, 10, 1, 10)), (0, (3, 10, 1, 10, 1, 10))]

  plots = []
  labels = []

  num = len(forward_data)

  # Plot
  fig = plt.figure(figsize=(10,4*num))

  gs = [[] for y in range(num)]
  axes = [[] for y in range(2*num)]

  ratios = [2,1]*num
  y_labels = [r'\textbf{Surface Flux}', r'\textbf{C/R}']*num

  for i in range(num):
    # We'll use two separate gridspecs to have different margins, hspace, etc
    top = 0.95 - i/500.0
    bottom = 0.1 - i/500.0
    ratios = [2,1]*num
    gs[i] = gridspec.GridSpec(2*num, 1, height_ratios=ratios, top=top, bottom=bottom, hspace=0)

  # Set the first plot
  axes[0] = fig.add_subplot(gs[0][0,:])
  axes[0].set_ylabel(y_labels[0], size=16)
  axes[0].grid(linestyle=':')
  # axes[0].set_title('Results for a 0.01 MeV Point Source in a H Sphere', size=18)

  for i in range(1,len(axes)):
    j = i/2
    # Shared axes with C/R
    axes[i] = fig.add_subplot(gs[j][i,:], sharex=axes[i-1])
    # Hide shared x-tick labels
    plt.setp(axes[i-1].get_xticklabels(), visible=False)
    # Set the y labels
    axes[i].set_ylabel(y_labels[i], size=16)

    # remove first tick label for the first subplot
    yticks = axes[i].yaxis.get_major_ticks()
    # yticks[0].label1.set_visible(False)
    axes[i].grid(linestyle=':')

  plot_titles = [ '(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

  # Set the x label
  x_label = r'\textbf{Energy (MeV)}'
  axes[5].set_xlabel(x_label, size=16)

  ax=plt.gca()
  ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

  # Get the forward energy bin boundaries
  energy_bins = np.array(energy_bins)
  # Get the x bin widths
  bin_widths = (energy_bins[1:] - energy_bins[:-1])

  for i in range(num):
    j = i*2

    # Make sure that the forward data is a dictionary
    if not isinstance(forward_data[i], dict):
        print "The forward data must be a dictionary with keys \"e_bins\", \"mean\" and \"re\""
        sys.exit(1)

    # Make sure that the tadjoint data is a dictionary
    if not isinstance(adjoint_data[i], dict):
        print "The adjoint data must be a dictionary with keys \"mean\" and \"re\""
        sys.exit(1)

    # Get forward data
    # Average the flux to the bin width
    forward_y = np.array(forward_data[i]['mean'])/bin_widths
    # Calculate the error for the bin averaged surface flux
    forward_error = np.array(forward_data[i]['re'])*forward_y

    # Get Adjoint Data
    NORM = 1.0

    # Average the flux to the bin width
    adjoint_y = np.array(adjoint_data[i]['mean'])/bin_widths*NORM
    # Calculate the error for the bin averaged surface flux
    adjoint_error = np.array(adjoint_data[i]['re'])*adjoint_y


    # the first subplot
    mfps = ('%f' % (radius[i]/2.0)).rstrip('0').rstrip('.')

    # place a text box in upper left in axes coords
    axes[j].text(0.5, 0.95,  r'\textbf{'+ plot_titles[i] + '}', transform=axes[j].transAxes, horizontalalignment='center', verticalalignment='top', fontsize='18' )
    axes[j].text(0.5, 0.80, r'\textbf{'+ mfps +' mfps}', transform=axes[j].transAxes, horizontalalignment='center', verticalalignment='top', fontsize='18' )

    if not top_ylims is None:
      axes[j].set_ylim(top_ylims[i][0],top_ylims[i][1])

    # Plot Adjoint Data

    # Plot histogram of results
    label = r"\textbf{Adjoint}"
    if not NORM == 1.0:
      label += "*" + str(NORM)
    m, bins, plt1 = axes[j].hist(energy_bins[:-1], bins=energy_bins, weights=adjoint_y, histtype='step', label=label, color='b', linestyle=linestyles[0], linewidth=1.8 )

    # Plot error bars
    mid = 0.5*(bins[1:] + bins[:-1])
    # plt2 = plt.errorbar(mid, m, yerr=adjoint_error, ecolor='b', fmt=None)

    # plt.errorbar(mid, adjoint_y, yerr=adjoint_error, label="adjoint", fmt="--s", markersize=6, color='b' )

    handle1 = Line2D([], [], c='b', linestyle='--', dashes=linestyles[0][1], linewidth=1.8)
    plots.append( handle1 )
    labels.append(label)

    # Plot Forward Data

    # Plot histogram of results
    label = r"\textbf{Forward}"
    m, bins, plt1 = axes[j].hist(energy_bins[:-1], bins=energy_bins, weights=forward_y, histtype='step', label=label, color='g', linestyle=linestyles[1], linewidth=1.8 )

    # Plot error bars
    mid = 0.5*(bins[1:] + bins[:-1])
    # plt2 = plt.errorbar(mid, m, yerr=forward_error, ecolor='g', fmt=None)

    # plt.errorbar(mid, forward_y, yerr=forward_error, label="forward", fmt="--s", markersize=6, color='b' )

    handle1 = Line2D([], [], c='g', linestyle='--', dashes=linestyles[1][1], linewidth=1.8)
    plots.append( handle1 )
    labels.append(label)

    # The C/R subplot (with shared x-axis)

    yerr = np.sqrt( ((1.0/forward_y)**2)*(adjoint_error)**2 + ((energy_bins[:-1]/forward_y**2)**2)*(forward_error)**2 )
    y = adjoint_y/forward_y

    if output_name is None:
      output_name = element + "_infinite_medium_" + str(radius)

    output_data_name = output_name + "_3_sigma.txt"

    f = open(output_data_name, 'w')
    f.write( "\n#Energy\tRatio\tUncertainty\n" )

    # calculate % of C/R values within 1,2,3 sigma
    sigma_bins = [1e-4, 2e-3, 1e-2]
    sigma_bins = [1e-4, 1e-2]
    sigma_k = [0]*len(sigma_bins)

    num_in_one_sigma = [0]*(len(sigma_bins)-1)
    num_in_two_sigma = [0]*(len(sigma_bins)-1)
    num_in_three_sigma = [0]*(len(sigma_bins)-1)
    num_below = [0]*(len(sigma_bins)-1)
    num_above = [0]*(len(sigma_bins)-1)
    num_equal = [0]*(len(sigma_bins)-1)

    for l in range(len(sigma_bins)-1):
      for k in range(len(energy_bins)):
        if energy_bins[k] > sigma_bins[l]:
          sigma_k[l] = k-1
          break

    sigma_k[-1] = len(y)-1
    for k in range(len(energy_bins)):
      if energy_bins[k] > sigma_bins[-1]:
        sigma_k[-1] = k
        break

    length = [0]*(len(sigma_bins)-1)

    max_k = 0
    for k in range(0, len(y)):
      # Print C/R results
      # print energy_bins[k+1], ": ", (1.0-y[k])*100, u"\u00B1", yerr[k]*100, "%"
      # print energy_bins[k+1], ": ", y[k], "\t",forward_y[k]

      m = 0
      for l in range(len(sigma_bins)-1):
        if k <= sigma_k[l+1] and k > sigma_k[l]:
          m = l

      if not np.isfinite( y[k] ):
        # print energy_binsx[k+1], ": ", y[k], "\t",forward_y[k], "\t",adjoint_y[k]
        if forward_y[k] == adjoint_y[k]:
          y[k] = 1.0
          yerr[k] = 0.0
        else:
          y[k] = 0
          yerr[k] = 1
      else:
        max_k = k

        # Calculate number above and below reference
        if y[k] < 1.0:
          num_below[m] += 1
        elif y[k] > 1.0:
          num_above[m] += 1
        else:
          num_equal[m] += 1
        length[m] += 1

        diff = abs( 1.0 - y[k] )
        # message = '%.4e' % energy_bins[k+1] + ": " + '%.6f' % (y[k]) + u"\u00B1" + '%.6f' % (yerr[k]) + "%"
        message = '%.4e' % energy_bins[k+1] + "\t" + '%.6f' % (y[k]) +"\t"+ '%.6f' % (yerr[k])

        sigma = bcolors.NO_SIGMA
        if diff <= 3*yerr[k]:
            num_in_three_sigma[m] += 1
            sigma = bcolors.SIGMA3
        if diff <= 2*yerr[k]:
            num_in_two_sigma[m] += 1
            sigma = bcolors.SIGMA2
        if diff <= yerr[k]:
            num_in_one_sigma[m] += 1
            sigma = bcolors.SIGMA1

        f.write( message +"\n")
        message = sigma + str(m) + ' ' + message + bcolors.ENDC
        print message

    min_energy = str(energy_bins[sigma_k[0]])
    for k in range(len(sigma_bins)-1):
      max_energy = str(energy_bins[sigma_k[k+1]+1])
      message = "\n----------------------------------------------------------------"
      print message
      max_energy = str(energy_bins[sigma_k[k+1]+1])
      message = "For energy in range " + min_energy + " to " + max_energy
      min_energy = max_energy
      print message
      message = "----------------------------------------------------------------"
      print message
      f.write( message +"\n")
      message = '%.3f' % (float(num_above[k])/length[k]*100) + "% above reference"
      print "  ", message
      f.write( message +"\n")
      message = '%.3f' % (float(num_below[k])/length[k]*100) + "% below reference"
      print "  ", message
      f.write( message +"\n")
      message = '%.3f' % (float(num_equal[k])/length[k]*100) + "% equal to reference"
      print "  ", message
      f.write( message +"\n")
      message = "----------------------------------------------------------------"
      print message
      f.write( message +"\n")
      message = '%.3f' % (float(num_in_one_sigma[k])/length[k]*100) + "% C/R within 1 sigma"
      print "  ", bcolors.SIGMA1, message, bcolors.ENDC
      f.write( message +"\n")
      message = "----------------------------------------------------------------"
      print message
      f.write( message +"\n")
      message = '%.3f' % (float(num_in_two_sigma[k])/length[k]*100) + "% C/R within 2 sigma"
      print "  ", bcolors.SIGMA2, message, bcolors.ENDC
      f.write( message +"\n")
      message = "----------------------------------------------------------------"
      print message
      f.write( message +"\n")
      message = '%.3f' % (float(num_in_three_sigma[k])/length[k]*100) + "% C/R within 3 sigma"
      print "  ", bcolors.SIGMA3, message, bcolors.ENDC
      f.write( message +"\n")
    f.close()

    max_k += 1

    # Plot histogram of results
    m, bins, _ = axes[j+1].hist(energy_bins[:max_k], bins=energy_bins[:max_k+1], weights=y[:max_k], histtype='step', label=r"\textbf{ratio}", color='b', linestyle=linestyles[0], linewidth=1.8 )
    # Plot error bars
    mid = 0.5*(bins[1:] + bins[:-1])

    axes[j+1].errorbar(mid, m, yerr=yerr[:max_k], ecolor='b', fmt=None)

    output_plot_names = []
    output_plot_names.append( output_name + ".eps" )
    output_plot_names.append( output_name + ".png" )
    output_plot_names.append( output_name + ".pdf" )

    # remove first tick label for the first subplot
    yticks = axes[j].yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)

    if not bottom_ylims is None:
      axes[j+1].set_ylim(bottom_ylims[i][0],bottom_ylims[i][1])

    axes[j].legend(loc=legend_pos[i], fontsize=14)

  if not xlim is None:
    plt.xlim(xlim[0],xlim[1])

  print "Plot outputted to: ",output_name
  # Save the figure
  for i in range(0,len(output_plot_names)):
    fig.savefig( output_plot_names[i], bbox_inches='tight', dpi=600)

  plt.show()