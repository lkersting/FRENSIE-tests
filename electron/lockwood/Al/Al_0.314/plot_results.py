#! /usr/bin/env python
# Luke Kersting
# This script asks for energy deposition data and run names which it then plots.
import csv
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
from pylab import figure, title, setp, close, clf
import argparse as ap
import inspect, os
import scipy.integrate as integrate

# Set up the argument parser
description = "This script asks for energy deposition data and run names which "\
              "which it then plots."

parser = ap.ArgumentParser(description=description)

experimental_msg = "Flag to add the experimental data to the generated plot."
parser.add_argument('-e', help=experimental_msg, action='store_true')

mcnp_comparison_msg = "Flag to compare against mcnp data run."
parser.add_argument('-m', help=mcnp_comparison_msg, action='store_true')

output_msg = "The output file name."
parser.add_argument('-o', help=output_msg, required=False)

parser.add_argument("input_files", nargs='*')

# Parse the user's arguments
user_args = parser.parse_args()
file_paths = user_args.input_files

# Set the polynomial fit order
poly = 5

# Set the y min and max value for dose
ymin = 0.0 # Always 0.0
ymax = 7 # 15.0 (lin-lin), 10.0 (lin-log), 7.0 (log-log)
# Set the y min and max value for C/E
ratio_min = 0.0 # 0.0 (lin-lin), 0.0 (lin-log), 0.0 (log-log)
ratio_max = 2.0 # 4.5 (lin-lin), 2.5 (lin-log), 2.0 (log-log)

# C/E ranges
ranges=[0.0, 0.6, 1.0]

# provide the x data in fraction of a mean range units
FMR =[0.0224, 0.0831, 0.1600, 0.2260, 0.2970, 0.3570, 0.4220, 0.5010, 0.5790, 0.6380, 0.7170, 0.8790]
CSDA_R = 1.13E-01

# Number of files
N = len(file_paths)

data_x = [[] for y in range(N)]
data_y = [[] for y in range(N)]
data_error = [[] for y in range(N)]
policy = np.empty(N, dtype=object)
interp = np.empty(N, dtype=object)
elastic = np.empty(N, dtype=object)
names = np.empty(N, dtype=object)
fit = np.empty(N, dtype=object)
integral = np.zeros(N)

# Get computational results
for n in range(N):

    with open(file_paths[n]) as input:
        fullname = input.readline().strip()[1:]
        interp[n] = fullname[0:7]
        policy[n] = fullname[8:-3]
        elastic[n] = fullname[-2:]

        data = zip(*(line.strip().split('\t') for line in input))
        data_name = data[0][0] + data[1][0] + data[2][0]
        data_x[n] = np.asfarray(data[0][1:])
        data_y[n] = np.asfarray(data[1][1:])
        # Error given as relative error
        data_error[n] = np.asfarray(data[2][1:])*data_y[n]

# Set the plot legend names and legend title
title = ''
if len(np.unique(interp)) == 1:
  title = np.unique(interp)[0]
else:
  names = interp
if len(np.unique(policy)) == 1:
  if len(title) > 0:
    title += ' ' + np.unique(policy)[0]
  else:
    title = np.unique(policy)[0]
else:
  if names[0] != None:
    names += ' ' + policy
  else:
    names = policy
if len(np.unique(elastic)) == 1:
  if len(title) > 0:
    title += ' '
  else:
    title = ''

  if np.unique(elastic)[0] == 'CE':
    title += ' Coupled Runs'
  elif np.unique(elastic)[0] == 'DE':
    title += ' Decoupled Runs'
  elif np.unique(elastic)[0] == 'HE':
    title += ' Hybrid Runs'
else:
  if names[0] != None:
    for i in range(len(elastic)):
      if elastic[i] == 'CE':
        names[i] += ' Coupled'
      elif elastic[i] == 'DE':
        names[i] += ' Decoupled'
      elif elastic[i] == 'HE':
        names[i] += ' Hybrid'
      else:
        names[i] += ' ' + elastic[i]
  else:
    for i in range(len(elastic)):
      if elastic[i] == 'CE':
        names[i] = 'Coupled'
      elif elastic[i] == 'DE':
        names[i] = 'Decoupled'
      elif elastic[i] == 'HE':
        names[i] = 'Hybrid'
      else:
        names[i] = elastic[i]


# Plot
fig = plt.figure(num=1, figsize=(10,6))

# set height ratios for sublots
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])

# the first subplot
ax0 = plt.subplot(gs[0])

# x_label = 'Range ($\mathrm{g/cm^2}$)'
x_label = 'Fraction of a Mean Range'

plt.xlabel(x_label, size=14)
plt.ylabel('Energy Deposition ($\mathrm{MeV\/cm^2/g}$)', size=14)
plt.title('$\mathrm{Energy\/Deposition\/from\/0.314\/MeV\/Electron\/in\/Aluminum}$', size=16)
ax=plt.gca()


plt.xlim(0.0,0.9)
plt.ylim(ymin,ymax)

if user_args.e:
    # Get experimental data
    directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    filename = directory + "/experimental_results.txt"
    with open(filename) as input:
        data = zip(*(line.strip().split('\t') for line in input))
        data_name = data[0][0] + data[1][0] + data[2][0]
        exp_x = np.asfarray(data[0][1:])
        exp_y = np.asfarray(data[1][1:])
        # Error is given in % Rel Error
        exp_error = np.asfarray(data[2][1:])*exp_y/100.0

    print '\nLockwood Experimental'
    exp_fit = np.poly1d(np.polyfit(FMR, exp_y, poly))

    P = np.polyint(exp_fit)
    exp_int = P(FMR[-1])*CSDA_R
    print 'Integral using polyint      = ', exp_int

    x = np.linspace(0, FMR[-1], 1000)
    y = np.zeros(1000)
    for i in range(1000):
      y[i] = exp_fit(x[i])
    # plt.plot(x,y, color='b' )
    print 'Integral using simpson rule = ', integrate.simps(y, x)*CSDA_R

    # Plot the data
    line0,err0,arg3, = ax0.errorbar(FMR, exp_y, yerr=exp_error, label="Lockwood (Exp.)", fmt="-s", markersize=5, linewidth=1.8 )


markers = ["--v","-.o",":^","--<","-.>",":+","--x","-.1",":2","--3","-.4",":8","--p","-.P",":*","--h","-.H",":X","--D","-.d"]
linestyle = ["--","-.",":","--","-.",":","--","-.",":","--3","-.",":","--","-.",":","--","-.",":","--","-."]
markerssizes = [6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6]
marker_color = ['g', 'r', 'm', 'k', 'y', 'c', 'g', 'r', 'm', 'k', 'y', 'c']

# linestyles: 'solid', 'dashed', 'dashdotted', 'densely dotted', 'dashdotdotted', 'densely dashed', 'densely dashdotted', 'densely dashdotdotted', 'dotted', 'loosely dashed', 'loosely dashdotted', 'loosely dashdotdotted')

linestyles = [(0, (5, 5)), (0, (3, 5, 1, 5)), (0, (1, 1)), (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (1, 5)), (0, (5, 10)), (0, (3, 10, 1, 10)), (0, (3, 10, 1, 10, 1, 10))]

if user_args.m:
  names = ['MCNP6.2','FRENSIE-ACE', 'FRENSIE-ENDL' ]
for n in range(N):
    print '\n',names[n]
    # Plot the computational results
    plt.errorbar(FMR, data_y[n], yerr=data_error[n], label=names[n], fmt=markers[n], markersize=markerssizes[n], color=marker_color[n], linewidth=1.8)

    fit[n] = np.poly1d(np.polyfit(FMR, data_y[n], poly))

    P = np.polyint(fit[n])
    integral[n] = P(FMR[-1])*CSDA_R
    print 'Integral using polyint      = ', integral[n]

    x = np.linspace(0, FMR[-1], 1000)
    y = np.zeros(1000)
    for i in range(1000):
      y[i] = fit[n](x[i])
    # plt.plot(x,y, color=marker_color[n] )
    print 'Integral using simpson rule = ', integrate.simps(y, x)*CSDA_R
    if user_args.e:
      print 'C/E                         = ', integral[n]/exp_int
      print 'Relative error              =', (exp_int - integral[n])/exp_int*100, '%'

if len(title) > 0 and N < 6:
  l = plt.legend(loc='best', title=title)
  # Change the font size of the legend
  setp(l.get_texts(), fontsize=12)
  setp(l.get_title(), fontsize=15)
elif N > 6:
  lgd = plt.legend(loc="upper left", bbox_to_anchor=(1,1))
else:
  l = plt.legend(loc='best')

ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))


markers = ["v","o","^","<",">","+","x","1","2","3","4","8","p","P","*","h","H","X","D","d"]
if user_args.e:

    # The C/E subplot (with shared x-axis)
    ax1 = plt.subplot(gs[1], sharex = ax0)
    plt.xlabel(x_label, size=14)
    plt.ylabel('C/E', size=14)

    for n in range(N):
        # Calculate the propagated error for C/E
        yerr = np.sqrt( ((1.0/exp_y)**2)*(data_error[n])**2 + ((data_y[n]/exp_y**2)**2)*(exp_error)**2 )
        # Calculate C/E
        y = data_y[n]/exp_y

        print names[n]
        # Print C/R results

        max_diff = [0 for x in range(len(ranges)-1)]
        diff_range = [0 for x in range(len(ranges)-1)]
        diff_error = [0 for x in range(len(ranges)-1)]
        sum_diff = [0 for x in range(len(ranges)-1)]
        sum_err = [0 for x in range(len(ranges)-1)]
        range_len = [0 for x in range(len(ranges)-1)]
        r_i = 0
        sum_tot_diff = np.average( np.abs(1.0 - y))*100
        sum_tot_err = np.average( np.abs(yerr) )*100

        for i in range(0, len(y)):
          diff = (1.0-y[i])*100.0
          print FMR[i], ": ",diff, u"\u00B1", yerr[i]*100, "%"
          for j in range(1,len(ranges)-1):
            if FMR[i] > ranges[j] and FMR[i] < ranges[j+1]:
              r_i = j
              break

          sum_diff[r_i] += abs(diff)
          sum_err[r_i] += abs(yerr[i])*100.0
          range_len[r_i] += 1
          if abs(diff) > abs(max_diff[r_i]):
            max_diff[r_i] = diff
            diff_range[r_i] = FMR[i]
            diff_error[r_i] = yerr[i]*100.0

        print '--------------------------------------------------------'
        for j in range(len(max_diff)):
          print "--- range [",ranges[j] ,",",ranges[j+1] ,"] ---"
          print "Maximum percent relative diff:", max_diff[j], u"\u00B1", diff_error[j], "% at", diff_range[j]
          print "Average percent relative diff:", sum_diff[j]/range_len[j], u"\u00B1", sum_err[j]/range_len[j], "%"
        print '--------------------------------------------------------'
        print "Average percent relative diff: ", sum_tot_diff, u"\u00B1", sum_tot_err, "%"

        # Plot the C/E with errorbars
        ax1.errorbar(FMR, y, yerr=yerr, label=names[n], fmt=markers[n], markersize=markerssizes[n], color=marker_color[n])

if user_args.m:

    # The C/E subplot (with shared x-axis)
    ax1 = plt.subplot(gs[1], sharex = ax0)
    plt.xlabel(x_label, size=14)
    plt.ylabel('C/E', size=14)

    for n in range(1,N):
        # Calculate the propagated error for C/E
        yerr = np.sqrt( ((1.0/data_y[0])**2)*(data_error[n])**2 + ((data_y[n]/data_y[0]**2)**2)*(data_error[0])**2 )
        # Calculate C/E
        y = data_y[n]/data_y[0]

        print names[n]
        # Print C/R results

        max_diff = [0 for x in range(len(ranges)-1)]
        diff_range = [0 for x in range(len(ranges)-1)]
        diff_error = [0 for x in range(len(ranges)-1)]
        sum_diff = [0 for x in range(len(ranges)-1)]
        sum_err = [0 for x in range(len(ranges)-1)]
        range_len = [0 for x in range(len(ranges)-1)]
        r_i = 0
        sum_tot_diff = np.average( np.abs(1.0 - y))*100
        sum_tot_err = np.average( np.abs(yerr) )*100

        for i in range(0, len(y)):
          diff = (1.0-y[i])*100.0
          print FMR[i], ": ",diff, u"\u00B1", yerr[i]*100, "%"
          for j in range(1,len(ranges)-1):
            if FMR[i] > ranges[j] and FMR[i] < ranges[j+1]:
              r_i = j
              break

          sum_diff[r_i] += abs(diff)
          sum_err[r_i] += abs(yerr[i])*100.0
          range_len[r_i] += 1
          if abs(diff) > abs(max_diff[r_i]):
            max_diff[r_i] = diff
            diff_range[r_i] = FMR[i]
            diff_error[r_i] = yerr[i]*100.0

        print '--------------------------------------------------------'
        for j in range(len(max_diff)):
          print "--- range [",ranges[j] ,",",ranges[j+1] ,"] ---"
          print "Maximum percent relative diff:", max_diff[j], u"\u00B1", diff_error[j], "% at", diff_range[j]
          print "Average percent relative diff:", sum_diff[j]/range_len[j], u"\u00B1", sum_err[j]/range_len[j], "%"
        print '--------------------------------------------------------'
        print "Average percent relative diff: ", sum_tot_diff, u"\u00B1", sum_tot_err, "%"

        # Plot the C/E with errorbars
        ax1.errorbar(FMR, y, yerr=yerr, label=names[n], fmt=markers[n], markersize=markerssizes[n], color=marker_color[n])

# make x ticks for first suplot invisible
plt.setp(ax0.get_xticklabels(), visible=False)

# remove first tick label for the first subplot
yticks = ax0.yaxis.get_major_ticks()
yticks[0].label1.set_visible(False)
ax0.grid(linestyle=':')
ax1.grid(linestyle=':')

plt.ylim(ratio_min,ratio_max)

# remove vertical gap between subplots
plt.subplots_adjust(hspace=.0)

plt.show()

output = "lockwood_results.pdf"
if user_args.o:
    output = user_args.o

print "Plot outputted to: ",output
fig.savefig(output, bbox_inches='tight', dpi=600)