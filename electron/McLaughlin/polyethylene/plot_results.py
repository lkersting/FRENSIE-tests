#! /usr/bin/env python
# Luke Kersting
# This script asks for energy deposition data and run names which it then plots.
import csv
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D
import argparse as ap
import inspect, os

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

# Number of files
N = len(file_paths)

# Number of data points in each file
M = 50
data_x = [[] for y in range(N)]
data_y = [[] for y in range(N)]
data_error = [[] for y in range(N)]
names = [0 for x in range(N)]

# Get computational results
for n in range(N):

    with open(file_paths[n]) as input:
        names[n] = input.readline()[1:].strip()
        print names[n]
        print input.readline().strip()[1:]
        data = zip(*(line.strip().split('\t') for line in input))
        data_x[n] = np.asfarray(data[0][0:M])*1000.0
        data_y[n] = np.asfarray(data[1][0:M])
        data_error[n]= np.asfarray(data[2][0:M])*data_y[n]

# Plot
if user_args.m:
  fig = plt.figure(num=1, figsize=(9,7))
  # set height ratios for sublots
  gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
else:
  fig = plt.figure(num=1, figsize=(10,5))
  # set height ratios for sublots
  gs = gridspec.GridSpec(1, 1)

# fig = plt.figure(num=1, figsize=(10,5))

# the first subplot
ax0 = plt.subplot(gs[0])

x_label = 'Range ($\mathrm{mg/cm^2}$)'
plt.xlabel(x_label, size=14)
plt.ylabel('Dose ($\mathrm{MeV\/cm^2/g}$)', size=14)
plt.title('$\mathrm{Energy\/Deposition\/from\/100\/keV\/Electron\/in\/Polyethylene}$', size=16)
ax=plt.gca()

plt.xlim(0.0,1.2)
plt.ylim(0.0,3.5)
# if user_args.m:
#     plt.ylim(0.0,0.05)

plots = []
labels = []
if user_args.e:
    directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

    # Get experimental data
    filename = directory + "/experimental_results.tsv"
    with open(filename) as input:
        data = zip(*(line.strip().split('\t') for line in input))
        exp_x = np.asfarray(data[0][1:])
        exp_y = np.asfarray(data[1][1:])

    # Plot the experimental data
    # plt1, = plt.scatter(exp_x, exp_y, label="McLaughlin (Exp.)", marker='s' )
    plt1, = plt.plot(exp_x, exp_y, linestyle='None', marker='s', markersize=5, color='b' )

    plots.append( plt1 )
    labels.append("McLaughlin (Exp.)" )

markers = ["--v",":^","--<","-.>",":+","--x","-.1",":2","--3","-.4",":8","--p","-.P",":*","--h","-.H",":X","--D","-.d"]
linestyle = ["--","-.",":","--","-.",":","--","-.",":","--3","-.",":","--","-.",":","--","-.",":","--","-."]
markerssizes = [6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6]
marker_color = ['g', 'r', 'm', 'k', 'y', 'c', 'g', 'r', 'm', 'k', 'y', 'c']

# linestyles: 'solid', 'dashed', 'dashdotted', 'densely dotted', 'dashdotdotted', 'densely dashed', 'densely dashdotted', 'densely dashdotdotted', 'dotted', 'loosely dashed', 'loosely dashdotted', 'loosely dashdotdotted')

linestyles = [(0, ()), (0, (5, 5)), (0, (3, 5, 1, 5)), (0, (1, 1)), (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (1, 5)), (0, (5, 10)), (0, (3, 10, 1, 10)), (0, (3, 10, 1, 10, 1, 10))]

if user_args.m:
    names = ['MCNP6.2','FRENSIE-ACE', 'FRENSIE-ENDL' ]
for n in range(N):
    x = np.insert( data_x[n], 0, 0.0)

    # Plot histogram of results
    m, bins, plt1 = plt.hist(x[:-1], bins=x, weights=data_y[n], histtype='step', label=names[n], color=marker_color[n], linestyle=linestyles[n], linewidth=1.8 )

    # Plot error bars
    mid = 0.5*(bins[1:] + bins[:-1])
    plt2 = plt.errorbar(mid, m, yerr=data_error[n], ecolor=marker_color[n], fmt=None)

    handle1 = Line2D([], [], c=marker_color[n], linestyle='--', dashes=linestyles[n][1], linewidth=1.8)
    plots.append( handle1 )
    labels.append(names[n])

# plt.legend(loc=1)
plt.legend(plots, labels, loc=1)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))

output = "mclaughlin_polyethylene_results.eps"
if user_args.o:
    output = user_args.o

ax.grid(linestyle=':')

print "Plot outputted to: ",output
fig.savefig(output, bbox_inches='tight', dpi=600)
plt.show()