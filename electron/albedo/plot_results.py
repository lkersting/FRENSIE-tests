#! /usr/bin/env python
# Luke Kersting
# This script asks for albedo data and run names which it then plots.
import csv
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib import gridspec
import argparse as ap
import inspect, os
import pylab

# Set up the argument parser
description = "This script asks for albedo data and run names which "\
              "which it then plots against experimental data."

parser = ap.ArgumentParser(description=description)


mcnp_msg = "Flag to plot a comparison with mcnp (the first file should be for mcnp)."
parser.add_argument('-m', help=mcnp_msg, action='store_true')

all_experimental_msg = "Flag to add all experimental data (with source name) to the generated plot."
parser.add_argument('-a', help=all_experimental_msg, action='store_true')

all_experimental_msg = "Flag to add experimental data to the generated plot."
parser.add_argument('-e', help=all_experimental_msg, action='store_true')

output_msg = "The output file name."
parser.add_argument('-o', help=output_msg, required=False)

parser.add_argument("input_files", nargs='*')

# Parse the user's arguments
user_args = parser.parse_args()
file_paths = user_args.input_files

# Number of files
N = len(file_paths)

# Number of data points in each file
data_x = [0 for x in range(N)]
data_y = [0 for x in range(N)]
data_error = [0 for x in range(N)]
names = [0 for x in range(N)]

# Get computational results
for n in range(len(file_paths)):
    with open(file_paths[n]) as input:
        names[n] = next(input)
        names[n] = names[n].replace("#","")
        names[n] = names[n].replace("\n","")
        next(input)
        data = zip(*(line.strip().split() for line in input))
        data_x[n] = np.asfarray(data[0][:])
        data_y[n] = np.asfarray(data[1][:])
        data_error[n] = np.asfarray(data[2][:])*data_y[n]

# Plot
if user_args.m:
  fig = plt.figure(num=1, figsize=(9,7))
else:
  fig = plt.figure(num=1, figsize=(10,7))

# set height ratios for sublots
if user_args.m:
  gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
else:
  gs = gridspec.GridSpec(1, 1)

# the first subplot
ax0 = plt.subplot(gs[0])

x_label = 'Energy (keV)'
plt.xlabel(x_label, size=14)
plt.ylabel('Reflection Coef.', size=14)
plt.title('Electron Albedos for an infinite slab of Al', size=16)
ax=plt.gca()

plt.xlim(9e-2,256.0)
plt.ylim(0.07,0.27)
if user_args.m:
    plt.xlim(9e-2,256.0)
    # plt.ylim(0.1,0.32)
    plt.ylim(0.07,0.25)

markers = ["o","*","v","^","<",">","+","x","1","2","3","4","p","s","h","D","d","H","8","o","*"]
if user_args.e:
  # exp_names = ['assad', 'bienlein','bishop', 'bongeler', 'bronshtein', 'cosslett', 'drescher', 'el_gomati', 'heinrich', 'kanter', 'kulenkampff', 'lockwood', 'neubert', 'palluel', 'philibert', 'reimer', 'shimizu', 'soum', 'trump', 'wittry' ]
  # exp_names = ['assad', 'bienlein','bishop', 'bongeler', 'bronshtein', 'cosslett', 'drescher', 'el_gomati', 'heinrich', 'kanter', 'kulenkampff', 'lockwood', 'neubert', 'reimer', 'shimizu', 'soum', 'trump', 'wittry' ]
  exp_names = ['assad', 'bienlein','bishop', 'bronshtein', 'cosslett', 'drescher', 'el_gomati', 'heinrich', 'kanter', 'kulenkampff', 'lockwood', 'neubert', 'reimer', 'shimizu', 'soum', 'trump', 'wittry' ]

  # Get experimental data
  directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

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
          x[j] = float(data[0][j])
          y[j] = float(data[1][j])

    if i == 1:
      plt.scatter(x, y, label="Experimental", marker=markers[1], s=50, facecolors='none', edgecolors='b' )
    else:
      plt.scatter(x, y, marker=markers[1], s=50, facecolors='none', edgecolors='b' )
ax.set_xscale('log')

if user_args.a:
  # exp_names = ['assad', 'bienlein','bishop', 'bongeler', 'bronshtein', 'cosslett', 'drescher', 'el_gomati', 'heinrich', 'kanter', 'kulenkampff', 'lockwood', 'neubert', 'palluel', 'philibert', 'reimer', 'shimizu', 'soum', 'trump', 'wittry' ]
  exp_names = ['assad', 'bienlein','bishop', 'bongeler', 'bronshtein', 'cosslett', 'drescher', 'el_gomati', 'heinrich', 'kanter', 'kulenkampff', 'lockwood', 'neubert', 'reimer', 'shimizu', 'soum', 'trump', 'wittry' ]
  # exp_names = ['assad', 'bienlein','bishop', 'bronshtein', 'cosslett', 'drescher', 'el_gomati', 'heinrich', 'kanter', 'kulenkampff', 'lockwood', 'neubert', 'shimizu', 'soum', 'trump', 'wittry' ]

  # Get experimental data
  directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

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
          x[j] = float(data[0][j])
          y[j] = float(data[1][j])

    plt.scatter(x, y, label=name, marker=markers[i], s=50, facecolors='none', edgecolors='b' )
ax.set_xscale('log')


# markers = ["v","^","<",">","s","+","x","1","2","3","4","8","p","*","h","H","X","D","d"]
# markerssizes = [7,7,7,7,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7]
# marker_color = ['g', 'r', 'c', 'm', 'y', 'k', 'w', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
markers = ["-s","--v","-.o",":^","--p",":*","-.h","-H","--X",":D","-.d","-8"]
markerssizes = [7,8,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8]
marker_color = ['g', 'r', 'm', 'k', 'y', 'c', 'g', 'r', 'm', 'k', 'y', 'c']

if user_args.m:
    names = ['MCNP6.2','FRENSIE-ACE', 'FRENSIE-ENDL' ]
for n in range(N):
    x = np.asfarray(data_x[n])*1000
    y = np.asfarray(data_y[n])
    yerr = np.asfarray(data_error[n])*y
    plt.errorbar(x, y, yerr=yerr, label=names[n], fmt=markers[n], markersize=markerssizes[n], color=marker_color[n] )

if user_args.a:
  lgd = plt.legend(loc="upper left", bbox_to_anchor=(1,1))
elif user_args.m:
  pylab.legend(loc='best')
  ax0.grid(linestyle=':')
else:
  pylab.legend(loc='best')

ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
# plt.ylim(0.12,0.22)
# leg = plt.legend(loc='best', ncol=2)
# leg.get_frame().set_alpha(0.5)

if user_args.m:

    # The C/E subplot (with shared x-axis)
    ax1 = plt.subplot(gs[1], sharex = ax0)
    plt.xlabel(x_label, size=14)
    plt.ylabel('C/R', size=14)

    # Get mcnp data
    mcnp_x = np.asfarray(data_x[0])*1000
    mcnp_y = np.asfarray(data_y[0])
    mcnp_yerr = np.asfarray(data_error[0])*mcnp_y

    for n in range(1,N):
        x = np.asfarray(data_x[n])*1000
        y = np.asfarray(data_y[n])
        yerr = np.asfarray(data_error[n])*y

        print '\n------ ',names[n], ' ------'
        max_diff = 0.0
        diff_energy = 0.0
        diff_x = []
        diff_y = []
        diff_error = []
        c_index = 0
        r_index = 0
        max_diff = 0
        while True:
          if mcnp_x[r_index] < x[c_index]:
            r_index += 1
          elif mcnp_x[r_index] > x[c_index]:
            c_index += 1
          else:
            # Calculate the relative difference and prop. error
            c_over_r = y[c_index]/mcnp_y[r_index]
            diff = (1.0 - c_over_r)*100.0
            error = np.sqrt( ((1.0/mcnp_y[r_index])**2)*(yerr[c_index])**2 + ((y[c_index]/mcnp_y[r_index]**2)**2)*(mcnp_yerr[r_index])**2 )

            # Add data to list
            diff_x.append( mcnp_x[r_index] )
            diff_y.append( c_over_r )
            diff_error.append( error )

            # Print results
            print mcnp_x[r_index], ": ", diff, "% ", error
            if abs(diff) > abs(max_diff):
              max_diff = diff
              max_energy = mcnp_x[r_index]
              max_error = error

            c_index += 1
            r_index += 1
          if r_index == len(mcnp_x) or c_index == len(x):
            break

        print '--------------------------'
        print 'Max Relative Diff: ',diff_energy, ": ", max_diff, "% ", error,"\n"

        plt.errorbar(diff_x, diff_y, yerr=diff_error, label=names[n], fmt=markers[n], markersize=markerssizes[n], color=marker_color[n] )

        # make x ticks for first suplot invisible
        plt.setp(ax0.get_xticklabels(), visible=False)

        # remove first tick label for the first subplot
        yticks = ax0.yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        ax1.grid(linestyle=':')

        plt.xlim(0.1,300.0)
        # plt.ylim(0.97,1.05)

        # remove vertical gap between subplots
        plt.subplots_adjust(hspace=.0)

ax0.grid(linestyle=':')
plt.show()

output = "albedo_results.pdf"
if user_args.o:
    output = user_args.o

print "Plot outputted to: ",output
if user_args.a:
  fig.savefig(output, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)
else:
  fig.savefig(output, bbox_inches='tight', dpi=300)