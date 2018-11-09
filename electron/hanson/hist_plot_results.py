#! /usr/bin/env python
# Luke Kersting
# This script asks for transmission (Frac/Deg2) data and run names which it then plots.
import numpy as np
from numpy import exp, linspace
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D
import argparse as ap
import inspect, os
from lmfit import Model

def gaussian(x, amp, wid):
  return amp * exp(-x*x /(wid*wid) )

# Set up the argument parser
description = "This script asks for transmission (Frac/Deg2) data and run names which "\
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
M = 15
data_x = [[] for y in range(N)]
data_y = [[] for y in range(N)]
data_error = [[] for y in range(N)]
names = [0 for x in range(N)]
amp = [0 for x in range(N)]
wid = [0 for x in range(N)]
# Get computational results
for n in range(N):

    with open(file_paths[n]) as input:
        names[n] = input.readline()[1:].strip()
        print names[n]
        print input.readline().strip()[1:]
        data = zip(*(line.strip().split('\t') for line in input))
        data_x[n] = np.asfarray(data[0][0:M])
        data_y[n] = np.asfarray(data[1][0:M])
        data_error[n] = np.asfarray(data[2][0:M])*data_y[n]

# Plot
if user_args.m:
  fig = plt.figure(num=1, figsize=(9,7))
  # set height ratios for sublots
  gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
else:
  fig = plt.figure(num=1, figsize=(9,7))
  # set height ratios for sublots
  gs = gridspec.GridSpec(1, 1)

# the first subplot
ax0 = plt.subplot(gs[0])

x_label = 'Angle (Degree)'
plt.xlabel(x_label, size=14)
plt.ylabel('#/Square Degrees', size=14)
plt.title('$\mathrm{15.7\/MeV\/Electron\/Angular\/Distribution\/from\/a\/9.658\/\mu m\/Gold\/Foil}$', size=16)
ax=plt.gca()

plt.xlim(0.0,7.0)
plt.ylim(0.0,0.05)
if user_args.m:
    plt.ylim(0.0,0.05)

plots = []
labels = []
# Experimental values given by Hanson
F0_exp = 0.0441
Theta_exp = 2.58
if user_args.e:
    # Get experimental data
    directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    filename = directory + "/experimental_results.tsv"
    with open(filename) as input:
        data = zip(*(line.strip().split('\t') for line in input))
        exp_x = np.asfarray(data[0][1:14])
        exp_y = np.asfarray(data[1][1:14])

    # Plot the experimental data
    plt1, = plt.plot(exp_x, exp_y, linestyle='None', marker='s', markersize=5 )

    # Get a least squares fit
    gmodel = Model(gaussian)
    result = gmodel.fit(exp_y, x=exp_x, amp=F0_exp, wid=Theta_exp)
    print 'Hanson Fit:\n--------------------------------------------------------'
    F0_fit = result.best_values['amp']
    Theta_fit = result.best_values['wid']
    print 'F(0)         = ',F0_fit,'\tTheta(1/e) = ',Theta_fit
    print 'Rel Diff Fit = ',(F0_exp-F0_fit)/F0_exp,'\tRel Diff   = ',(Theta_exp-Theta_fit)/Theta_exp,'\n'

    x = linspace(0, exp_x[len(exp_x)-1])
    # y = result.eval(x=x)
    # plt2, = plt.plot(x, y, 'b--' )

    # Plot the Gaussian fit given by Hanson
    y = [None] * len(x)
    for i in range(len(y)):
      y[i] = gaussian(x[i], F0_exp, Theta_exp)
    plt3, = plt.plot(x, y, color='b' )

    plots.append( (plt1,plt3) )
    labels.append("Hanson (Exp.)" )
    # plots.append((plt2) )
    # labels.append("Least Square Fit (Hanson)" )


markers = ["--v","-.o",":^","--<","-.>",":+","--x","-.1",":2","--3","-.4",":8","--p","-.P",":*","--h","-.H",":X","--D","-.d"]
linestyle = ["--","-.",":","--","-.",":","--","-.",":","--3","-.",":","--","-.",":","--","-.",":","--","-."]
markerssizes = [6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6]
marker_color = ['g', 'r', 'm', 'k', 'y', 'c', 'g', 'r', 'm', 'k', 'y', 'c']

# linestyles: 'solid', 'dashed', 'dashdotted', 'densely dotted', 'dashdotdotted', 'densely dashed', 'densely dashdotted', 'densely dashdotdotted', 'dotted', 'loosely dashed', 'loosely dashdotted', 'loosely dashdotdotted')

linestyles = [(0, ()), (0, (5, 5)), (0, (3, 5, 1, 5)), (0, (1, 1)), (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (1, 5)), (0, (5, 10)), (0, (3, 10, 1, 10)), (0, (3, 10, 1, 10, 1, 10))]

if user_args.m:
    names = ['MCNP6.2','FRENSIE-ACE', 'FRENSIE-ENDL' ]
# names = ['MCNP6.2','FRENSIE-ACE', 'FRENSIE-ENDL' ]
for n in range(N):
    # Insert first bin lower bounds as an angle of 0
    x = np.insert( data_x[n], 0, 0.0)

    # Plot histogram of results
    m, bins, plt1 = plt.hist(x[:-1], bins=x, weights=data_y[n], histtype='step', label=names[n], color=marker_color[n], linestyle=linestyles[n], linewidth=1.8 )

    # Plot error bars
    mid = 0.5*(bins[1:] + bins[:-1])
    plt2 = plt.errorbar(mid, m, yerr=data_error[n], ecolor=marker_color[n], fmt=None)

    # Get a least squares fit
    gmodel = Model(gaussian)
    result = gmodel.fit(data_y[n], x=mid, amp=F0_exp, wid=Theta_exp)

    print '\n',names[n], ':\n--------------------------------------------------------'
    amp[n] = result.best_values['amp']
    wid[n] = result.best_values['wid']
    print 'F(0)         = ',amp[n],'\tTheta(1/e) = ',wid[n]
    if user_args.e:
        answer = (F0_exp-amp[n])/F0_exp
        print 'Rel Diff Exp = ','%.6e' % answer,'\tRel Diff   = ',(Theta_exp-wid[n])/Theta_exp
        print 'Rel Diff Fit = ',(F0_fit-amp[n])/F0_fit,'\tRel Diff   = ',(Theta_fit-wid[n])/Theta_fit

    handle1 = Line2D([], [], c=marker_color[n], linestyle='--', dashes=linestyles[n][1], linewidth=1.8)
    plots.append( handle1 )
    labels.append(names[n])

plt.legend(plots, labels, loc=1)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
if user_args.m:
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

markers = ["v","o","^","<",">","+","x","1","2","3","4","8","p","P","*","h","H","X","D","d"]
if user_args.m:

    # The C/E subplot (with shared x-axis)
    ax1 = plt.subplot(gs[1], sharex = ax0)
    plt.xlabel(x_label, size=14)
    plt.ylabel('C/R', size=14)

    for n in range(1,N):
        print "\n", names[n]
        # Insert first bin lower bounds as an angle of 0
        x = np.insert( data_x[n][:-2], 0, 0.0)

        # Calculate C/R
        yerr = np.sqrt( ((1.0/data_y[0][:-2])**2)*(data_error[n][:-2])**2 + ((data_y[n][:-2]/data_y[0][:-2]**2)**2)*(data_error[0][:-2])**2 )
        y = data_y[n][:-2]/data_y[0][:-2]

        # Print C/R results
        for i in range(0, len(y)):
          print x[i+1], ": ", (1.0-y[i])*100, u"\u00B1", yerr[i]*100, "%"

        # Plot histogram of results
        m, bins, _ = ax1.hist(x[:-1], bins=x, weights=y, histtype='step', label=names[n], color=marker_color[n], linestyle=linestyles[n], linewidth=1.8 )
        # Plot error bars
        mid = 0.5*(bins[1:] + bins[:-1])
        ax1.errorbar(mid, m, yerr=yerr, ecolor=marker_color[n], fmt=None)

    # make x ticks for first suplot invisible
    plt.setp(ax0.get_xticklabels(), visible=False)

    # remove first tick label for the first subplot
    yticks = ax0.yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    ax0.grid(linestyle=':')
    ax1.grid(linestyle=':')

    plt.xlim(0.0,6.78)
    plt.ylim(0.99,1.01)

    # remove vertical gap between subplots
    plt.subplots_adjust(hspace=.0)

output = "hanson_results.pdf"
if user_args.o:
    output = user_args.o

name = output[:-3] + 'txt'
# Write title to file
out_file = open(name, 'w')
out_file.write( "# " + "Least Squares Fit" +"\n")

# Write data header to file
header = "# Grid Policy\tF(0)\tTheta(1/e)\n"
out_file.write(header)

# Write out fit values
for n in range(N):
    data = names[n] + "\t" + \
            '%.16e' % amp[n] + "\t" + \
            '%.16e' % wid[n] + "\n"
    out_file.write( data )

print "Plot outputted to: ",output
print "Fits outputted to: ",name
fig.savefig(output, bbox_inches='tight', dpi=600)
plt.show()