#!/usr/bin/python
import sys, os
from optparse import *
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from dyson_sphere_simulation_plot import plotDysonSphereSimulationSpectrum

if __name__ == "__main__":

    # Parse the command line arguments
    parser = OptionParser()
    parser.add_option("--rendezvous_file", type="string", dest="rendezvous_file",
                      help="the rendezvous file to load")
    parser.add_option("--estimator_id", type="int", dest="estimator_id",
                      help="the estimator id to use")
    parser.add_option("--entity_id", type="int", dest="entity_id",
                      help="the entity id to use")
    parser.add_option("--mcnp_file", type="string", dest="mcnp_file",
                      help="the mcnp output file to load")
    parser.add_option("--mcnp_file_start", type="int", dest="mcnp_file_start",
                      help="the mcnp output file start line")
    parser.add_option("--mcnp_file_end", type="int", dest="mcnp_file_end",
                      help="the mcnp output file end line")
    parser.add_option("--emin", type="float", dest="emin",
                      help="the min energy to use")
    parser.add_option("--emax", type="float", dest="emax",
                      help="the max energy to use")
    parser.add_option("--top_ymax", type="float", dest="top_ymax",
                      help="the top plot max y value")
    parser.add_option("--bottom_ymin", type="float", dest="bottom_ymin",
                      help="the bottom plot min y value")
    parser.add_option("--bottom_ymax", type="float", dest="bottom_ymax",
                      help="the bottom plot max y value")
    parser.add_option("--legend_xpos", type="float", dest="legend_xpos",
                      help="the legend x position")
    parser.add_option("--legend_ypos", type="float", dest="legend_ypos",
                      help="the legend y position")
    options,args = parser.parse_args()

    # Plot the spectrum
    plotDysonSphereSimulationSpectrum( options.rendezvous_file,
                                       options.estimator_id,
                                       options.entity_id,
                                       options.mcnp_file,
                                       options.mcnp_file_start,
                                       options.mcnp_file_end,
                                       False,
                                       top_ylims = [0, options.top_ymax],
                                       bottom_ylims = [options.bottom_ymin, options.bottom_ymax],
                                       xlims = [options.emin, options.emax],
                                       legend_pos = [options.legend_xpos,options.legend_ypos] )
