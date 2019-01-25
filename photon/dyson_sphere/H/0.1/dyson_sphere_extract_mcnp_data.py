#!/usr/bin/python
import sys, os
from optparse import *
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from dyson_sphere_simulation_plot import extractMCNPEstimatorData

if __name__ == "__main__":

    # Parse the command line arguments
    parser = OptionParser()
    parser.add_option("--mcnp_file", type="string", dest="mcnp_file",
                      help="the mcnp file to load")
    parser.add_option("--mcnp_file_start", type="int", dest="mcnp_file_start",
                      help="the mcnp output file start line")
    parser.add_option("--mcnp_file_end", type="int", dest="mcnp_file_end",
                      help="the mcnp output file end line")
    parser.add_option("--output_file", type="string", dest="output_file",
                      help="the output file name")
    options,args = parser.parse_args()

    # Extract the estimator data
    extractMCNPEstimatorData( options.mcnp_file,
                              options.mcnp_file_start,
                              options.mcnp_file_end,
                              options.output_file )
