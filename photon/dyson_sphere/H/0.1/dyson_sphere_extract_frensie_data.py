#!/usr/bin/python
import sys, os
from optparse import *
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from dyson_sphere_simulation_plot import extractFrensieEstimatorData

if __name__ == "__main__":

    # Parse the command line arguments
    parser = OptionParser()
    parser.add_option("--rendezvous_file", type="string", dest="rendezvous_file",
                      help="the rendezvous file to load")
    parser.add_option("--estimator_id", type="int", dest="estimator_id",
                      help="the estimator id to use")
    parser.add_option("--entity_id", type="int", dest="entity_id",
                      help="the entity id to use")
    parser.add_option("--output_file", type="string", dest="output_file",
                      help="the output file name")
    options,args = parser.parse_args()

    # Extract the estimator data
    extractFrensieEstimatorData( options.rendezvous_file,
                                 options.estimator_id,
                                 options.entity_id,
                                 options.output_file )
    
