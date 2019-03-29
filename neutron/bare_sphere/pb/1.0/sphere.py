#!/usr/bin/python
import sys, os
from optparse import *
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from sphere_simulation import sphereSimulation
import PyFrensie.Utility as Utility


if __name__ == "__main__":

        # Parse the command line options
    parser = OptionParser()
    parser.add_option("--threads", type="int", dest="threads", default=1,
                      help="the number of threads to use")
    parser.add_option("--db_path", type="string", dest="db_path",
                      help="the database name (with extension)")
    parser.add_option("--sim_name", type="string", dest="sim_name", default="sphere",
                      help="the simulation name")
    parser.add_option("--log_file", type="string", dest="log_file",
                      help="the file that will be used for logging")
    parser.add_option("--num_particles", type="float", dest="num_particles", default=1e3,
                      help="the number of particles to run")
    parser.add_option("--source_energy", type="float", dest="source_energy",default=1.0,
                      help="the source energy in MeV")
    options,args = parser.parse_args()

    if options.db_path is None:
        print "The database path must be specified!"
        sys.exit(1)

    # Run the simulation
    sphereSimulation( options.sim_name,
                      options.db_path,
                      options.num_particles,
                      options.source_energy,
                      Utility.doubleArrayFromString("{1e-9,100l,1.0}"),
                      options.threads,
                      options.log_file )
