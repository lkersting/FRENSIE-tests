#!/usr/bin/python
from os import path, environ
from optparse import *
import sys

# Add the parent directory to the path
sys.path.insert(1,path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__))))))
import infinite_medium_simulation as simulation
import PyFrensie.Data as Data
import PyFrensie.Utility as Utility
import PyFrensie.MonteCarlo as MonteCarlo

pyfrensie_path =path.dirname( path.dirname(path.abspath(MonteCarlo.__file__)))
database_path = environ['DATABASE_PATH']

if __name__ == "__main__":

    # Parse the command line options
    parser = OptionParser()
    parser.add_option("--elastic_mode", type="string", dest="mode", default="coupled",
                      help="the elastic distribution mode ( decoupled, coupled, hybrid )")
    parser.add_option("--elastic_method", type="string", dest="method", default="2D",
                      help="the elastic coupled sampling method ( 2D, 1D, modified 2D )")
    parser.add_option("--grid_policy", type="string", dest="grid_policy", default="unit_correlated",
                      help="the bivariate Grid Policy ( unit_correlated, unit vase, correlated )")
    parser.add_option("--threads", type="int", dest="threads", default=4,
                      help="the number of threads to use")
    parser.add_option("--num_particles", type="float", dest="num_particles", default=1e6,
                      help="the number of particles to run")
    parser.add_option("--wall_time", type="float", dest="wall_time", default=1350,
                      help="the wall time for the run")
    parser.add_option("--log_file", type="string", dest="log_file",
                      help="the file that will be used for logging")
    options,args = parser.parse_args()

    # Set the elastic distribution mode ( DECOUPLED, COUPLED, HYBRID )
    mode = simulation.getElasticModeFromString(options.mode)

    # Set the elastic coupled sampling method ( TWO_D, ONE_D, MODIFIED_TWO_D )
    method = simulation.getElasticMethodFromString(options.method)

    # Set the bivariate Grid Policy ( UNIT_BASE_CORRELATED, CORRELATED, UNIT_BASE )
    grid_policy = simulation.getGridPolicyFromString(options.grid_policy)

    # Set the max and cutoff energy
    energy=0.01
    cutoff_energy=1e-4

    # Set the simulation properties
    properties = simulation.setAdjointSimulationProperties( options.num_particles, options.wall_time, mode, method, cutoff_energy, energy )

    # Turn certain reactions off
    properties.setAdjointElasticModeOff()
    properties.setAdjointElectroionizationModeOff()
    properties.setAdjointBremsstrahlungModeOff()
    # properties.setAdjointAtomicExcitationModeOff()

    # Set the simulation name
    sim_name = simulation.setAdjointSimulationName( properties, energy, "Pb", grid_policy )

    # Create the results directory
    simulation.createResultsDirectory(sim_name)

    geometry_path = path.dirname(path.realpath(__file__)) + "/geom.h5m"

    # Set the energy bins
    bins = list(Utility.doubleArrayFromString( "{ 1e-4, 149l, 1e-2}" ))

    # Set the file version
    if grid_policy == MonteCarlo.UNIT_BASE_GRID:
      version = 0
    elif grid_policy == MonteCarlo.UNIT_BASE_CORRELATED_GRID:
      version = 1
    elif grid_policy == MonteCarlo.CORRELATED_GRID:
      version = 2

    # Run the simulation
    simulation.runAdjointUniformEnergyInfiniteMediumSimulation(
          sim_name,
          database_path,
          geometry_path,
          properties,
          cutoff_energy,
          energy,
          bins,
          82000,
          version,
          options.threads,
          options.log_file )