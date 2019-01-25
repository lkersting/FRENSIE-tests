#!/usr/bin/python
from os import path, makedirs
import sys
import numpy
import datetime
import socket

# Add the parent directory to the path
sys.path.insert(1,path.dirname(path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))))
import simulation_setup as setup
import PyFrensie.Data as Data
import PyFrensie.Data.Native as Native
import PyFrensie.Geometry.DagMC as DagMC
import PyFrensie.Geometry as Geometry
import PyFrensie.Utility as Utility
import PyFrensie.Utility.MPI as MPI
import PyFrensie.Utility.Prng as Prng
import PyFrensie.Utility.Coordinate as Coordinate
import PyFrensie.Utility.Distribution as Distribution
import PyFrensie.MonteCarlo as MonteCarlo
import PyFrensie.MonteCarlo.Collision as Collision
import PyFrensie.MonteCarlo.ActiveRegion as ActiveRegion
import PyFrensie.MonteCarlo.Event as Event
import PyFrensie.MonteCarlo.Manager as Manager

pyfrensie_path =path.dirname( path.dirname(path.abspath(MonteCarlo.__file__)))
# from os import path, makedirs
# import sys
# import numpy
# import datetime
from optparse import *
# import socket

# # Add the parent directory to the path
# sys.path.insert(1,path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))
# sys.path.insert(1,path.dirname(path.dirname(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))))
# import simulation_setup as setup
# from infinite_medium_simulation import runAdjointInfiniteMediumSimulation, processAdjointDataFromRendezvous
# import PyFrensie.Data as Data
# import PyFrensie.Data.Native as Native
# import PyFrensie.Geometry.DagMC as DagMC
# import PyFrensie.Geometry as Geometry
# import PyFrensie.Utility as Utility
# import PyFrensie.Utility.MPI as MPI
# import PyFrensie.Utility.Prng as Prng
# import PyFrensie.Utility.Coordinate as Coordinate
# import PyFrensie.Utility.Distribution as Distribution
# import PyFrensie.MonteCarlo as MonteCarlo
# import PyFrensie.MonteCarlo.Collision as Collision
# import PyFrensie.MonteCarlo.ActiveRegion as ActiveRegion
# import PyFrensie.MonteCarlo.Event as Event
# import PyFrensie.MonteCarlo.Manager as Manager

##----------------------------------------------------------------------------##
## ---------------------- GLOBAL SIMULATION VARIABLES ----------------------- ##
##----------------------------------------------------------------------------##

# Set the source and cutoff energy
energy=0.01
cutoff_energy = 1e-4

# Set the elastic distribution mode ( DECOUPLED, COUPLED, HYBRID )
mode=MonteCarlo.DECOUPLED_DISTRIBUTION

# Set the elastic coupled sampling method
# ( TWO_D_UNION, ONE_D_UNION, MODIFIED_TWO_D_UNION )
method=MonteCarlo.MODIFIED_TWO_D_UNION

# Set the bivariate Grid Policy ( UNIT_BASE_CORRELATED, CORRELATED, UNIT_BASE )
grid_policy=MonteCarlo.UNIT_BASE_CORRELATED_GRID

##----------------------------------------------------------------------------##
## ------------------------- SIMULATION PROPERTIES -------------------------- ##
##----------------------------------------------------------------------------##
def setSimulationProperties( histories, time ):

  properties = setup.setAdjointSimulationProperties( histories, time, mode, method )

  ## -------------------------- ELECTRON PROPERTIES ------------------------- ##

  # Set the min electron energy in MeV (Default is 100 eV)
  properties.setMinAdjointElectronEnergy( cutoff_energy )

  # Set the max electron energy in MeV (Default is 20 MeV)
  properties.setMaxAdjointElectronEnergy( energy )

  # Set the critical line energies
  properties.setCriticalAdjointElectronLineEnergies( [energy] )

  # Set the cutoff weight properties for rouletting
  properties.setAdjointElectronRouletteThresholdWeight( 1e-8 )
  properties.setAdjointElectronRouletteSurvivalWeight( 1e-6 )

  # Turn certain reactions off
  # properties.setAdjointElasticModeOff()
  # properties.setAdjointElectroionizationModeOff()
  # properties.setAdjointBremsstrahlungModeOff()
  # properties.setAdjointAtomicExcitationModeOff()

##----------------------------------------------------------------------------##
## ------------------------ Create Results Directory ------------------------ ##
##----------------------------------------------------------------------------##
def createResultsDirectory():

  date = str(datetime.datetime.today()).split()[0]
  directory = "results/adjoint/" + date

  if not path.exists(directory):
    makedirs(directory)

  print directory
  return directory

##----------------------------------------------------------------------------##
## -------------------------- setSimulationName ------------------------------##
##----------------------------------------------------------------------------##
# Define a function for naming an electron simulation
def setSimulationName( properties ):
  extension, title = setup.setAdjointSimulationNameExtention( properties )
  name = "adjoint_" + str(energy) + "_" + grid_policy
  if nudge_past_max_energy:
    name += '_nudged_past_max'
  name += extension
  date = str(datetime.datetime.today()).split()[0]

  output = "results/adjoint/" + date + "/" + name

  return (output, title)

##----------------------------------------------------------------------------##
## -------------------------- getSimulationName ------------------------------##
##----------------------------------------------------------------------------##
# Define a function for naming an electron simulation
def getSimulationName():

  properties = setSimulationProperties( 1, 1.0 )

  name, title = setSimulationName( properties )

  return name


if __name__ == "__main__":

    # Parse the command line options
    parser = OptionParser()
    parser.add_option("--threads", type="int", dest="threads", default=4,
                      help="the number of threads to use")
    parser.add_option("--num_particles", type="float", dest="num_particles", default=1e6,
                      help="the number of particles to run")
    parser.add_option("--wall_time", type="float", dest="wall_time", default=1350,
                      help="the wall time for the run")
    parser.add_option("--log_file", type="string", dest="log_file",
                      help="the file that will be used for logging")
    options,args = parser.parse_args()


    # Set database directory path (for Denali)
    if socket.gethostname() == "Denali":
      database_path = "/home/software/mcnpdata/database.xml"
    else: # Set database directory path (for Cluster)
      database_path = "/home/lkersting/software/mcnp6.2/MCNP_DATA/database.xml"

    geometry_path = path.dirname(path.realpath(__file__)) + "/geom.h5m"

    # Set the bivariate interpolation (LOGLOGLOG, LINLINLIN)
    interpolation=MonteCarlo.LOGLOGLOG_INTERPOLATION

    # Set the data file type (ACE_EPR_FILE, Native_EPR_FILE)
    file_type=Data.ElectroatomicDataProperties.Native_EPR_FILE

    # Set the file version
    version = 1

    properties = setSimulationProperties( options.num_particles, options.wall_time )
    print properties.isAdjointElasticModeOn()
    # Set the simulation name and title
    name, title = setSimulationName( properties )

    geometry_path = path.dirname(path.realpath(__file__)) + "/geom.h5m"

    surface_ids = [1, 16, 18]

    # Set the energy bins
    bins = list(Utility.doubleArrayFromString( "{ 1e-4, 88i, 9e-3, 199i, 1e-2}" ))

    # Set the source critical lines
    source_critical_line = [ 1.0e-2, 9.98014149e-03, 9.96028344e-03, 9.94042584e-03, 9.92056868e-03, 9.90071198e-03, 9.88085572e-03, 9.86099992e-03, 9.84114457e-03, 9.82128967e-03, 9.80143523e-03, 9.78158124e-03, 9.76172770e-03, 9.74187464e-03, 9.72202203e-03, 9.70216987e-03, 9.68231810e-03, 9.66246694e-03, 9.64261610e-03, 9.62276580e-03, 9.60291601e-03, 9.58306660e-03, 9.56321765e-03, 9.54336920e-03, 9.52352121e-03, 9.50367370e-03 ]


    # Run the simulation
    runAdjointInfiniteMediumSimulation( name,
                                        database_path,
                                        geometry_path,
                                        options.num_particles,
                                        properties,
                                        cutoff_energy,
                                        energy,
                                        source_critical_line,
                                        surface_ids,
                                        bins,
                                        1000,
                                        version,
                                        options.threads,
                                        options.log_file )

    print "Processing the results:"
    N = properties.getMinNumberOfRendezvous()
    rendezvous = name + "_rendezvous_" + str(N) + ".xml"
    processAdjointDataFromRendezvous( rendezvous )
    print "Results will be in ", path.dirname(name)