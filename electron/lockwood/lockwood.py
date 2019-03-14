#! /usr/bin/env python
from os import path, makedirs, environ
import sys
import numpy
import errno
import datetime

# Add the parent directory to the path
sys.path.insert(1,path.dirname(path.dirname(path.abspath(__file__))))
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
database_path = environ['DATABASE_PATH']

##----------------------------------------------------------------------------##
## ---------------------- GLOBAL SIMULATION VARIABLES ----------------------- ##
##----------------------------------------------------------------------------##

# Set the element
atom=Data.Al_ATOM; element="Al"; zaid=13000
# Set the source energy
energy=0.314
# Set the test number
test_number=0

# Set the bivariate interpolation (LOGLOGLOG, LINLINLIN, LINLINLOG)
interpolation=MonteCarlo.LOGLOGLOG_INTERPOLATION

# Set the bivariate Grid Policy (UNIT_BASE_CORRELATED, CORRELATED, UNIT_BASE)
grid_policy=MonteCarlo.UNIT_BASE_CORRELATED_GRID

# Set the elastic distribution mode ( DECOUPLED, COUPLED, HYBRID )
mode=MonteCarlo.COUPLED_DISTRIBUTION

# Set the elastic coupled sampling method
# ( TWO_D_UNION, ONE_D_UNION, MODIFIED_TWO_D_UNION )
method=MonteCarlo.MODIFIED_TWO_D_UNION

# Set the data file type (ACE_EPR_FILE, Native_EPR_FILE)
file_type=Data.ElectroatomicDataProperties.Native_EPR_FILE

# Set if a refined grid should be used ( True, False )
use_refined_grid=False

# Set the calorimeter thickness (g/cm2)
calorimeter_thickness=5.050E-03

# Set the range (g/cm2)
test_range=0.0025

geometry_path = path.dirname(path.realpath(__file__)) + "/"
geometry_path += element + "/" + element + "_" + str(energy) + "/dagmc/geom_" + str(test_number) + ".h5m"

# Run the simulation
def runSimulation( threads, histories, time ):

  ##--------------------------------------------------------------------------##
  ## ------------------------------ MPI Session ----------------------------- ##
  ##--------------------------------------------------------------------------##
  session = MPI.GlobalMPISession( len(sys.argv), sys.argv )
  Utility.removeAllLogs()
  session.initializeLogs( 0, True )

  if session.rank() == 0:
    print "The PyFrensie path is set to: ", pyfrensie_path

  properties = setSimulationProperties( histories, time )

  ##--------------------------------------------------------------------------##
  ## ---------------------------- GEOMETRY SETUP ---------------------------- ##
  ##--------------------------------------------------------------------------##


  # Set geometry path and type
  model_properties = DagMC.DagMCModelProperties( geometry_path )
  model_properties.useFastIdLookup()

  # Construct model
  geom_model = DagMC.DagMCModel( model_properties )

  ##--------------------------------------------------------------------------##
  ## -------------------------- EVENT HANDLER SETUP ------------------------- ##
  ##--------------------------------------------------------------------------##

  # Set event handler
  event_handler = Event.EventHandler( properties )

  ## -------------------- Energy Deposition Calorimeter --------------------- ##

  # Setup a cell pulse height estimator
  estimator_id = 1
  cell_ids = [2]
  energy_deposition_estimator = Event.WeightAndEnergyMultipliedCellPulseHeightEstimator( estimator_id, 1.0, cell_ids )

  # Add the estimator to the event handler
  event_handler.addEstimator( energy_deposition_estimator )

  ##--------------------------------------------------------------------------##
  ## ----------------------- SIMULATION MANAGER SETUP ----------------------- ##
  ##--------------------------------------------------------------------------##

  # Initialized database
  database = Data.ScatteringCenterPropertiesDatabase(database_path)
  scattering_center_definition_database = Collision.ScatteringCenterDefinitionDatabase()

  # Set element properties
  element_properties = database.getAtomProperties( atom )

  element_definition = scattering_center_definition_database.createDefinition( element, Data.ZAID(zaid) )

  version = 0
  if use_refined_grid:
    if grid_policy == MonteCarlo.UNIT_BASE_GRID:
      version = 1
    elif grid_policy == MonteCarlo.UNIT_BASE_CORRELATED_GRID:
      version = 2
    elif grid_policy == MonteCarlo.CORRELATED_GRID:
      version = 3
  elif file_type == Data.ElectroatomicDataProperties.ACE_EPR_FILE:
    version = 14

  element_definition.setElectroatomicDataProperties(
            element_properties.getSharedElectroatomicDataProperties( file_type, version ) )

  material_definition_database = Collision.MaterialDefinitionDatabase()
  material_definition_database.addDefinition( element, 1, (element,), (1.0,) )

  # Fill model
  model = Collision.FilledGeometryModel( database_path, scattering_center_definition_database, material_definition_database, properties, geom_model, True )

  # Set particle distribution
  particle_distribution = ActiveRegion.StandardParticleDistribution( "source distribution" )

  # Set the energy dimension distribution
  delta_energy = Distribution.DeltaDistribution( energy )
  energy_dimension_dist = ActiveRegion.IndependentEnergyDimensionDistribution( delta_energy )
  particle_distribution.setDimensionDistribution( energy_dimension_dist )

  # Set the direction dimension distribution
  particle_distribution.setDirection( 0.0, 0.0, 1.0 )

  # Set the spatial dimension distribution
  particle_distribution.setPosition( 0.0, 0.0, 0.0 )

  particle_distribution.constructDimensionDistributionDependencyTree()

  # Set source components
  source_component = [ActiveRegion.StandardElectronSourceComponent( 0, 1.0, geom_model, particle_distribution )]

  # Set source
  source = ActiveRegion.StandardParticleSource( source_component )

  # Set the archive type
  archive_type = "xml"

  # Set the simulation name and title
  name = setSimulationName( properties, use_refined_grid )
  title = setup.getSimulationPlotTitle( name )

  factory = Manager.ParticleSimulationManagerFactory( model,
                                                      source,
                                                      event_handler,
                                                      properties,
                                                      name,
                                                      archive_type,
                                                      threads )

  manager = factory.getManager()

  Utility.removeAllLogs()
  session.initializeLogs( 0, False )

  manager.runSimulation()

  if session.rank() == 0:

    print "Processing the results:"
    processData( energy_deposition_estimator, name, title, test_range, calorimeter_thickness )

    print "Results will be in ", path.dirname(name)

##----------------------------------------------------------------------------##
## --------------------- Run Simulation From Rendezvous --------------------- ##
##----------------------------------------------------------------------------##
def runSimulationFromRendezvous( threads,
                                 histories,
                                 time,
                                 rendezvous,
                                 range_for_test,
                                 cal_thickness,
                                 log_file = None,
                                 num_rendezvous = None ):

  ##--------------------------------------------------------------------------##
  ## ------------------------------ MPI Session ----------------------------- ##
  ##--------------------------------------------------------------------------##
  session = MPI.GlobalMPISession( len(sys.argv), sys.argv )

  # Suppress logging on all procs except for the master (proc=0)
  Utility.removeAllLogs()
  session.initializeLogs( 0, True )

  if session.rank() == 0:
      print "The PyFrensie path is set to: ", pyfrensie_path

  if not log_file is None:
      session.initializeLogs( log_file, 0, True )

  # Set the data path
  Collision.FilledGeometryModel.setDefaultDatabasePath( database_path )

  time_sec = time*60
  if not num_rendezvous is None:
      new_simulation_properties = MonteCarlo.SimulationGeneralProperties()
      new_simulation_properties.setNumberOfHistories( int(histories) )
      new_simulation_properties.setMinNumberOfRendezvous( int(num_rendezvous) )
      new_simulation_properties.setSimulationWallTime( float(time_sec) )

      factory = Manager.ParticleSimulationManagerFactory( rendezvous,
                                                          new_simulation_properties,
                                                          threads )
  else:
      factory = Manager.ParticleSimulationManagerFactory( rendezvous,
                                                          int(histories),
                                                          float(time_sec),
                                                          threads )

  manager = factory.getManager()

  manager.initialize()

  # Allow logging on all procs
  session.restoreOutputStreams()

  ## Run the simulation
  if session.size() == 1:
      manager.runInterruptibleSimulation()
  else:
      manager.runSimulation()

  if session.rank() == 0:

      # Get the event handler
      event_handler = manager.getEventHandler()
      sim_name = rendezvous.split("_rendezvous_")[0]

      # Get the estimator data
      estimator_1 = event_handler.getEstimator( 1 )


      title = setup.getSimulationPlotTitle( sim_name )

      print "Processing the results:"
      processData( estimator_1, sim_name, title, range_for_test, cal_thickness )

      print "Results will be in ", path.dirname(sim_name)

##----------------------------------------------------------------------------##
## ------------------------- SIMULATION PROPERTIES -------------------------- ##
##----------------------------------------------------------------------------##
def setSimulationProperties( histories, time ):

  properties = setup.setSimulationProperties( histories, time, interpolation, grid_policy, mode, method )


  ## -------------------------- ELECTRON PROPERTIES ------------------------- ##

  # Turn certain reactions off
  # properties.setElasticModeOff()
  # properties.setElectroionizationModeOff()
  # properties.setBremsstrahlungModeOff()
  # properties.setAtomicExcitationModeOff()

  return properties

##----------------------------------------------------------------------------##
## ------------------------ Create Results Directory ------------------------ ##
##----------------------------------------------------------------------------##
def createResultsDirectory():

  directory = setup.getResultsDirectory(file_type, interpolation)

  directory = element + "/" + directory + "/"

  if grid_policy == MonteCarlo.UNIT_BASE_CORRELATED_GRID:
    policy="unit_correlated"
  elif grid_policy == MonteCarlo.UNIT_BASE_GRID:
    policy="unit_base"
  elif grid_policy == MonteCarlo.CORRELATED_GRID:
    policy="correlated"

  directory += policy

  if not path.exists(directory):
    try:
        makedirs(directory)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass

  print directory
  return directory

##----------------------------------------------------------------------------##
## -------------------------- setSimulationName ------------------------------##
##----------------------------------------------------------------------------##
# Define a function for naming an electron simulation
def setSimulationName( properties, refined ):
  extension = setup.setSimulationNameExtention( properties, file_type )
  name = "lockwood_" + element + "_" + str(energy) + "_" + str(test_number)
  if refined:
    name += "_refined"
  name += extension

  if grid_policy == MonteCarlo.UNIT_BASE_CORRELATED_GRID:
    policy="unit_correlated"
  elif grid_policy == MonteCarlo.UNIT_BASE_GRID:
    policy="unit_base"
  elif grid_policy == MonteCarlo.CORRELATED_GRID:
    policy="correlated"

  output = element + "/" + setup.getResultsDirectory(file_type, interpolation) + "/" + policy + "/" + name

  return output

##----------------------------------------------------------------------------##
## -------------------------- getSimulationName ------------------------------##
##----------------------------------------------------------------------------##
# Define a function for naming an electron simulation
def getSimulationName():

  properties = setSimulationProperties( 1, 1.0 )

  name = setSimulationName( properties, use_refined_grid )
  title = setup.getSimulationPlotTitle( name )

  return name

##----------------------------------------------------------------------------##
##----------------------- processDataFromRendezvous --------------------------##
##----------------------------------------------------------------------------##

# This function pulls pulse height estimator data outputs it to a separate file.
def processDataFromRendezvous( rendezvous_file, range, calorimeter_thickness ):

  # Activate just-in-time initialization to prevent automatic loading of the
  # geometry and data tables
  Utility.activateJustInTimeInitialization()

  # Set the database path
  Collision.FilledGeometryModel.setDefaultDatabasePath( database_path )

  # Load data from file
  manager = Manager.ParticleSimulationManagerFactory( rendezvous_file ).getManager()
  event_handler = manager.getEventHandler()

  # Get the estimator data
  estimator_1 = event_handler.getEstimator( 1 )

  # Get the simulation name and title
  properties = manager.getSimulationProperties()

  filename = rendezvous_file.split("_rendezvous_")[0]
  title = setup.getSimulationPlotTitle( filename )

  print "Processing the results:"
  processData( estimator_1, filename, title, range, calorimeter_thickness )

  print "Results will be in ", path.dirname(filename)

##----------------------------------------------------------------------------##
##------------------------------- processData --------------------------------##
##----------------------------------------------------------------------------##

# This function pulls pulse height estimator data outputs it to a separate file.
def processData( estimator, filename, title, range, calorimeter_thickness ):

  ids = list(estimator.getEntityIds())

  processed_data = estimator.getEntityBinProcessedData( ids[0] )
  energy_dep_mev = processed_data['mean'][0]
  rel_error = processed_data['re'][0]

  today = datetime.date.today()

  # Read the data file for surface tallies
  name = filename+"_energy_dep.txt"
  out_file = open(name, 'w')

  # Write the header to the file
  header = "# Range (g/cm2)\tEnergy Deposition (MeV cm2/g)\tError\t"+str(today)+"\n"
  out_file.write(header)

  # Write the energy deposition to the file
  data = str(range) + '\t' + str(energy_dep_mev/calorimeter_thickness) + '\t' + str(rel_error/calorimeter_thickness)
  out_file.write(data)
  out_file.close()
