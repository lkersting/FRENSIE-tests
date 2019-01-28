#! /usr/bin/env python
from os import path, makedirs
import sys
import numpy
import datetime
import socket

# Add the parent directory to the path
sys.path.insert(1,path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))
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

##----------------------------------------------------------------------------##
## ---------------------- GLOBAL SIMULATION VARIABLES ----------------------- ##
##----------------------------------------------------------------------------##

# Set the element
atom=Data.H_ATOM; element="H"; zaid=1000
# Set the forward source max energy
max_energy=0.01

# Set the min energy (default is 100 eV )
energy_cutoff=1e-4

# Set the elastic distribution mode ( DECOUPLED, COUPLED, HYBRID )
mode=MonteCarlo.DECOUPLED_DISTRIBUTION

# Set the elastic coupled sampling method
# ( TWO_D_UNION, ONE_D_UNION, MODIFIED_TWO_D_UNION )
method=MonteCarlo.MODIFIED_TWO_D_UNION

# Set the ionization sampling method
# ( KNOCK_ON_SAMPLING, OUTGOING_ENERGY_SAMPLING )
ionization=MonteCarlo.KNOCK_ON_SAMPLING

# Set the bivariate Grid Policy ( 'UNIT_BASE_CORRELATED', 'UNIT_BASE' )
grid_policy='UNIT_BASE_CORRELATED'

# Set the nudge past max energy mode on/off (true/false)
nudge_past_max_energy = True

# Set database directory path (for Denali)
if socket.gethostname() == "Denali":
  database_path = "/home/software/mcnpdata/database.xml"
# Set database directory path (for Elbrus)
elif socket.gethostname() == "Elbrus":
  database_path = "/home/software/mcnpdata/database.xml"
# Set database directory path (for Cluster)
else:
  database_path = "/home/lkersting/software/mcnp6.2/MCNP_DATA/database.xml"

geometry_path = path.dirname(path.realpath(__file__)) + "/geom.h5m"

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

  # Set the energy bins
  bins = list(Utility.doubleArrayFromString( "{ 1e-4, 399l, 1e-2}" ))

  ## ------------------------ Surface Flux Estimator ------------------------ ##

  # Setup an adjoint surface flux estimator
  estimator_id = 2
  surface_ids = [1, 16, 18]
  surface_flux_estimator = Event.WeightMultipliedSurfaceFluxEstimator( estimator_id, 1.0, surface_ids, geom_model )

  # Set the particle type
  surface_flux_estimator.setParticleTypes( [MonteCarlo.ADJOINT_ELECTRON] )

  # Set the energy bins
  surface_flux_estimator.setSourceEnergyDiscretization( bins )

  # Create response function
  uniform_energy = Distribution.UniformDistribution( energy_cutoff, max_energy, 1.0 )
  particle_response_function = ActiveRegion.EnergyParticleResponseFunction( uniform_energy )
  response_function = ActiveRegion.StandardParticleResponse( particle_response_function )

  # Set the response function
  surface_flux_estimator.setResponseFunctions( [response_function] )

  # Add the estimator to the event handler
  event_handler.addEstimator( surface_flux_estimator )

  # ## -------------------------- Particle Tracker ---------------------------- ##

  # particle_tracker = Event.ParticleTracker( 0, 5 )

  # # Add the particle tracker to the event handler
  # event_handler.addParticleTracker( particle_tracker )

  ##--------------------------------------------------------------------------##
  ## ----------------------- SIMULATION MANAGER SETUP ----------------------- ##
  ##--------------------------------------------------------------------------##

  # Initialized database
  database = Data.ScatteringCenterPropertiesDatabase(database_path)
  scattering_center_definition_database = Collision.ScatteringCenterDefinitionDatabase()

  # Set element properties
  element_properties = database.getAtomProperties( atom )

  element_definition = scattering_center_definition_database.createDefinition( element, Data.ZAID(zaid) )

  if grid_policy == 'UNIT_BASE_CORRELATED':
    version = 0
  elif grid_policy == 'UNIT_BASE':
    version = 2

  if not nudge_past_max_energy:
    version += 1

  if ionization == MonteCarlo.OUTGOING_ENERGY_SAMPLING:
    version += 4

  file_type = Data.AdjointElectroatomicDataProperties.Native_EPR_FILE

  element_definition.setAdjointElectroatomicDataProperties(
            element_properties.getSharedAdjointElectroatomicDataProperties( file_type, version ) )

  material_definition_database = Collision.MaterialDefinitionDatabase()
  material_definition_database.addDefinition( element, 1, (element,), (1.0,) )

  # Fill model
  model = Collision.FilledGeometryModel( database_path, scattering_center_definition_database, material_definition_database, properties, geom_model, True )

  # Set particle distribution
  particle_distribution = ActiveRegion.StandardParticleDistribution( "source distribution" )

  # Set the energy dimension distribution
  uniform_energy = Distribution.UniformDistribution( energy_cutoff, max_energy )
  energy_dimension_dist = ActiveRegion.IndependentEnergyDimensionDistribution( uniform_energy )
  particle_distribution.setDimensionDistribution( energy_dimension_dist )

  # Set the spatial dimension distribution
  particle_distribution.setPosition( 0.0, 0.0, 0.0 )

  particle_distribution.constructDimensionDistributionDependencyTree()

  # Set source components
  source_component = [ActiveRegion.StandardAdjointElectronSourceComponent( 0, 1.0, model, particle_distribution )]

  # Set source
  source = ActiveRegion.StandardParticleSource( source_component )

  # Set the archive type
  archive_type = "xml"

  # Set the simulation name and title
  name, title = setSimulationName( properties )

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
    processData( event_handler, name, title )
    print "Results will be in ", path.dirname(name)

##----------------------------------------------------------------------------##
## --------------------- Run Simulation From Rendezvous --------------------- ##
##----------------------------------------------------------------------------##
def runSimulationFromRendezvous( threads, histories, time, rendezvous ):

  ##--------------------------------------------------------------------------##
  ## ------------------------------ MPI Session ----------------------------- ##
  ##--------------------------------------------------------------------------##
  session = MPI.GlobalMPISession( len(sys.argv), sys.argv )
  Utility.removeAllLogs()
  session.initializeLogs( 0, True )

  if session.rank() == 0:
    print "The PyFrensie path is set to: ", pyfrensie_path

  # Set the data path
  Collision.FilledGeometryModel.setDefaultDatabasePath( database_path )

  time_sec = time*60
  factory = Manager.ParticleSimulationManagerFactory( rendezvous, histories, time_sec, threads )

  manager = factory.getManager()

  Utility.removeAllLogs()
  session.initializeLogs( 0, False )

  manager.runSimulation()

  if session.rank() == 0:

    rendezvous_number = manager.getNumberOfRendezvous()

    components = rendezvous.split("rendezvous_")
    archive_name = components[0] + "rendezvous_"
    archive_name += str( rendezvous_number - 1 )
    archive_name += "."
    archive_name += components[1].split(".")[1]

    # Get the event handler
    event_handler = manager.getEventHandler()

    # Get the simulation name and title
    properties = manager.getSimulationProperties()

    filename, title = setSimulationName( properties )

    print "Processing the results:"
    processData( event_handler, filename, title )

    print "Results will be in ", path.dirname(filename)

##----------------------------------------------------------------------------##
## ------------------------- SIMULATION PROPERTIES -------------------------- ##
##----------------------------------------------------------------------------##
def setSimulationProperties( histories, time ):

  properties = setup.setAdjointSimulationProperties( histories, time, mode, method )

  ## -------------------------- ELECTRON PROPERTIES ------------------------- ##

  # Set the min electron energy in MeV (Default is 100 eV)
  properties.setMinAdjointElectronEnergy( energy_cutoff )

  # Set the max electron energy in MeV (Default is 20 MeV)
  properties.setMaxAdjointElectronEnergy( max_energy )

  # Set the cutoff weight properties for rouletting
  properties.setAdjointElectronRouletteThresholdWeight( 1e-8 )
  properties.setAdjointElectronRouletteSurvivalWeight( 1e-6 )

  # Turn certain reactions off
  # properties.setAdjointElasticModeOff()
  # properties.setAdjointElectroionizationModeOff()
  # properties.setAdjointBremsstrahlungModeOff()
  # properties.setAdjointAtomicExcitationModeOff()

  return properties

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
  name = "adjoint_" + str(max_energy) + "_" + grid_policy

  if ionization == MonteCarlo.OUTGOING_ENERGY_SAMPLING:
    name += '_outgoing_energy'
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

##----------------------------------------------------------------------------##
##------------------------ processDataFromRendezvous -------------------------##
##----------------------------------------------------------------------------##

# This function pulls data from the rendezvous file
def processDataFromRendezvous( rendezvous_file ):

  Collision.FilledGeometryModel.setDefaultDatabasePath( database_path )

  # Load data from file
  manager = Manager.ParticleSimulationManagerFactory( rendezvous_file ).getManager()
  event_handler = manager.getEventHandler()

  # Get the simulation name and title
  properties = manager.getSimulationProperties()

  filename, title = setSimulationName( properties )

  print "Processing the results:"
  processData( event_handler, filename, title )

  print "Results will be in ", path.dirname(filename)

##----------------------------------------------------------------------------##
##------------------------------- processData --------------------------------##
##----------------------------------------------------------------------------##
def processData( event_handler, filename, title ):

  # Process surface flux data
  surface_flux = event_handler.getEstimator( 2 )

  file1 = filename + "_1"
  setup.processSurfaceFluxSourceEnergyBinData( surface_flux, 1, file1, title )

  file2 = filename + "_2"
  setup.processSurfaceFluxSourceEnergyBinData( surface_flux, 18, file2, title )

  file3 = filename + "_5"
  setup.processSurfaceFluxSourceEnergyBinData( surface_flux, 16, file3, title )

##----------------------------------------------------------------------------##
##------------------------ printParticleTrackInfo -------------------------##
##----------------------------------------------------------------------------##

# This function pulls data from the rendezvous file
def printParticleTrackInfo( rendezvous_file ):

  Collision.FilledGeometryModel.setDefaultDatabasePath( database_path )

  # Load data from file
  manager = Manager.ParticleSimulationManagerFactory( rendezvous_file ).getManager()
  event_handler = manager.getEventHandler()

  # Get the simulation name and title
  properties = manager.getSimulationProperties()

  if "epr14" not in rendezvous_file:
    file_type = Data.ElectroatomicDataProperties.Native_EPR_FILE
  else:
    file_type = Data.ElectroatomicDataProperties.ACE_EPR_FILE

  filename, title = setSimulationName( properties )

  # Process surface flux data
  particle_tracker = event_handler.getParticleTracker( 0 )

  history_map = particle_tracker.getHistoryData()
  print "len(history_map) = ", len(history_map)
  print "len(history_map[0]) = ", len(list(history_map[0]))

  print particle_tracker.getTrackedHistories()
  print list(history_map)

  cached_particle_state = None
  for i in history_map:
    print "\nHistory number:", i
    if MonteCarlo.ADJOINT_ELECTRON in history_map[i]:
      map_i = history_map[i][MonteCarlo.ADJOINT_ELECTRON]
      for j in range(len(map_i)):
        print "  j:",j
        for k in range(len(map_i[j])):
          print "    k:",k
          print "state:\tenergy\t\t\t\tweight\t\tlocation\t\t\t\tdirection"
          cached_particle_state = map_i[j][k]

          for l in range(len(cached_particle_state)):
            if l < 5:
              location = list(cached_particle_state[l][0])
              direction = list(cached_particle_state[l][1])
              energy = cached_particle_state[l][2]
              time = cached_particle_state[l][3]
              weight = cached_particle_state[l][4]
              collision = cached_particle_state[l][5]
              print l,":\t",'%.20e' % energy,"\t",'%.6e' % weight,"\t",location,"\t",direction,"\t",collision
