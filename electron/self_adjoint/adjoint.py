#! /usr/bin/env python
from os import path, makedirs
import sys
import numpy
import datetime
import socket

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

##----------------------------------------------------------------------------##
## ---------------------- GLOBAL SIMULATION VARIABLES ----------------------- ##
##----------------------------------------------------------------------------##

# Set the element
atom=Data.H_ATOM; element="H"; zaid=1000
# Set the forward source energy ( 0.001, 0.01, 0.1 )
energy=0.01

# Set the min energy (default is 100 eV )
energy_cutoff=1e-4

# Set the elastic distribution mode ( DECOUPLED, COUPLED, HYBRID )
mode=MonteCarlo.DECOUPLED_DISTRIBUTION

# Set the elastic coupled sampling method
# ( TWO_D_UNION, ONE_D_UNION, MODIFIED_TWO_D_UNION )
method=MonteCarlo.MODIFIED_TWO_D_UNION

# Set the bivariate Grid Policy ( 'UNIT_BASE_CORRELATED', 'CORRELATED', 'UNIT_BASE' )
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

geometry_path = path.dirname(path.realpath(__file__)) + "/geom_"

# Set the energy bins
if energy == 0.1:
  geometry_path += "100keV.h5m"
elif energy == 0.01:
  geometry_path += "10keV.h5m"
elif energy == 0.001:
  geometry_path += "1keV.h5m"
else:
  print "ERROR: energy ", energy, " not supported!"

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
  if energy == 0.1:
    bins = list(Utility.doubleArrayFromString( "{ 1e-4, 5e-4, 198i, 1e-1}" ))
  elif energy == 0.01:
    bins = list(Utility.doubleArrayFromString( "{ 1e-4, 58i, 6e-3, 99i, 1e-2}" ))
  elif energy == 0.001:
    bins = list(Utility.doubleArrayFromString( "{ 1e-4, 197i, 1e-3}" ))
  else:
    print "ERROR: energy ", energy, " not supported!"

  ## ------------------------ Surface Flux Estimator ------------------------ ##

  # Setup an adjoint surface flux estimator
  estimator_id = 2
  surface_ids = [1]
  surface_flux_estimator = Event.WeightMultipliedSurfaceFluxEstimator( estimator_id, 1.0, surface_ids, geom_model )

  # Set the particle type
  surface_flux_estimator.setParticleTypes( [MonteCarlo.ADJOINT_ELECTRON] )

  # Set the energy bins
  surface_flux_estimator.setSourceEnergyDiscretization( bins )

  # Create response function
  delta_energy = Distribution.DeltaDistribution( energy )
  particle_response_function = ActiveRegion.EnergyParticleResponseFunction( delta_energy )
  response_function = ActiveRegion.StandardParticleResponse( particle_response_function )

  # Set the response function
  surface_flux_estimator.setResponseFunctions( [response_function] )

  # Add the estimator to the event handler
  event_handler.addEstimator( surface_flux_estimator )

  ## -------------------------- Particle Tracker ---------------------------- ##

  # particle_tracker = Event.ParticleTracker( 0, 1000 )

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
  elif grid_policy == 'CORRELATED':
    version = 4

  if not nudge_past_max_energy:
    version += 1

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
  uniform_energy = Distribution.UniformDistribution( energy_cutoff, energy )
  energy_dimension_dist = ActiveRegion.IndependentEnergyDimensionDistribution( uniform_energy )
  particle_distribution.setDimensionDistribution( energy_dimension_dist )

  # Set the spatial dimension distribution
  particle_distribution.setPosition( 0.0, 0.0, 0.0 )

  particle_distribution.constructDimensionDistributionDependencyTree()


  source_critical_line = [ 1.0e-2, 9.98014149e-03, 9.96028344e-03, 9.94042584e-03, 9.92056868e-03, 9.90071198e-03, 9.88085572e-03, 9.86099992e-03, 9.84114457e-03, 9.82128967e-03, 9.80143523e-03, 9.78158124e-03, 9.76172770e-03, 9.74187464e-03, 9.72202203e-03, 9.70216987e-03, 9.68231810e-03, 9.66246694e-03, 9.64261610e-03, 9.62276580e-03, 9.60291601e-03, 9.58306660e-03, 9.56321765e-03, 9.54336920e-03, 9.52352121e-03, 9.50367370e-03 ]

  # Set source components
  source_component = [ActiveRegion.StandardAdjointElectronSourceComponent( 0, 1.0, geom_model, particle_distribution, source_critical_line )]

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
  ids = list( surface_flux.getEntityIds() )
  setup.processSurfaceFluxSourceEnergyBinData( surface_flux, ids[0], filename, title )
