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

# Set database directory path (for Denali)
if socket.gethostname() == "Denali":
  database_path = "/home/software/mcnpdata/database.xml"
else: # Set database directory path (for Cluster)
  database_path = "/home/lkersting/software/mcnp6.2/MCNP_DATA/database.xml"

geometry_path = path.dirname(path.realpath(__file__)) + "/geom.h5m"

##----------------------------------------------------------------------------##
## ----------------------------- RUN SIMULATION ----------------------------- ##
##----------------------------------------------------------------------------##
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
  geometry_type = "DagMC" #(ROOT or DAGMC)

  # Set element zaid and name
  atom=Data.Au_ATOM
  zaid=79000
  element="Au"

  # Set geometry path and type
  model_properties = DagMC.DagMCModelProperties( geometry_path )
  model_properties.useFastIdLookup()

  # Set model
  geom_model = DagMC.DagMCModel( model_properties )

  ##--------------------------------------------------------------------------##
  ## -------------------------- EVENT HANDLER SETUP ------------------------- ##
  ##--------------------------------------------------------------------------##

  # Set event handler
  event_handler = Event.EventHandler( properties )

  ## -------------------- Transmission Current Estimator -------------------- ##

  # Setup a surface current estimator for the transmission current
  estimator_id = 1
  surface_ids = [1]
  transmission_current_estimator = Event.WeightMultipliedSurfaceCurrentEstimator( estimator_id, 1.0, surface_ids )

  # Set the particle type
  transmission_current_estimator.setParticleTypes( [MonteCarlo.ELECTRON] )

  # Set the cosine bins
  cosine_bins_1 = [ -1.000000000000000, 0.000000000000000, 0.939692620785908, 0.965925826289068, 0.984807753012208, 0.990268068741570, 0.994521895368273, 0.995396198367179, 0.996194698091746, 0.996917333733128, 0.997564050259824, 0.998134798421867, 0.998629534754574, 0.999048221581858, 0.999390827019096, 0.999657324975557, 0.999847695156391, 0.999961923064171, 1.000000000000000 ]

  transmission_current_estimator.setCosineDiscretization( cosine_bins_1 )

  # Add the estimator to the event handler
  event_handler.addEstimator( transmission_current_estimator )

  ## --------------------- Reflection Current Estimator --------------------- ##

  # Setup a surface current estimator for the reflection current
  estimator_id = 2
  surface_ids = [2]
  reflection_current_estimator = Event.WeightMultipliedSurfaceCurrentEstimator( estimator_id, 1.0, surface_ids )

  # Set the particle type
  reflection_current_estimator.setParticleTypes( [MonteCarlo.ELECTRON] )

  # Set the cosine bins
  cosine_bins_2 = [ -1.0, -0.999999, 0.0, 1.0 ]
  reflection_current_estimator.setCosineDiscretization( cosine_bins_2 )

  # Add the estimator to the event handler
  event_handler.addEstimator( reflection_current_estimator )

  ## ---------------------- Track Length Flux Estimator --------------------- ##

  # Setup a track length flux estimator
  estimator_id = 3
  cell_ids = [1]
  track_flux_estimator = Event.WeightMultipliedCellTrackLengthFluxEstimator( estimator_id, 1.0, cell_ids, geom_model )

  # Set the particle type
  track_flux_estimator.setParticleTypes( [MonteCarlo.ELECTRON] )

  # Set the energy bins
  energy_bins = numpy.logspace(numpy.log10(1.5e-5), numpy.log10(15.7), num=101) #[ 1.5e-5, 99l, 15.7 ]
  track_flux_estimator.setEnergyDiscretization( energy_bins )

  # Add the estimator to the event handler
  event_handler.addEstimator( track_flux_estimator )

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
  if file_type == Data.ElectroatomicDataProperties.ACE_EPR_FILE:
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
  delta_energy = Distribution.DeltaDistribution( 15.7 )
  energy_dimension_dist = ActiveRegion.IndependentEnergyDimensionDistribution( delta_energy )
  particle_distribution.setDimensionDistribution( energy_dimension_dist )

  # Set the direction dimension distribution
  particle_distribution.setDirection( 0.0, 0.0, 1.0 )

  # Set the spatial dimension distribution
  particle_distribution.setPosition( 0.0, 0.0, -0.1 )

  particle_distribution.constructDimensionDistributionDependencyTree()

  # Set source components
  source_component = [ActiveRegion.StandardElectronSourceComponent( 0, 1.0, geom_model, particle_distribution )]

  # Set source
  source = ActiveRegion.StandardParticleSource( source_component )

  # Set the archive type
  archive_type = "xml"

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
    processCosineBinData( transmission_current_estimator, cosine_bins_1, name, title )

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

  factory = Manager.ParticleSimulationManagerFactory( rendezvous, histories, time, threads )

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

    # Call destructor for manager and factory
    manager = 0
    factory = 0

    print "Processing the results:"
    processData( archive_name )

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

  if not path.exists(directory):
    makedirs(directory)

  return directory

##----------------------------------------------------------------------------##
## -------------------------- setSimulationName ------------------------------##
##----------------------------------------------------------------------------##
# Define a function for naming an electron simulation
def setSimulationName( properties ):

  extension, title = setup.setSimulationNameExtention( properties, file_type )
  name = "hanson" + extension
  output = setup.getResultsDirectory(file_type, interpolation) + "/" + name

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
##------------------------------- processData --------------------------------##
##----------------------------------------------------------------------------##

# This function pulls data from the rendezvous file
def processData( rendezvous_file ):

  Collision.FilledGeometryModel.setDefaultDatabasePath( database_path )

  # Load data from file
  manager = Manager.ParticleSimulationManagerFactory( rendezvous_file ).getManager()
  event_handler = manager.getEventHandler()

  # Get the estimator data
  estimator_1 = event_handler.getEstimator( 1 )
  cosine_bins = estimator_1.getCosineDiscretization()

  # Get the simulation name and title
  properties = manager.getSimulationProperties()

  if "epr14" not in rendezvous_file:
    file_type = Data.ElectroatomicDataProperties.Native_EPR_FILE
  else:
    file_type = Data.ElectroatomicDataProperties.ACE_EPR_FILE

  filename, title = setSimulationName( properties )

  print "Processing the results:"
  processCosineBinData( estimator_1, cosine_bins, filename, title )

  print "Results will be in ", path.dirname(filename)

##----------------------------------------------------------------------------##
##--------------------------- processCosineBinData ---------------------------##
##----------------------------------------------------------------------------##

# This function pulls cosine estimator data outputs it to a separate file.
def processCosineBinData( estimator, cosine_bins, filename, title ):

  ids = list(estimator.getEntityIds() )

  today = datetime.date.today()

  degree = numpy.pi/180.0
  square_degree = degree*degree

  # Read the data file for surface tallies
  name = filename+"_spectrum.txt"
  out_file = open(name, 'w')

  # Get the current and relative error
  processed_data = estimator.getEntityBinProcessedData( ids[0] )
  current = list(processed_data['mean'])
  current_rel_error = list(processed_data['re'])

  # Convert to #/Square Degree
  num_square_degree = [None] * len(current)
  num_square_degree_rel_error = [None] * len(current_rel_error)
  angle_bins = [None] * len(cosine_bins)

  size = len(current)
  for i in range(0, size ):
    k = size - i
    j = k - 1

    # Calculate the angle from the cosine_bins
    angle_bins[i] = numpy.arccos(float(cosine_bins[k]))/degree

    # Calculate the current in 1/square degrees
    cosine_diff = float(cosine_bins[k]) - float(cosine_bins[j])
    sterradians = 2.0*numpy.pi*cosine_diff
    num_per_ster = float(current[j])/sterradians
    num_square_degree[i] = num_per_ster*square_degree

    # Set the relative error
    num_square_degree_rel_error[i] = float(current_rel_error[j])

  # Set the last angle bin boundary
  angle_bins[size] = numpy.arccos(float(cosine_bins[0]))/degree

  # Write title to file
  out_file.write( "# " + title +"\n")
  # Write data header to file
  header = "# Degrees\tTransmission (Frac/Deg2)\tError\t"+str(today)+"\n"
  out_file.write(header)

  # Write data to file
  for i in range(0, size):
      output = '%.4e' % angle_bins[i] + "\t" + \
              '%.16e' % num_square_degree[i] + "\t" + \
              '%.16e' % num_square_degree_rel_error[i] + "\n"
      out_file.write( output )

  # Write the last angle bin boundary
  output = '%.4e' % angle_bins[size] + "\n"
  out_file.write( output )
  out_file.close()