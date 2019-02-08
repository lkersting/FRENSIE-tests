#! /usr/bin/env python
from os import path, makedirs
import sys
import numpy as np
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

##---------------------------------------------------------------------------##
## Set up and run the forward simulation
def runForwardAlbedoSimulation( sim_name,
                                db_path,
                                geom_name,
                                properties,
                                source_energy,
                                zaid,
                                file_type,
                                version,
                                threads,
                                log_file = None ):


    ## Initialize the MPI session
    session = MPI.GlobalMPISession( len(sys.argv), sys.argv )

    # Suppress logging on all procs except for the master (proc=0)
    Utility.removeAllLogs()
    session.initializeLogs( 0, True )

    if not log_file is None:
        session.initializeLogs( log_file, 0, True )

    if session.rank() == 0:
      print "The PyFrensie path is set to: ", pyfrensie_path

    simulation_properties = properties

  ##--------------------------------------------------------------------------##
  ## ---------------------------- MATERIALS SETUP --------------------------- ##
  ##--------------------------------------------------------------------------##

    # Set element name
    element_name="albedo"

    ## Set up the materials
    database = Data.ScatteringCenterPropertiesDatabase( db_path )

    # Extract the properties for the zaid from the database
    atom_properties = database.getAtomProperties( Data.ZAID(zaid) )

    # Set the definition for the zaid for this simulation
    scattering_center_definitions = Collision.ScatteringCenterDefinitionDatabase()
    atom_definition = scattering_center_definitions.createDefinition( element_name, Data.ZAID(zaid) )

    atom_definition.setElectroatomicDataProperties(
      atom_properties.getSharedElectroatomicDataProperties( file_type, version ) )

    # Set the definition for material 1
    material_definitions = Collision.MaterialDefinitionDatabase()
    material_definitions.addDefinition( element_name, 1, [element_name], [1.0] )

  ##--------------------------------------------------------------------------##
  ## ---------------------------- GEOMETRY SETUP ---------------------------- ##
  ##--------------------------------------------------------------------------##

    ## Set up the geometry
    model_properties = DagMC.DagMCModelProperties( geom_name )
    model_properties.useFastIdLookup()

    # Load the model
    model = DagMC.DagMCModel( model_properties )

    # Fill the model with the defined material
    filled_model = Collision.FilledGeometryModel( db_path, scattering_center_definitions, material_definitions, simulation_properties, model, True )

  ##--------------------------------------------------------------------------##
  ## ----------------------------- SOURCE SETUP ----------------------------- ##
  ##--------------------------------------------------------------------------##

    ## Set up the source
    particle_distribution = ActiveRegion.StandardParticleDistribution( "mono-energetic beam dist" )

    particle_distribution.setEnergy( source_energy )
    particle_distribution.setPosition( 0.0, 0.0, -0.1 )
    particle_distribution.setDirection( 0.0, 0.0, 1.0 )
    particle_distribution.constructDimensionDistributionDependencyTree()

    # The generic distribution will be used to generate electrons
    electron_distribution = [ActiveRegion.StandardElectronSourceComponent( 0, 1.0, model, particle_distribution )]

    # Assign the electron source component to the source
    source = ActiveRegion.StandardParticleSource( electron_distribution )

  ##--------------------------------------------------------------------------##
  ## -------------------------- EVENT HANDLER SETUP ------------------------- ##
  ##--------------------------------------------------------------------------##

    ## Set up the event handler
    event_handler = Event.EventHandler( model, simulation_properties )

  ##----------------------- Surface Current Estimators -----------------------##

    current_estimator = event_handler.getEstimator( 1 )

    # Set the cosine bins
    cosine_bins = [ -1.0, -0.99, 0.0, 1.0 ]
    current_estimator.setCosineDiscretization( cosine_bins )


  ##--------------------------------------------------------------------------##
  ## ----------------------- SIMULATION MANAGER SETUP ----------------------- ##
  ##--------------------------------------------------------------------------##

    # Set the archive type
    archive_type = "xml"

    ## Set up the simulation manager
    factory = Manager.ParticleSimulationManagerFactory( filled_model,
                                                        source,
                                                        event_handler,
                                                        simulation_properties,
                                                        sim_name,
                                                        archive_type,
                                                        threads )

    # Create the simulation manager
    manager = factory.getManager()

    # Allow logging on all procs
    session.restoreOutputStreams()

    ## Run the simulation
    if session.size() == 1:
        manager.runInterruptibleSimulation()
    else:
        manager.runSimulation()

    if session.rank() == 0:

      # Get the plot title and filename
      title = setup.getSimulationPlotTitle( sim_name )

      print "Processing the results:"
      processCosineBinData( current_estimator, source_energy, sim_name, title )

      print "Results will be in ", path.dirname(path.abspath(sim_name))

##---------------------------------------------------------------------------##
## Set up and run the adjoint simulation
def runAdjointAlbedoSimulation( sim_name,
                                db_path,
                                geom_name,
                                properties,
                                cutoff_energy,
                                max_energy,
                                zaid,
                                version,
                                threads,
                                log_file = None ):


    ## Initialize the MPI session
    session = MPI.GlobalMPISession( len(sys.argv), sys.argv )

    # Suppress logging on all procs except for the master (proc=0)
    Utility.removeAllLogs()
    session.initializeLogs( 0, True )

    if not log_file is None:
        session.initializeLogs( log_file, 0, True )

    if session.rank() == 0:
      print "The PyFrensie path is set to: ", pyfrensie_path

    simulation_properties = properties

  ##--------------------------------------------------------------------------##
  ## ---------------------------- MATERIALS SETUP --------------------------- ##
  ##--------------------------------------------------------------------------##

    # Set element name
    element_name="albedo"

    ## Set up the materials
    database = Data.ScatteringCenterPropertiesDatabase( db_path )

    # Extract the properties for the zaid from the database
    atom_properties = database.getAtomProperties( Data.ZAID(zaid) )

    # Set the definition for the zaid for this simulation
    scattering_center_definitions = Collision.ScatteringCenterDefinitionDatabase()
    atom_definition = scattering_center_definitions.createDefinition( element_name, Data.ZAID(zaid) )

    file_type = Data.AdjointElectroatomicDataProperties.Native_EPR_FILE

    atom_definition.setAdjointElectroatomicDataProperties(
      atom_properties.getSharedAdjointElectroatomicDataProperties( file_type, version ) )

    # Set the definition for material 1
    material_definitions = Collision.MaterialDefinitionDatabase()
    material_definitions.addDefinition( element_name, 1, [element_name], [1.0] )

  ##--------------------------------------------------------------------------##
  ## ---------------------------- GEOMETRY SETUP ---------------------------- ##
  ##--------------------------------------------------------------------------##

    ## Set up the geometry
    model_properties = DagMC.DagMCModelProperties( geom_name )
    model_properties.useFastIdLookup()

    # Load the model
    model = DagMC.DagMCModel( model_properties )

    # Fill the model with the defined material
    filled_model = Collision.FilledGeometryModel( db_path, scattering_center_definitions, material_definitions, simulation_properties, model, True )

  ##--------------------------------------------------------------------------##
  ## ----------------------------- SOURCE SETUP ----------------------------- ##
  ##--------------------------------------------------------------------------##

    ## Set up the source
    particle_distribution = ActiveRegion.StandardParticleDistribution( "adjoint uniform energy dist" )

    # Uniform distribution from cutoff energy to max problem energy
    uniform_energy = Distribution.UniformDistribution( cutoff_energy, max_energy )
    energy_dimension_dist = ActiveRegion.IndependentEnergyDimensionDistribution( uniform_energy )
    particle_distribution.setDimensionDistribution( energy_dimension_dist )

    # Slightly to the left (negative z-direction) of the semi-infinite slab
    particle_distribution.setPosition( 0.0, 0.0, -0.1 )

    # Uniform distribution for all angles in the positive z direction
    uniform_positive_mu = Distribution.UniformDistribution( 0.0, 1.0, 1.0 )
    mu_dimension_dist = ActiveRegion.IndependentSecondaryDirectionalDimensionDistribution( uniform_positive_mu )
    particle_distribution.setDimensionDistribution( mu_dimension_dist )

    particle_distribution.constructDimensionDistributionDependencyTree()

    # The generic distribution will be used to generate electrons
    electron_distribution = [ActiveRegion.StandardAdjointElectronSourceComponent( 0, 1.0, model, particle_distribution )]

    # Assign the electron source component to the source
    source = ActiveRegion.StandardParticleSource( electron_distribution )

  ##--------------------------------------------------------------------------##
  ## -------------------------- EVENT HANDLER SETUP ------------------------- ##
  ##--------------------------------------------------------------------------##

    ## Set up the event handler
    event_handler = Event.EventHandler( model, simulation_properties )

  ##----------------------- Surface Current Estimators -----------------------##

    current_estimator = event_handler.getEstimator( 2 )

    # Set the energy bins (for each cosine bin)
    bin_string = "{ " + str(cutoff_energy) + ", 149l, " + str(max_energy) +" }"
    energy_bins = list(Utility.doubleArrayFromString( bin_string ))
    current_estimator.setEnergyDiscretization( energy_bins )

    # Set the cosine bins
    cosine_bins = [ 1.0, np.cos(np.deg2rad(10)), np.cos(np.deg2rad(50)), np.cos(np.deg2rad(70)), np.cos(np.deg2rad(110)), np.cos(np.deg2rad(130)), np.cos(np.deg2rad(170)), -1.0 ]
    current_estimator.setCosineDiscretization( cosine_bins )

  ##--------------------------------------------------------------------------##
  ## ----------------------- SIMULATION MANAGER SETUP ----------------------- ##
  ##--------------------------------------------------------------------------##

    # Set the archive type
    archive_type = "xml"

    ## Set up the simulation manager
    factory = Manager.ParticleSimulationManagerFactory( filled_model,
                                                        source,
                                                        event_handler,
                                                        simulation_properties,
                                                        sim_name,
                                                        archive_type,
                                                        threads )

    # Create the simulation manager
    manager = factory.getManager()

    # Allow logging on all procs
    session.restoreOutputStreams()

    ## Run the simulation
    if session.size() == 1:
        manager.runInterruptibleSimulation()
    else:
        manager.runSimulation()

    if session.rank() == 0:

      # Get the plot title and filename
      title = "FRENSIE - Adjoint"

      print "Processing the results:"
      processCosineEnergyBinData( current_estimator, sim_name, title )

      print "Results will be in ", path.dirname(path.abspath(sim_name))

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

  # Set database directory path (for Denali)
  if socket.gethostname() == "Denali":
    database_path = "/home/software/mcnpdata/database.xml"
  # Set database directory path (for Elbrus)
  elif socket.gethostname() == "Elbrus":
    database_path = "/home/software/mcnpdata/database.xml"
  else: # Set database directory path (for Cluster)
    database_path = "/home/lkersting/software/mcnp6.2/MCNP_DATA/database.xml"

  # Set the data path
  Collision.FilledGeometryModel.setDefaultDatabasePath( database_path )

  print rendezvous

  time_sec = time*60
  factory = Manager.ParticleSimulationManagerFactory( rendezvous, histories, time_sec, threads )

  manager = factory.getManager()
  manager.setSimulationName( rendezvous )

  Utility.removeAllLogs()
  session.initializeLogs( 0, False )

  manager.runSimulation()

  if session.rank() == 0:

    # Get the event handler
    event_handler = manager.getEventHandler()

    # Get the estimator data
    estimator_1 = event_handler.getEstimator( 1 )

    # Get the simulation name and title
    properties = manager.getSimulationProperties()

    # Get the plot title and filename
    title = setup.getSimulationPlotTitle( rendezvous )
    filename = rendezvous.split("_rendezvous_")[0]
    energy = float(rendezvous.split("_")[2])
    print str(energy)

    print "Processing the results:"
    processCosineBinData( estimator_1, energy, filename, title )

    print "Results will be in ", path.dirname(filename)

##----------------------------------------------------------------------------##
## ------------------------ Create Results Directory ------------------------ ##
##----------------------------------------------------------------------------##
def createResultsDirectory(file_type, interpolation, element):

  directory = setup.getResultsDirectory(file_type, interpolation)

  directory = element + "/" + directory

  if not path.exists(directory):
    makedirs(directory)

##---------------------------------------------------------------------------##
## -------------------------- setSimulationName -----------------------------##
##---------------------------------------------------------------------------##
# Define a function for naming an electron simulation
def setSimulationName( properties, file_type, element, energy, refined ):
  extension = setup.setSimulationNameExtention( properties, file_type )
  name = "albedo_" + element + "_" + str(energy)
  if refined:
    name += "_refined"
  name += extension
  interpolation = properties.getElectronTwoDInterpPolicy()
  output = setup.getResultsDirectory(file_type, interpolation) + "/" + name

  return output

##---------------------------------------------------------------------------##
## ---------------------- setAdjointSimulationName --------------------------##
##---------------------------------------------------------------------------##
# Define a function for naming an electron simulation
def setAdjointSimulationName( properties, element, grid_policy, ionization_sampling, nudge_past_max_energy ):
  extension = setup.setAdjointSimulationNameExtention( properties )
  name = "adjoint_albedo_" + element + "_"

  # Add the grid policy to the name
  if grid_policy == MonteCarlo.UNIT_BASE_CORRELATED_GRID:
      name += "unit_correlated"
  elif grid_policy == MonteCarlo.CORRELATED_GRID:
      name += "correlated"
  else:
      name += "unit_base"

  # Add the ionization sampling to the name
  if ionization == MonteCarlo.OUTGOING_ENERGY_SAMPLING:
    name += '_outgoing_energy'

  # Add the nudge past max energy mode to the name
  if not nudge_past_max_energy:
    name += '_no_nudge'

  name += extension

  date = str(datetime.datetime.today()).split()[0]
  directory = "results/adjoint/" + date + "/"

  output = directory + "/" + name

  return output

##----------------------------------------------------------------------------##
##------------------------------- processData --------------------------------##
##----------------------------------------------------------------------------##

# This function pulls data from the rendezvous file
def processData( rendezvous_file ):

  # Set database directory path (for Denali)
  if socket.gethostname() == "Denali":
    database_path = "/home/software/mcnpdata/database.xml"
  # Set database directory path (for Elbrus)
  elif socket.gethostname() == "Elbrus":
    database_path = "/home/software/mcnpdata/database.xml"
  else: # Set database directory path (for Cluster)
    database_path = "/home/lkersting/software/mcnp6.2/MCNP_DATA/database.xml"

  Collision.FilledGeometryModel.setDefaultDatabasePath( database_path )

  # Load data from file
  manager = Manager.ParticleSimulationManagerFactory( rendezvous_file ).getManager()
  event_handler = manager.getEventHandler()

  # Get the estimator data
  estimator_1 = event_handler.getEstimator( 1 )

  # Get the plot title and filename
  title = setup.getSimulationPlotTitle( rendezvous_file )
  filename = rendezvous_file.split("_rendezvous_")[0]
  energy = float(rendezvous_file.split("_")[2])

  print "Processing the results:\n"
  processCosineBinData( estimator_1, energy, filename, title )

  print "\nResults will be in ", path.dirname(filename)

##----------------------------------------------------------------------------##
##--------------------------- processCosineBinData ---------------------------##
##----------------------------------------------------------------------------##

# This function pulls cosine estimator data outputs it to a separate file.
def processCosineBinData( estimator, energy, filename, title ):

  ids = list(estimator.getEntityIds() )
  if not 2 in ids:
    print "ERROR: estimator does not contain entity 2!"
    raise ValueError(message)

  today = datetime.date.today()

  # Read the data file for surface tallies
  name = filename+"_albedo.txt"
  out_file = open(name, 'w')

  # Get the current and relative error
  processed_data = estimator.getEntityBinProcessedData( 2 )
  current = processed_data['mean']
  current_rel_error = processed_data['re']
  cosine_bins = estimator.getCosineDiscretization()

  print current
  print cosine_bins
  # Write title to file
  out_file.write( "# " + title +"\n")
  # Write data header to file
  header = "# Energy (MeV)\tAlbedo\tError\t"+str(today)+"\n"
  out_file.write(header)

  # Write data to file
  output = '%.6e' % energy + "\t" + \
           '%.16e' % current[-1] + "\t" + \
           '%.16e' % current_rel_error[-1] + "\n"
  out_file.write( output )
  out_file.close()

##----------------------------------------------------------------------------##
##------------------------ processCosineEnergyBinData ------------------------##
##----------------------------------------------------------------------------##

# This function pulls cosine energy bin estimator data outputs it to a separate file.
def processCosineEnergyBinData( estimator, filename, title ):

  ids = list(estimator.getEntityIds() )
  if not 2 in ids:
    print "ERROR: estimator does not contain entity 2!"
    raise ValueError(message)

  today = datetime.date.today()

  # Read the data file for surface tallies
  name = filename+"_albedo.txt"
  out_file = open(name, 'w')

  # Get the current and relative error
  processed_data = estimator.getEntityBinProcessedData( 2 )
  current = processed_data['mean']
  current_rel_error = processed_data['re']
  cosine_bins = estimator.getCosineDiscretization()
  energy_bins = estimator.getEnergyDiscretization()

  print current
  print cosine_bins
  print energy_bins
  # Write title to file
  out_file.write( "# " + title +"\n")
  # Write data header to file
  header = "# Angle\tEnergy (MeV)\tCurrent\tError\t"+str(today)+"\n"
  out_file.write(header)

  # Write data to file
  for i in range(1, len(cosine_bins) ):
    angle = '%.6e' % cosine_bins[i-1] + " - ", '%.6e' % cosine_bins[i]
    out_file.write( angle )
    for j in range(0, len(energy_bins) ):
      output = '%.6e' % energy_bins[j] + "\t" + \
              '%.16e' % current[i-1+j] + "\t" + \
              '%.16e' % current_rel_error[i-1+j] + "\n"
      out_file.write( output )
  out_file.close()
