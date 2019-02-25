#! /usr/bin/env python
from os import path, makedirs, environ
import sys
import numpy
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

##---------------------------------------------------------------------------##
## Set up and run the forward simulation
def runForwardDeltaEnergyInfiniteMediumSimulation( sim_name,
                                                   db_path,
                                                   geom_name,
                                                   properties,
                                                   source_energy,
                                                   energy_bins,
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
    element_name="Infinite"

    ## Set up the materials
    database = Data.ScatteringCenterPropertiesDatabase( db_path )

    # Extract the properties for the atom from the database
    atom_properties = database.getAtomProperties( Data.ZAID(zaid) )

    # Set the definition for the atom for this simulation
    scattering_center_definitions = Collision.ScatteringCenterDefinitionDatabase()
    atom_definition = scattering_center_definitions.createDefinition( element_name, Data.ZAID(zaid) )

    file_type = Data.ElectroatomicDataProperties.Native_EPR_FILE

    atom_definition.setElectroatomicDataProperties( atom_properties.getSharedElectroatomicDataProperties( file_type, version ) )

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
    particle_distribution = ActiveRegion.StandardParticleDistribution( "isotropic mono-energetic dist" )

    particle_distribution.setEnergy( source_energy )
    particle_distribution.setPosition( 0.0, 0.0, 0.0 )
    particle_distribution.constructDimensionDistributionDependencyTree()

    # The generic distribution will be used to generate electrons
    electron_distribution = ActiveRegion.StandardElectronSourceComponent( 0, 1.0, model, particle_distribution )

    # Assign the electron source component to the source
    source = ActiveRegion.StandardParticleSource( [electron_distribution] )

  ##--------------------------------------------------------------------------##
  ## -------------------------- EVENT HANDLER SETUP ------------------------- ##
  ##--------------------------------------------------------------------------##

    ## Set up the event handler
    event_handler = Event.EventHandler( model, simulation_properties )

  ## ------------------------ Surface Flux Estimator ------------------------ ##

    surface_flux_estimator = event_handler.getEstimator( 1 )

    # Set the energy bins
    surface_flux_estimator.setEnergyDiscretization( energy_bins )


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

##---------------------------------------------------------------------------##
## Set up and run the adjoint simulation
def runAdjointDeltaEnergyInfiniteMediumSimulation( sim_name,
                                                   db_path,
                                                   geom_name,
                                                   properties,
                                                   energy_cutoff,
                                                   source_energy,
                                                   source_critical_line,
                                                   energy_bins,
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
    element_name="Infinite"

    ## Set up the materials
    database = Data.ScatteringCenterPropertiesDatabase( db_path )

    # Extract the properties for H from the database
    atom_properties = database.getAtomProperties( Data.ZAID(zaid) )

    # Set the definition for H for this simulation
    scattering_center_definitions = Collision.ScatteringCenterDefinitionDatabase()
    atom_definition = scattering_center_definitions.createDefinition( element_name, Data.ZAID(zaid) )

    file_type = Data.AdjointElectroatomicDataProperties.Native_EPR_FILE

    atom_definition.setAdjointElectroatomicDataProperties( atom_properties.getSharedAdjointElectroatomicDataProperties( file_type, version ) )

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
    particle_distribution = ActiveRegion.StandardParticleDistribution( "adjoint isotropic uniform energy dist" )

    uniform_energy = Distribution.UniformDistribution( energy_cutoff, source_energy )
    energy_dimension_dist = ActiveRegion.IndependentEnergyDimensionDistribution( uniform_energy )
    particle_distribution.setDimensionDistribution( energy_dimension_dist )
    particle_distribution.setPosition( 0.0, 0.0, 0.0 )
    particle_distribution.constructDimensionDistributionDependencyTree()

    # The generic distribution will be used to generate electron
    if source_critical_line == None:
      adjoint_electron_distribution = ActiveRegion.StandardAdjointElectronSourceComponent( 0, 1.0, filled_model, particle_distribution )
    else:
      adjoint_electron_distribution = ActiveRegion.StandardAdjointElectronSourceComponent( 0, 1.0, model, particle_distribution, source_critical_line )

    # Assign the electron source component to the source
    source = ActiveRegion.StandardParticleSource( [adjoint_electron_distribution] )

  ##--------------------------------------------------------------------------##
  ## -------------------------- EVENT HANDLER SETUP ------------------------- ##
  ##--------------------------------------------------------------------------##

    ## Set up the event handler
    event_handler = Event.EventHandler( model, simulation_properties )

  ## ------------------------ Surface Flux Estimator ------------------------ ##

    # Setup an adjoint surface flux estimator
    surface_flux_estimator = event_handler.getEstimator( 2 )

    # Create and set the estimator response function
    response_function = ActiveRegion.EnergyParticleResponseFunction( Distribution.DeltaDistribution( source_energy, source_energy - energy_cutoff ) )
    response = ActiveRegion.StandardParticleResponse( response_function )
    surface_flux_estimator.setResponseFunctions( [response] )

    # Set the energy bin discretization
    surface_flux_estimator.setSourceEnergyDiscretization( energy_bins )

  ## -------------------------- Particle Tracker ---------------------------- ##

    # particle_tracker = Event.ParticleTracker( 0, 20 )

    # # Add the particle tracker to the event handler
    # event_handler.addParticleTracker( particle_tracker )

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

##---------------------------------------------------------------------------##
## Set up and run the forward simulation
def runForwardUniformEnergyInfiniteMediumSimulation( sim_name,
                                                     db_path,
                                                     geom_name,
                                                     properties,
                                                     min_energy,
                                                     max_energy,
                                                     energy_bins,
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
    element_name="Infinite"

    ## Set up the materials
    database = Data.ScatteringCenterPropertiesDatabase( db_path )

    # Extract the properties for the atom from the database
    atom_properties = database.getAtomProperties( Data.ZAID(zaid) )

    # Set the definition for the atom for this simulation
    scattering_center_definitions = Collision.ScatteringCenterDefinitionDatabase()
    atom_definition = scattering_center_definitions.createDefinition( element_name, Data.ZAID(zaid) )

    file_type = Data.ElectroatomicDataProperties.Native_EPR_FILE

    atom_definition.setElectroatomicDataProperties( atom_properties.getSharedElectroatomicDataProperties( file_type, version ) )

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

    # Set up the source
    particle_distribution = ActiveRegion.StandardParticleDistribution( "isotropic uniform source dist" )

    # Set the energy dimension distribution
    uniform_energy = Distribution.UniformDistribution( min_energy, max_energy )
    energy_dimension_dist = ActiveRegion.IndependentEnergyDimensionDistribution( uniform_energy )
    particle_distribution.setDimensionDistribution( energy_dimension_dist )
    particle_distribution.setPosition( 0.0, 0.0, 0.0 )
    particle_distribution.constructDimensionDistributionDependencyTree()

    # Set source components
    electron_distribution = ActiveRegion.StandardElectronSourceComponent( 0, 1.0, model, particle_distribution )

    # Assign the electron source component to the source
    source = ActiveRegion.StandardParticleSource( [electron_distribution] )

  ##--------------------------------------------------------------------------##
  ## -------------------------- EVENT HANDLER SETUP ------------------------- ##
  ##--------------------------------------------------------------------------##

    ## Set up the event handler
    event_handler = Event.EventHandler( model, simulation_properties )

  ## ------------------------ Surface Flux Estimator ------------------------ ##

    surface_flux_estimator = event_handler.getEstimator( 1 )

    # Set the energy bins
    surface_flux_estimator.setEnergyDiscretization( energy_bins )


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

##---------------------------------------------------------------------------##
## Set up and run the adjoint simulation
def runAdjointUniformEnergyInfiniteMediumSimulation( sim_name,
                                                     db_path,
                                                     geom_name,
                                                     properties,
                                                     energy_cutoff,
                                                     max_energy,
                                                     energy_bins,
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
    element_name="Infinite"

    ## Set up the materials
    database = Data.ScatteringCenterPropertiesDatabase( db_path )

    # Extract the properties for H from the database
    atom_properties = database.getAtomProperties( Data.ZAID(zaid) )

    # Set the definition for H for this simulation
    scattering_center_definitions = Collision.ScatteringCenterDefinitionDatabase()
    atom_definition = scattering_center_definitions.createDefinition( element_name, Data.ZAID(zaid) )

    file_type = Data.AdjointElectroatomicDataProperties.Native_EPR_FILE

    atom_definition.setAdjointElectroatomicDataProperties( atom_properties.getSharedAdjointElectroatomicDataProperties( file_type, version ) )

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
    particle_distribution = ActiveRegion.StandardParticleDistribution( "adjoint isotropic uniform energy dist" )

    uniform_energy = Distribution.UniformDistribution( energy_cutoff, max_energy )
    energy_dimension_dist = ActiveRegion.IndependentEnergyDimensionDistribution( uniform_energy )
    particle_distribution.setDimensionDistribution( energy_dimension_dist )
    particle_distribution.setPosition( 0.0, 0.0, 0.0 )
    particle_distribution.constructDimensionDistributionDependencyTree()

    # The generic distribution will be used to generate electron
    adjoint_electron_distribution = ActiveRegion.StandardAdjointElectronSourceComponent( 0, 1.0, filled_model, particle_distribution )

    # Assign the electron source component to the source
    source = ActiveRegion.StandardParticleSource( [adjoint_electron_distribution] )

  ##--------------------------------------------------------------------------##
  ## -------------------------- EVENT HANDLER SETUP ------------------------- ##
  ##--------------------------------------------------------------------------##

    ## Set up the event handler
    event_handler = Event.EventHandler( model, simulation_properties )

  ## ------------------------ Surface Flux Estimator ------------------------ ##

    # Setup an adjoint surface flux estimator
    surface_flux_estimator = event_handler.getEstimator( 2 )

    # Create and set the estimator response function
    response_function = ActiveRegion.EnergyParticleResponseFunction( Distribution.UniformDistribution( energy_cutoff, max_energy, 1.0 ) )
    response = ActiveRegion.StandardParticleResponse( response_function )
    surface_flux_estimator.setResponseFunctions( [response] )

    # Set the energy bin discretization
    surface_flux_estimator.setSourceEnergyDiscretization( energy_bins )

  ## -------------------------- Particle Tracker ---------------------------- ##

    # particle_tracker = Event.ParticleTracker( 0, 20 )

    # # Add the particle tracker to the event handler
    # event_handler.addParticleTracker( particle_tracker )

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

##---------------------------------------------------------------------------##
def restartInfiniteMediumSimulation( rendezvous_file_name,
                                     num_particles,
                                     threads,
                                     time,
                                     log_file = None,
                                     num_rendezvous = None ):

    ## Initialize the MPI session
    session = MPI.GlobalMPISession( len(sys.argv), sys.argv )

    # Suppress logging on all procs except for the master (proc=0)
    Utility.removeAllLogs()
    session.initializeLogs( 0, True )

    if session.rank() == 0:
      print "The PyFrensie path is set to: ", pyfrensie_path

    if not log_file is None:
        session.initializeLogs( log_file, 0, True )

    # Set the database path
    Collision.FilledGeometryModel.setDefaultDatabasePath( database_path )

    time_sec = time*60

    if not num_rendezvous is None:
        new_simulation_properties = MonteCarlo.SimulationGeneralProperties()
        new_simulation_properties.setNumberOfHistories( int(num_particles) )
        new_simulation_properties.setSimulationWallTime( float(time_sec) )
        new_simulation_properties.setMinNumberOfRendezvous( int(num_rendezvous) )

        factory = Manager.ParticleSimulationManagerFactory( rendezvous_file_name,
                                                            new_simulation_properties,
                                                            threads )
    else:
        factory = Manager.ParticleSimulationManagerFactory( rendezvous_file_name,
                                                            int(num_particles),
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

##----------------------------------------------------------------------------##
## --------------------- setForwardSimulationProperties --------------------- ##
##----------------------------------------------------------------------------##
def setForwardSimulationProperties( histories, time, interpolation, grid_policy, mode, method, cutoff_energy, source_energy ):

  properties = setup.setSimulationProperties( histories, time, interpolation, grid_policy, mode, method, cutoff_energy, source_energy )

  ## -------------------------- ELECTRON PROPERTIES ------------------------- ##

  # Turn off Atomic Relaxation
  properties.setAtomicRelaxationModeOff( MonteCarlo.ELECTRON )

  # Turn certain reactions off
  # properties.setElasticModeOff()
  # properties.setElectroionizationModeOff()
  # properties.setBremsstrahlungModeOff()
  # properties.setAtomicExcitationModeOff()

  return properties

##----------------------------------------------------------------------------##
## -------------------- setAdjointSimulationProperties ---------------------- ##
##----------------------------------------------------------------------------##
def setAdjointSimulationProperties( histories, time, mode, method, min_energy, max_energy ):

  properties = setup.setAdjointSimulationProperties( histories, time, mode, method, min_energy, max_energy )

  ## -------------------------- ELECTRON PROPERTIES ------------------------- ##

  # Set the critical line energies
  properties.setCriticalAdjointElectronLineEnergies( [max_energy] )

  # Turn certain reactions off
  # properties.setAdjointElasticModeOff()
  # properties.setAdjointElectroionizationModeOff()
  # properties.setAdjointBremsstrahlungModeOff()
  # properties.setAdjointAtomicExcitationModeOff()

  return properties

##----------------------------------------------------------------------------##
## ---------------------- setForwardSimulationName ---------------------------##
##----------------------------------------------------------------------------##
# Define a function for naming an electron simulation
def setForwardSimulationName( properties, energy, file_type, element ):
  extension = setup.setSimulationNameExtention( properties, file_type )
  sim_name = "forward_" + element + "_" + str(energy) + extension

  date = str(datetime.datetime.today()).split()[0]
  directory = "results/forward/" + date

  return directory + "/" + sim_name

##----------------------------------------------------------------------------##
## ---------------------- setAdjointSimulationName ---------------------------##
##----------------------------------------------------------------------------##
# Define a function for naming an electron simulation
def setAdjointSimulationName( properties, energy, element, grid_policy ):
  extension = setup.setAdjointSimulationNameExtention( properties )
  sim_name = "adjoint_" + element + "_" + str(energy) + extension

  # Add the grid policy to the name
  if grid_policy == MonteCarlo.UNIT_BASE_CORRELATED_GRID:
    sim_name += "_unit_correlated"
  elif grid_policy == MonteCarlo.UNIT_BASE_GRID:
    sim_name += "_unit_base"
  elif grid_policy == MonteCarlo.CORRELATED_GRID:
    sim_name += "_correlated"
  else:
    message = 'The grid policy ' + grid_policy + ' is currently not available!'
    raise Exception(message)

  sim_name += extension

  date = str(datetime.datetime.today()).split()[0]
  directory = "results/adjoint/" + date

  return directory + "/" + sim_name

##----------------------------------------------------------------------------##
## ----------------------- Create Results Directory ------------------------- ##
##----------------------------------------------------------------------------##
def createResultsDirectory(sim_name):

  directory = path.dirname(sim_name)

  if not path.exists(directory):
    makedirs(directory)

##----------------------------------------------------------------------------##
## ----------------------- getGridPolicyFromString -------------------------- ##
##----------------------------------------------------------------------------##
def getGridPolicyFromString(raw_grid_policy):

  # Set the bivariate Grid Policy ( UNIT_BASE_CORRELATED, CORRELATED, UNIT_BASE )
  if raw_grid_policy == "unit correlated":
    return MonteCarlo.UNIT_BASE_CORRELATED_GRID
  elif raw_grid_policy == "unit base":
    return MonteCarlo.UNIT_BASE_GRID
  elif raw_grid_policy == "correlated":
    return MonteCarlo.CORRELATED_GRID
  else:
    message = 'The grid policy ' + raw_grid_policy + ' is currently not available!'
    raise Exception(message)

##----------------------------------------------------------------------------##
## ----------------------- getElasticModeFromString ------------------------- ##
##----------------------------------------------------------------------------##
def getElasticModeFromString(raw_mode):

    # Set the elastic distribution mode ( DECOUPLED, COUPLED, HYBRID )
    if raw_mode == "decoupled":
      return MonteCarlo.DECOUPLED_DISTRIBUTION
    elif raw_mode == "coupled":
      return MonteCarlo.COUPLED_DISTRIBUTION
    elif raw_mode == "hybrid":
      return MonteCarlo.HYBRID_DISTRIBUTION
    else:
      message = 'The elastic distribution mode ' + raw_mode + ' is currently not available!'
      raise Exception(message)

##----------------------------------------------------------------------------##
## -------------------- getElasticMethodFromString -------------------------- ##
##----------------------------------------------------------------------------##
def getElasticMethodFromString(raw_method):

    # Set the elastic coupled sampling method ( TWO_D, ONE_D, MODIFIED_TWO_D )
    if raw_method == "modified 2D":
      return MonteCarlo.MODIFIED_TWO_D_UNION
    elif raw_method == "2D":
      return MonteCarlo.TWO_D_UNION
    elif raw_method == "1D":
      return MonteCarlo.ONE_D_UNION
    else:
      message = 'The elastic coupled sampling method ' + raw_method + ' is currently not available!'
      raise Exception(message)
