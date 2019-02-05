import numpy
from os import path, makedirs
import sys
import PyFrensie.Geometry as Geometry
import PyFrensie.Geometry.DagMC as DagMC
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
import PyFrensie.Data as Data
import PyFrensie.Data.Native as Native

pyfrensie_path =path.dirname( path.dirname(path.abspath(MonteCarlo.__file__)))

##---------------------------------------------------------------------------##
## Set up and run the forward simulation
def runForwardInfiniteMediumSimulation( sim_name,
                                        db_path,
                                        geom_name,
                                        num_particles,
                                        simulation_properties,
                                        source_energy,
                                        surface_ids,
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
    electron_distribution = [ActiveRegion.StandardElectronSourceComponent( 0, 1.0, model, particle_distribution )]

    # Assign the electron source component to the source
    source = ActiveRegion.StandardParticleSource( [electron_distribution] )

  ##--------------------------------------------------------------------------##
  ## -------------------------- EVENT HANDLER SETUP ------------------------- ##
  ##--------------------------------------------------------------------------##

    ## Set up the event handler
    event_handler = Event.EventHandler( model, simulation_properties )

  ## ------------------------ Surface Flux Estimator ------------------------ ##

    # Create the surface flux estimator
    surface_flux_estimator = Event.WeightMultipliedSurfaceFluxEstimator( estimator_id, 1.0, surface_ids, model )
    surface_flux_estimator.setEnergyDiscretization( energy_bins )
    surface_flux_estimator.setParticleTypes( [MonteCarlo.ELECTRON] )

    event_handler.addEstimator( surface_flux_estimator )

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
def runAdjointInfiniteMediumSimulation( sim_name,
                                        db_path,
                                        geom_name,
                                        num_particles,
                                        simulation_properties,
                                        energy_cutoff,
                                        source_energy,
                                        source_critical_line,
                                        surface_ids,
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
    particle_distribution = ActiveRegion.StandardParticleDistribution( "isotropic mono-energetic dist" )

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
    estimator_id = 1

    # Create the estimator response function
    response_function = ActiveRegion.EnergyParticleResponseFunction( Distribution.DeltaDistribution( source_energy, source_energy - energy_cutoff ) )
    response = ActiveRegion.StandardParticleResponse( response_function )

    # Create the surface flux estimator
    surface_flux_estimator = Event.WeightMultipliedSurfaceFluxEstimator( estimator_id, 1.0, surface_ids, model )
    surface_flux_estimator.setSourceEnergyDiscretization( energy_bins )
    surface_flux_estimator.setResponseFunctions( [response] )
    surface_flux_estimator.setParticleTypes( [MonteCarlo.ADJOINT_ELECTRON] )

    event_handler.addEstimator( surface_flux_estimator )

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
                                                     num_particles,
                                                     simulation_properties,
                                                     min_source_energy,
                                                     max_source_energy,
                                                     surface_ids,
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
    particle_distribution = ActiveRegion.StandardParticleDistribution( "isotropic uniform dist" )

    uniform_energy = Distribution.UniformDistribution( min_source_energy, max_source_energy )
    energy_dimension_dist = ActiveRegion.IndependentEnergyDimensionDistribution( uniform_energy )
    particle_distribution.setDimensionDistribution( energy_dimension_dist )
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

    # Create the surface flux estimator
    surface_flux_estimator = Event.WeightMultipliedSurfaceFluxEstimator( 1, 1.0, surface_ids, model )
    surface_flux_estimator.setEnergyDiscretization( energy_bins )
    surface_flux_estimator.setParticleTypes( [MonteCarlo.ELECTRON] )

    event_handler.addEstimator( surface_flux_estimator )

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
## Set up and run the adjoint simulation
def runAdjointUniformEnergyInfiniteMediumSimulation( sim_name,
                                                     db_path,
                                                     geom_name,
                                                     num_particles,
                                                     simulation_properties,
                                                     energy_cutoff,
                                                     energy_max,
                                                     source_critical_line,
                                                     surface_ids,
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
    particle_distribution = ActiveRegion.StandardParticleDistribution( "isotropic uniform dist" )

    uniform_energy = Distribution.UniformDistribution( energy_cutoff, energy_max )
    energy_dimension_dist = ActiveRegion.IndependentEnergyDimensionDistribution( uniform_energy )
    particle_distribution.setDimensionDistribution( energy_dimension_dist )
    particle_distribution.setPosition( 0.0, 0.0, 0.0 )
    particle_distribution.constructDimensionDistributionDependencyTree()

    # The generic distribution will be used to generate electrons
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

    # Create the estimator response function
    response_function = ActiveRegion.EnergyParticleResponseFunction( Distribution.UniformDistribution( energy_cutoff, energy_max, 1.0/(energy_max - energy_cutoff) ) )
    response = ActiveRegion.StandardParticleResponse( response_function )

    # Create the surface flux estimator
    surface_flux_estimator = Event.WeightMultipliedSurfaceFluxEstimator( 1, 1.0, surface_ids, model )
    surface_flux_estimator.setSourceEnergyDiscretization( energy_bins )
    surface_flux_estimator.setResponseFunctions( [response] )
    surface_flux_estimator.setParticleTypes( [MonteCarlo.ADJOINT_ELECTRON] )

    event_handler.addEstimator( surface_flux_estimator )

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
                                     db_path,
                                     num_particles,
                                     threads,
                                     log_file = None,
                                     num_rendezvous = None ):

    ## Initialize the MPI session
    session = MPI.GlobalMPISession( len(sys.argv), sys.argv )

    # Suppress logging on all procs except for the master (proc=0)
    Utility.removeAllLogs()
    session.initializeLogs( 0, True )

    if not log_file is None:
        session.initializeLogs( log_file, 0, True )

    # Set the database path
    Collision.FilledGeometryModel.setDefaultDatabasePath( db_path )

    if not num_rendevous is None:
        new_simulation_properties = MonteCarlo.SimulationGeneralProperties()
        new_simulation_properties.setNumberOfHistories( int(num_particles) )
        new_simulation_properties.setMinNumberOfRendezvous( int(num_rendezvous) )

        factory = Manager.ParticleSimulationManagerFactory( rendezvous_file_name,
                                                            new_simulation_properties,
                                                            threads )
    else:
        factory = Manger.ParticleSimulationManagerFactory( rendezvous_file_name,
                                                           int(num_particles),
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
##------------------------ processAdjointDataFromRendezvous -------------------------##
##----------------------------------------------------------------------------##

# This function pulls data from the rendezvous file
def processAdjointDataFromRendezvous( rendezvous_file ):

  Collision.FilledGeometryModel.setDefaultDatabasePath( database_path )

  # Load data from file
  manager = Manager.ParticleSimulationManagerFactory( rendezvous_file ).getManager()
  event_handler = manager.getEventHandler()

  # Get the simulation name and title
  properties = manager.getSimulationProperties()

  filename = setSimulationName( properties )
  title = setup.getSimulationPlotTitle( filename )

  print "Processing the results:"
  processAdjointData( event_handler, filename, title )

  print "Results will be in ", path.dirname(filename)

##----------------------------------------------------------------------------##
##------------------------------- processData --------------------------------##
##----------------------------------------------------------------------------##
def processAdjointData( event_handler, filename, title ):

  # Process surface flux data
  surface_flux = event_handler.getEstimator( 2 )

  file1 = filename + "_1"
  setup.processSurfaceFluxSourceEnergyBinData( surface_flux, 1, file1, title )

  file2 = filename + "_2"
  setup.processSurfaceFluxSourceEnergyBinData( surface_flux, 18, file2, title )

  file3 = filename + "_5"
  setup.processSurfaceFluxSourceEnergyBinData( surface_flux, 16, file3, title )