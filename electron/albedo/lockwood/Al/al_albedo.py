#! /usr/bin/env python
from os import path, makedirs, environ
import sys
from optparse import *

# Add the parent and grandparent directory to the path
parent_dir=path.dirname(path.dirname(path.dirname(path.abspath(__file__))))
sys.path.insert(1,parent_dir)
sys.path.insert(1,path.dirname(parent_dir))

import simulation_setup as setup
import albedo_simulation as simulation
import PyFrensie.Data as Data
import PyFrensie.MonteCarlo as MonteCarlo

pyfrensie_path =path.dirname( path.dirname(path.abspath(MonteCarlo.__file__)))
database_path = environ['DATABASE_PATH']

# Set the element
element="Al"; zaid=13000

# Set the cutoff energy
cutoff_energy=1e-4

# Set the bivariate Grid Policy (UNIT_BASE_CORRELATED, CORRELATED, UNIT_BASE)
grid_policy=MonteCarlo.UNIT_BASE_CORRELATED_GRID

# Set the elastic distribution mode ( DECOUPLED, COUPLED, HYBRID )
mode=MonteCarlo.DECOUPLED_DISTRIBUTION

# Set the elastic coupled sampling method
# ( TWO_D_UNION, ONE_D_UNION, MODIFIED_TWO_D_UNION )
method=MonteCarlo.MODIFIED_TWO_D_UNION

## ------- FORWARD OPTIONS ------- ##

# set spectrum source mode
spectrum_source=True

# set isotropic source mode
isotropic_source=True

# set atomic relaxation mode
atomic_relaxation=False

# Set the source energy (1e-4 - 1.033)
source_energy=1.033

# Set the source angle in degrees ( 0.0, 60.0 )
source_angle=0.0

# Set the bivariate interpolation (LOGLOGLOG, LINLINLIN, LINLINLOG)
interpolation=MonteCarlo.LOGLOGLOG_INTERPOLATION

# Set the data file type (ACE_EPR_FILE, Native_EPR_FILE)
file_type=Data.ElectroatomicDataProperties.Native_EPR_FILE

# Set if a refined grid should be used ( True, False )
use_refined_grid=False

geometry_path = path.dirname(path.realpath(__file__)) + "/geom.h5m"


def printSimulationName():
  # Set the simulation properties
  properties = setup.setSimulationProperties( 1, 1, interpolation, grid_policy, mode, method )

  # Print the simulation name
  sim_name = simulation.setSimulationName( properties, file_type, element, source_energy, source_angle, use_refined_grid )

  print sim_name

def printAdjointSimulationName():

  # Set the adjoint simulation properties
  properties = setup.setAdjointSimulationProperties( 1, 1, mode, method )

  sim_name = simulation.setAdjointSimulationName( properties, element, grid_policy, ionization, nudge_past_max )

  print sim_name

if __name__ == "__main__":

    # Parse the command line options
    parser = OptionParser()
    parser.add_option("--threads", type="int", dest="threads", default=1,
                      help="the number of threads to use")
    parser.add_option("--log_file", type="string", dest="log_file",
                      help="the file that will be used for logging")
    parser.add_option("--num_particles", type="float", dest="num_particles", default=1e3,
                      help="the number of particles to run")
    parser.add_option("--time", type="float", dest="time", default=1350.0,
                      help="the simultion wall time in minutes")
    parser.add_option("--transport", type="string", dest="transport", default="forward",
                      help="the simultion transport mode (forward/adjoint)")
    options,args = parser.parse_args()

    if options.transport == "forward":
      # Set the simulation properties
      properties = setup.setSimulationProperties( options.num_particles, options.time, interpolation, grid_policy, mode, method )

      # Set the min electron energy
      properties.setMinElectronEnergy( cutoff_energy )

      # Turn certain reactions off
      # properties.setElasticModeOff()
      # properties.setElectroionizationModeOff()
      # properties.setBremsstrahlungModeOff()
      # properties.setAtomicExcitationModeOff()

      if not atomic_relaxation:
        # Turn off Atomic Relaxation
        properties.setAtomicRelaxationModeOff( MonteCarlo.ELECTRON )

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

      if isotropic_source and spectrum_source:
        # Set the simulation name and title
        sim_name = simulation.setSimulationName( properties, file_type, element, "cosine_spectrum", source_angle, use_refined_grid )

        # Create the results directory
        simulation.createResultsDirectory(sim_name)

        # Run the simulation
        simulation.runForwardIsotrpoicSpectrumAlbedoSimulation(
                                                       sim_name,
                                                       database_path,
                                                       geometry_path,
                                                       properties,
                                                       cutoff_energy,
                                                       source_energy,
                                                       source_angle,
                                                       zaid,
                                                       file_type,
                                                       version,
                                                       options.threads,
                                                       options.log_file )
      elif spectrum_source:
        # Set the simulation name and title
        sim_name = simulation.setSimulationName( properties, file_type, element, "spectrum", source_angle, use_refined_grid )

        # Create the results directory
        simulation.createResultsDirectory(sim_name)

        # Run the simulation
        simulation.runForwardSpectrumAlbedoSimulation( sim_name,
                                                       database_path,
                                                       geometry_path,
                                                       properties,
                                                       cutoff_energy,
                                                       source_energy,
                                                       source_angle,
                                                       zaid,
                                                       file_type,
                                                       version,
                                                       options.threads,
                                                       options.log_file )

      else:
        # Set the simulation name and title
        sim_name = simulation.setSimulationName( properties, file_type, element, source_energy, source_angle, use_refined_grid )

        # Create the results directory
        simulation.createResultsDirectory(sim_name)

        # Run the simulation
        simulation.runForwardAlbedoSimulation( sim_name,
                                               database_path,
                                               geometry_path,
                                               properties,
                                               source_energy,
                                               source_angle,
                                               zaid,
                                               file_type,
                                               version,
                                               options.threads,
                                               options.log_file )

    elif options.transport == "adjoint":
      max_source_energy = 1.033
      # Set the adjoint simulation properties
      properties = setup.setAdjointSimulationProperties( options.num_particles, options.time, mode, method, cutoff_energy, max_source_energy )

      # Set the cutoff weight properties for rouletting
      # properties.setAdjointElectronRouletteThresholdWeight( 1e-8 )
      # properties.setAdjointElectronRouletteSurvivalWeight( 1e-6 )

      # Turn certain reactions off
      # properties.setAdjointElasticModeOff()
      # properties.setAdjointElectroionizationModeOff()
      # properties.setAdjointBremsstrahlungModeOff()
      # properties.setAdjointAtomicExcitationModeOff()

      # Create the results directory
      simulation.createAdjointResultsDirectory()

      # Set the simulation name and title
      sim_name = simulation.setAdjointSimulationName( properties, element, grid_policy, ionization, nudge_past_max )

      if grid_policy == MonteCarlo.UNIT_BASE_GRID:
        version = 3
      elif grid_policy == MonteCarlo.UNIT_BASE_CORRELATED_GRID:
        version = 4
      elif grid_policy == MonteCarlo.CORRELATED_GRID:
        version = 5

      # Run the simulation
      simulation.runAdjointAlbedoSimulation( sim_name,
                                             database_path,
                                             geometry_path,
                                             properties,
                                             cutoff_energy,
                                             max_source_energy,
                                             zaid,
                                             version,
                                             options.threads,
                                             options.log_file )
