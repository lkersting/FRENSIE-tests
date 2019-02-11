#! /usr/bin/env python
from os import path
import sys
import numpy as np
import datetime
import getpass

frensie_install=''
# Set frensie install for the lkersting (always the same directory as frensie-tests)
if getpass.getuser() == 'lkersting':
  frensie_install = path.dirname(path.dirname(path.dirname(path.abspath(__file__))))
  sys.path.insert(1, frensie_install + '/bin/')
  sys.path.insert(1, frensie_install + '/lib/python2.7/site-packages/')

# NOTE: If a specific version of FRENSIE is desired, the path below can be
# uncommented and the desired path to the frensie/lib can be used.
# frensie_install = path.dirname(path.dirname(path.dirname(path.abspath(__file__))))
# sys.path.insert(1, frensie_install + '/bin/')
# sys.path.insert(1, frensie_install + '/lib/python2.7/site-packages/')

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

##----------------------------------------------------------------------------##
## ------------------------- SIMULATION PROPERTIES -------------------------- ##
##----------------------------------------------------------------------------##
def setSimulationProperties( histories, time, interpolation, grid_policy, elastic_mode, elastic_sampling_method ):

  properties = MonteCarlo.SimulationProperties()

  ## -------------------------- GENERAL PROPERTIES -------------------------- ##

  # Set the particle mode
  properties.setParticleMode( MonteCarlo.ELECTRON_MODE )

  # Set the number of histories
  properties.setNumberOfHistories( histories )

  # Set the minimum number of rendezvous
  if histories > 100:
    properties.setMinNumberOfRendezvous( 10 )

  # Change time from minutes to seconds
  time_sec = time*60

  # Set the wall time
  properties.setSimulationWallTime( time_sec )

  ## -------------------------- NEUTRON PROPERTIES -------------------------- ##

  ## -------------------------- PHOTON PROPERTIES --------------------------- ##

  ## ------------------------- ELECTRON PROPERTIES -------------------------- ##

  # Set the min electron energy in MeV (Default is 100 eV)
  properties.setMinElectronEnergy( 1e-4 )

  # Set the max electron energy in MeV (Default is 20 MeV)
  properties.setMaxElectronEnergy( 20.0 )

  # Set the bivariate interpolation (LOGLOGLOG, LINLINLIN, LINLINLOG)
  properties.setElectronTwoDInterpPolicy( interpolation )

  # Set the bivariate Grid Policy (UNIT_BASE_CORRELATED, CORRELATED, UNIT_BASE)
  properties.setElectronTwoDGridPolicy( grid_policy )

  # Set the electron evaluation tolerance (Default is 1e-6)
  properties.setElectronEvaluationTolerance( 1e-6 )

  ## --- Elastic Properties ---

  # Turn elastic electron scattering off
  # properties.setElasticModeOff()

  # Set the elastic distribution mode ( DECOUPLED, COUPLED, HYBRID )
  properties.setElasticElectronDistributionMode( elastic_mode )

  # Set the elastic coupled sampling method
  # ( TWO_D_UNION, ONE_D_UNION, MODIFIED_TWO_D_UNION )
  properties.setCoupledElasticSamplingMode( elastic_sampling_method )

  # Set the elastic cutoff angle cosine ( -1.0 < mu < 1.0 )
  properties.setElasticCutoffAngleCosine( 1.0 )

  ## --- Electroionization Properties ---

  # Turn the electro-ionization reaction off
  # properties.setElectroionizationModeOff()

  ## --- Bremsstrahlung Properties ---

  # Turn electron bremsstrahlung reaction off
  # properties.setBremsstrahlungModeOff()

  # Set the bremsstrahlung angular distribution function
  # ( DIPOLE, TABULAR, ??? )
  properties.setBremsstrahlungAngularDistributionFunction( MonteCarlo.DIPOLE_DISTRIBUTION )

  ## --- Atomic Excitation Properties ---

  # Turn electron atomic excitation reaction off
  # properties.setAtomicExcitationModeOff()

  return properties

##----------------------------------------------------------------------------##
## ------------------------- SIMULATION PROPERTIES -------------------------- ##
##----------------------------------------------------------------------------##
def setAdjointSimulationProperties( histories, time, elastic_mode, elastic_sampling_method ):

  properties = MonteCarlo.SimulationProperties()

  ## -------------------------- GENERAL PROPERTIES -------------------------- ##

  # Set the particle mode
  properties.setParticleMode( MonteCarlo.ADJOINT_ELECTRON_MODE )

  # Set the number of histories
  properties.setNumberOfHistories( histories )

  # Set the minimum number of rendezvous
  if histories > 100:
    properties.setMinNumberOfRendezvous( 10 )

  # Change time from minutes to seconds
  time_sec = time*60

  # Set the wall time
  properties.setSimulationWallTime( time_sec )

  ## ---------------------- ADJOINT NEUTRON PROPERTIES ---------------------- ##

  ## ---------------------- ADJOINT PHOTON PROPERTIES ----------------------- ##

  ## --------------------- ADJOINT ELECTRON PROPERTIES ---------------------- ##

  # Set the min electron energy in MeV (Default is 100 eV)
  properties.setMinAdjointElectronEnergy( 1e-4 )

  # Set the max electron energy in MeV (Default is 20 MeV)
  properties.setMaxAdjointElectronEnergy( 20.0 )

  # Set the electron evaluation tolerance (Default is 1e-6)
  properties.setAdjointElectronEvaluationTolerance( 1e-6 )

  ## --- Adjoint Elastic Properties ---

  # Set the elastic distribution mode ( DECOUPLED, COUPLED, HYBRID )
  properties.setAdjointElasticElectronDistributionMode( elastic_mode )

  # Set the elastic coupled sampling method
  # ( TWO_D_UNION, ONE_D_UNION, MODIFIED_TWO_D_UNION )
  properties.setAdjointCoupledElasticSamplingMode( elastic_sampling_method )

  # Set the elastic cutoff angle cosine ( -1.0 < mu < 1.0 )
  properties.setAdjointElasticCutoffAngleCosine( 1.0 )

  return properties

##----------------------------------------------------------------------------##
## ---------------------- setSimulationNameExtention -------------------------##
##----------------------------------------------------------------------------##
# Define a function for naming an electron simulation
def setSimulationNameExtention( properties, file_type ):

  if file_type == Data.ElectroatomicDataProperties.ACE_EPR_FILE:
    # Use ACE EPR14 data
    name = "epr14"
  else:
    # Use Native analog data
    name = ""

  # Set the interp
  if properties.getElectronTwoDInterpPolicy() == MonteCarlo.LOGLOGLOG_INTERPOLATION:
      interp = "loglog"
  elif properties.getElectronTwoDInterpPolicy() == MonteCarlo.LINLINLIN_INTERPOLATION:
      interp = "linlin"
  else:
      interp = "linlog"

  # Set the sampling name
  sample_name=""
  if properties.getElectronTwoDGridPolicy() == MonteCarlo.UNIT_BASE_CORRELATED_GRID:
      sample_name = "unit_correlated"
  elif properties.getElectronTwoDGridPolicy() == MonteCarlo.CORRELATED_GRID:
      sample_name = "correlated"
  else:
      sample_name = "unit_base"

  # Set the name reaction and extention
  name_extention = ""
  name_reaction = ""
  if properties.isElasticModeOn():
    if properties.getElasticElectronDistributionMode() == MonteCarlo.COUPLED_DISTRIBUTION:
      if properties.getCoupledElasticSamplingMode() == MonteCarlo.MODIFIED_TWO_D_UNION:
        name_extention += "_m2d"
      elif properties.getCoupledElasticSamplingMode() == MonteCarlo.TWO_D_UNION:
        name_extention += "_2d"
      else:
        name_extention += "_1d"
    elif properties.getElasticElectronDistributionMode() == MonteCarlo.DECOUPLED_DISTRIBUTION:
      name_extention += "_decoupled"
    elif properties.getElasticElectronDistributionMode() == MonteCarlo.HYBRID_DISTRIBUTION:
      name_extention += "_hybrid"
  else:
    name_reaction = name_reaction + "_no_elastic"

  if not properties.isBremsstrahlungModeOn():
    name_reaction += "_no_brem"

  if not properties.isAtomicExcitationModeOn():
      name_reaction += "_no_excitation"

  if not properties.isElectroionizationModeOn():
      name_reaction += "_no_ionization"
  elif properties.getElectroionizationSamplingMode() == MonteCarlo.OUTGOING_ENERGY_SAMPLING:
    name_reaction += "_outgoing_energy"

  date = str(datetime.datetime.today()).split()[0]
  if name == "epr14":
    name = "_" + name + name_reaction
  else:
    name = "_" + interp + "_" + sample_name + name_extention + name_reaction

  return name

##----------------------------------------------------------------------------##
## ------------------ setAdjointSimulationNameExtention ----------------------##
##----------------------------------------------------------------------------##
# Define a function for naming an electron simulation
def setAdjointSimulationNameExtention( properties ):

  # Set the name reaction and extention
  name_extention = ""
  name_reaction = ""
  if properties.isAdjointElasticModeOn():
    if properties.getAdjointElasticElectronDistributionMode() == MonteCarlo.COUPLED_DISTRIBUTION:
      if properties.getAdjointCoupledElasticSamplingMode() == MonteCarlo.MODIFIED_TWO_D_UNION:
        name_extention += "_m2d"
      elif properties.getAdjointCoupledElasticSamplingMode() == MonteCarlo.TWO_D_UNION:
        name_extention += "_2d"
      else:
        name_extention += "_1d"
    elif properties.getAdjointElasticElectronDistributionMode() == MonteCarlo.DECOUPLED_DISTRIBUTION:
      name_extention += "_decoupled"
    elif properties.getAdjointElasticElectronDistributionMode() == MonteCarlo.HYBRID_DISTRIBUTION:
      name_extention += "_hybrid"
  else:
    name_reaction = name_reaction + "_no_elastic"

  if not properties.isAdjointBremsstrahlungModeOn():
    name_reaction += "_no_brem"
  if not properties.isAdjointElectroionizationModeOn():
      name_reaction += "_no_ionization"
  if not properties.isAdjointAtomicExcitationModeOn():
      name_reaction += "_no_excitation"

  date = str(datetime.datetime.today()).split()[0]
  name = name_extention + name_reaction

  return name

##----------------------------------------------------------------------------##
## ------------------------ getSimulationPlotTitle ---------------------------##
##----------------------------------------------------------------------------##
# Define a function for creating a plot title for an electron simulation
def getSimulationPlotTitle( filename ):

  if "epr14" in filename:
    title = "FRENSIE-ACE"
  else:
    # Set the interp in title
    title = ""

    # Set the interp in title
    if "loglog" in filename:
        title = "Log-Log"
    elif "linlin" in filename:
        title = "Lin-Lin"
    elif "linlog" in filename:
        title = "Lin-Log"
    else:
      message = 'The filename ' + filename + ' does not include an interp type!'
      raise Exception(message)

    # Add a space
    title += " "

    # Set the sampling routine in title
    if "unit_correlated" in filename:
        title += "Unit Base Correlated"
    elif "correlated" in filename:
        title += "Correlated"
    elif "unit_base" in filename:
        title += "Unit Base"
    else:
      message = 'The filename ' + filename + ' does not include a sampling routine type!'
      raise Exception(message)

    # Add a space
    title += " "

    # Set the elastic reaction in title
    if not "_no_elastic" in filename:
      if "_m2d" in filename:
        title += "M2D"
      elif "_2d" in filename:
        title += "2D"
      elif "_1d" in filename:
        title += "1D"
      elif "_decoupled" in filename:
        title += "DE"
      elif "_hybrid" in filename:
        title += "HE"
      else:
        message = 'The filename ' + filename + ' does not include an elastic reaction type!'
        raise Exception(message)

  return title

##----------------------------------------------------------------------------##
## ------------------------ getSimulationPlotTitle ---------------------------##
##----------------------------------------------------------------------------##
# Define a function for creating a plot title for an electron simulation
def getAdjointSimulationPlotTitle( filename ):

  # Set the interp in title
  title = ""

  # Set the sampling routine in title
  if "unit_correlated" in filename:
      title += "Unit Base Correlated"
  elif "correlated" in filename:
      title += "Correlated"
  elif "unit_base" in filename:
      title += "Unit Base"
  else:
    message = 'The filename ' + filename + ' does not include a sampling routine type!'
    raise Exception(message)

  # Add a space
  title += " "

  # Set the elastic reaction in title
  if not "_no_elastic" in filename:
    if "_m2d" in filename:
      title += "M2D"
    elif "_2d" in filename:
      title += "2D"
    elif "_1d" in filename:
      title += "1D"
    elif "_decoupled" in filename:
      title += "DE"
    elif "_hybrid" in filename:
      title += "HE"
    else:
      message = 'The filename ' + filename + ' does not include an elastic reaction type!'
      raise Exception(message)

  # Set the ionization sampling in title
  if "_outgoing_energy" in filename:
    title += " Outgoing Energy Ionization Sampling"

  return title

##----------------------------------------------------------------------------##
## ------------------------ Create Results Directory ------------------------ ##
##----------------------------------------------------------------------------##
def getResultsDirectory(file_type, interpolation):

  if file_type == Data.ElectroatomicDataProperties.ACE_EPR_FILE:
    # Use ACE EPR14 data
    name = "epr14"
  else:
    # Use Native analog data
    name = ""

  # Set the interp in results directory
  if interpolation == MonteCarlo.LOGLOGLOG_INTERPOLATION:
      interp = "loglog"
  elif interpolation == MonteCarlo.LINLINLIN_INTERPOLATION:
      interp = "linlin"
  else:
      interp = "linlog"

  date = str(datetime.datetime.today()).split()[0]
  if name == "epr14":
    directory = "results/" + name + "/" + date
  else:
    directory = "results/" + interp + "/" + date

  return directory

##----------------------------------------------------------------------------##
## ------------------------ Create Results Directory ------------------------ ##
##----------------------------------------------------------------------------##
def getResultsDirectoryFromString(file_type, interpolation):

  if file_type == "ACE":
    # Use ACE EPR14 data
    name = "epr14"
  else:
    # Use Native analog data
    name = ""

  # Set the interp in results directory
  if interpolation == "LOGLOGLOG":
      interp = "loglog"
  elif interpolation == "LINLINLIN":
      interp = "linlin"
  else:
      interp = "linlog"

  date = str(datetime.datetime.today()).split()[0]
  if name == "epr14":
    directory = "results/" + name + "/" + date
  else:
    directory = "results/" + interp + "/" + date

  return directory

##----------------------------------------------------------------------------##
##---------------------- processTrackFluxEnergyBinData -----------------------##
##----------------------------------------------------------------------------##
def processTrackFluxEnergyBinData( estimator, est_id, filename, title ):

  processed_data = estimator.getEntityBinProcessedData( est_id )
  flux = processed_data['mean']
  rel_error = processed_data['re']
  energy_bins = estimator.getEnergyDiscretization()

  today = datetime.date.today()

  # Write the flux data to a file
  name = filename+"_track_flux.txt"
  out_file = open(name, 'w')

  # Write title to file
  out_file.write( "# " + title +"\n")

  # Write the header to the file
  header = "# Energy (MeV)\tTrack Flux (#/cm$^2$)\tError\t"+str(today)+"\n"
  out_file.write(header)

  # Insert a zero flux for below the firest bin boundary
  flux = np.insert( flux, 0, 0.0)
  rel_error = np.insert( rel_error, 0, 0.0)

  for i in range(0, len(flux)):
    data = str(energy_bins[i]) + '\t' + str(flux[i]) + '\t' + str(rel_error[i]) + '\n'
    out_file.write(data)
  out_file.close()

##----------------------------------------------------------------------------##
##--------------------- processSurfaceFluxEnergyBinData ----------------------##
##----------------------------------------------------------------------------##
def processSurfaceFluxEnergyBinData( estimator, est_id, filename, title ):

  processed_data = estimator.getEntityBinProcessedData( est_id )
  flux = processed_data['mean']
  rel_error = processed_data['re']
  energy_bins = estimator.getEnergyDiscretization()

  today = datetime.date.today()

  # Write the flux data to a file
  name = filename+"_flux.txt"
  out_file = open(name, 'w')

  # Write title to file
  out_file.write( "# " + title +"\n")

  # Write the header to the file
  header = "# Energy (MeV)\tSurface Flux (#/cm$^2$)\tError\t"+str(today)+"\n"
  out_file.write(header)

  # Insert a zero flux for below the firest bin boundary
  flux = np.insert( flux, 0, 0.0)
  rel_error = np.insert( rel_error, 0, 0.0)

  for i in range(0, len(flux)):
    data = str(energy_bins[i]) + '\t' + str(flux[i]) + '\t' + str(rel_error[i]) + '\n'
    out_file.write(data)
  out_file.close()

##----------------------------------------------------------------------------##
##-------------------- processSurfaceCurrentEnergyBinData --------------------##
##----------------------------------------------------------------------------##
def processSurfaceCurrentEnergyBinData( estimator, est_id, filename, title ):

  processed_data = estimator.getEntityBinProcessedData( est_id )
  current = processed_data['mean']
  rel_error = processed_data['re']
  energy_bins = estimator.getEnergyDiscretization()

  today = datetime.date.today()

  # Write the current data to a file
  name = filename+"_current.txt"
  out_file = open(name, 'w')

  # Write title to file
  out_file.write( "# " + title +"\n")

  # Write the header to the file
  header = "# Energy (MeV)\tSurface Current (#)\tError\t"+str(today)+"\n"
  out_file.write(header)

  # Insert a zero current for below the firest bin boundary
  current = np.insert( current, 0, 0.0)
  rel_error = np.insert( rel_error, 0, 0.0)

  for i in range(0, len(current)):
    data = str(energy_bins[i]) + '\t' + str(current[i]) + '\t' + str(rel_error[i]) + '\n'
    out_file.write(data)
  out_file.close()

##----------------------------------------------------------------------------##
##-------------------- processSurfaceCurrentCosineBinData --------------------##
##----------------------------------------------------------------------------##
def processSurfaceCurrentCosineBinData( estimator, est_id, filename, title ):

  processed_data = estimator.getEntityBinProcessedData( est_id )
  current = processed_data['mean']
  rel_error = processed_data['re']
  cosine_bins = estimator.getCosineDiscretization()

  today = datetime.date.today()

  # Write the current data to a file
  name = filename+"_current.txt"
  out_file = open(name, 'w')

  # Write title to file
  out_file.write( "# " + title +"\n")

  # Write the header to the file
  header = "# Cosine \tSurface Current (#)\tError\t"+str(today)+"\n"
  out_file.write(header)

  # Insert a zero current for below the firest bin boundary
  current= np.insert( current, 0, 0.0)
  rel_error = np.insert( rel_error, 0, 0.0)

  for i in range(0, len(current)):
    data = str(cosine_bins[i]) + '\t' + str(current[i]) + '\t' + str(rel_error[i]) + '\n'
    out_file.write(data)
  out_file.close()

##----------------------------------------------------------------------------##
##------------------- processTrackFluxSourceEnergyBinData --------------------##
##----------------------------------------------------------------------------##
def processTrackFluxSourceEnergyBinData( estimator, est_id, filename, title ):

  processed_data = estimator.getEntityBinProcessedData( est_id )
  flux = processed_data['mean']
  rel_error = processed_data['re']
  energy_bins = estimator.getSourceEnergyDiscretization()

  today = datetime.date.today()

  # Write the flux data to a file
  name = filename+"_source_track_flux.txt"
  out_file = open(name, 'w')

  # Write title to file
  out_file.write( "# " + title +"\n")

  # Write the header to the file
  header = "# Source Energy (MeV)\tTrack Flux (#/cm$^2$)\tError\t"+str(today)+"\n"
  out_file.write(header)

  # Insert a zero flux for below the firest bin boundary
  flux = np.insert( flux, 0, 0.0)
  rel_error = np.insert( rel_error, 0, 0.0)

  for i in range(0, len(flux)):
    data = str(energy_bins[i]) + '\t' + str(flux[i]) + '\t' + str(rel_error[i]) + '\n'
    out_file.write(data)
  out_file.close()

##----------------------------------------------------------------------------##
##------------------ processSurfaceFluxSourceEnergyBinData -------------------##
##----------------------------------------------------------------------------##
def processSurfaceFluxSourceEnergyBinData( estimator, est_id, filename, title ):

  processed_data = estimator.getEntityBinProcessedData( est_id )
  flux = processed_data['mean']
  rel_error = processed_data['re']
  energy_bins = estimator.getSourceEnergyDiscretization()

  today = datetime.date.today()

  # Write the flux data to a file
  name = filename+"_source_flux.txt"
  out_file = open(name, 'w')

  # Write title to file
  out_file.write( "# " + title +"\n")

  # Write the header to the file
  header = "# Source Energy (MeV)\tSurface Flux (#/cm$^2$)\tError\t"+str(today)+"\n"
  out_file.write(header)

  # Insert a zero flux for below the firest bin boundary
  flux = np.insert( flux, 0, 0.0)
  rel_error = np.insert( rel_error, 0, 0.0)

  # Write data to file
  for i in range(0, len(flux)):
    data = str(energy_bins[i]) + '\t' + str(flux[i]) + '\t' + str(rel_error[i]) + '\n'
    out_file.write(data)
  out_file.close()

##----------------------------------------------------------------------------##
##---------------- processSurfaceCurrentSourceEnergyBinData ------------------##
##----------------------------------------------------------------------------##
def processSurfaceCurrentSourceEnergyBinData( estimator, est_id, filename, title ):

  processed_data = estimator.getEntityBinProcessedData( est_id )
  current = processed_data['mean']
  rel_error = processed_data['re']
  energy_bins = estimator.getSourceEnergyDiscretization()

  today = datetime.date.today()

  # Write the current data to a file
  name = filename+"_source_current.txt"
  out_file = open(name, 'w')

  # Write title to file
  out_file.write( "# " + title +"\n")

  # Write the header to the file
  header = "# Source Energy (MeV)\tSurface Current (#)\tError\t"+str(today)+"\n"
  out_file.write(header)

  # Insert a zero current for below the firest bin boundary
  current = np.insert( current, 0, 0.0)
  rel_error = np.insert( rel_error, 0, 0.0)

  for i in range(0, len(current)):
    data = str(energy_bins[i]) + '\t' + str(current[i]) + '\t' + str(rel_error[i]) + '\n'
    out_file.write(data)
  out_file.close()

##----------------------------------------------------------------------------##
##------------------------ printParticleTrackInfo -------------------------##
##----------------------------------------------------------------------------##

# This function print the particle tracker info
def printParticleTrackInfo( particle_tracker ):

  history_map = particle_tracker.getHistoryData()
  print "len(history_map) = ", len(history_map)
  # print "history_map = ", history_map
  # print "len(history_map[0]) = ", len(list(history_map[0]))

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
          print "state:\tenergy\t\tweight\t\tlocation\t\t\t\tdirection\t\t\t\tcollision number"
          cached_particle_state = map_i[j][k]

          for l in range(len(cached_particle_state)):
            location = list(cached_particle_state[l][0])
            direction = list(cached_particle_state[l][1])
            energy = cached_particle_state[l][2]
            time = cached_particle_state[l][3]
            weight = cached_particle_state[l][4]
            collision = cached_particle_state[l][5]
            if collision == 0 and l > 1:
              print ''
            if l < 6:
              print l,":\t",'%.8e' % energy,"\t",'%.6e' % weight,"\t",location,"\t",direction,"\t",collision
