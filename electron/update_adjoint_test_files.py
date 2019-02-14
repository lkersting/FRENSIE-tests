#!/usr/bin/python2.7
##---------------------------------------------------------------------------##
## Adjoint test data updater
##---------------------------------------------------------------------------##


import sys
from os import path, makedirs
import datetime
from optparse import *

# Add the parent directory to the path
sys.path.insert(1, path.dirname(path.dirname(path.abspath(__file__))))
import simulation_setup as setup
from native_epr_to_native_aepr import generateData, addToDatabase
import PyFrensie.Data as Data

# Simple class for color output
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

if __name__ == "__main__":

    # Parse the command-line arguments
    parser = OptionParser()
    parser.add_option("-d", "--db_name", type="string", dest="db_name", default="database.xml",
                      help="the database file with path")
    parser.add_option("-z", "--zaid", type="int", dest="zaid",
                      help="the atom number of the element")
    parser.add_option("-e", "--max_energy", type="float", dest="max_energy", default=0.01,
                      help="the max electron energy")
    parser.add_option("-g", "--grid_policy", type="string", dest="grid_policy", default="UnitBaseCorrelated",
                      help="the electron two d grid policy")
    parser.add_option("-v", "--version", type="int", dest="version", default=0,
                      help="the data file version number")
    parser.add_option("--scatter_above_max_mode_off", action="store_false", dest="above_max_mode", default=True,
                      help="Don't allow adjoint electrons to scatter above the max energy.")
    parser.add_option("-i", "--ionization_sampling_mode", type="string", dest="ionization_sampling_mode", default="Knock-On",
                      help="The forward electroionization sampling mode")
    parser.add_option("--ion_eval_tol", type="float", dest="ionization_eval_tol", default=1e-6,
                      help="The electroionization evaluation tolerance")
    parser.add_option("--ion_grid_convergence", type="float", dest="ionization_grid_convergence", default=1e-4,
                      help="The electroionization grid convergence")
    parser.add_option("--brem_eval_tol", type="float", dest="brem_eval_tol", default=1e-6,
                      help="The bremsstrahlung evaluation tolerance")
    parser.add_option("--brem_grid_convergence", type="float", dest="brem_grid_convergence", default=1e-4,
                      help="The bremsstrahlung grid convergence")
    parser.add_option("--xs_grid_convergence", type="float", dest="xs_grid_convergence", default=1e-4,
                      help="The electron cross section grid convergence")
    options,args = parser.parse_args()

    if path.exists( options.db_name ):
      # Open the database
      database = Data.ScatteringCenterPropertiesDatabase( options.db_name )
    else:
      print "ERROR: The file", options.db_name, "doesn't exist!"
      sys.exit(1)

    print "The PyFrensie path is set to: ", path.dirname( path.dirname(path.abspath(Data.__file__)))

    # Get the atom properties for the zaid
    zaid = Data.ZAID(options.zaid)
    atom_name = zaid.toName()
    if not database.doAtomPropertiesExist( zaid ):
        print "The database does not contain the " + atom_name + " data"
        sys.exit(1)

    element_properties = database.getAtomProperties( zaid )

    # Check if the native file exists
    if not element_properties.photoatomicDataAvailable( Data.PhotoatomicDataProperties.Native_EPR_FILE ):
        print "The database does not contain the " + atom_name + " native data"
        sys.exit(1)

    if options.grid_policy == 'UnitBase':
      epr_version = 1
    elif options.grid_policy == 'UnitBaseCorrelated':
      epr_version = 2
    elif options.grid_policy == 'Correlated':
      epr_version = 3

    if not element_properties.photoatomicDataAvailable( Data.PhotoatomicDataProperties.Native_EPR_FILE, epr_version ):
        print "The database does not contain version ", epr_version, " of " + atom_name + " native data"
        sys.exit(1)

    data_properties = element_properties.getSharedPhotoatomicDataProperties( Data.PhotoatomicDataProperties.Native_EPR_FILE, epr_version )

    epr_file_name = path.dirname(options.db_name) + "/" + data_properties.filePath()

    # Get/create the aepr directory path
    aepr_directory = path.dirname(path.dirname( epr_file_name )) + "/aepr"

    if not path.exists(aepr_directory):
      makedirs( aepr_directory )

    # Get the date for the table notes
    today = str(datetime.datetime.today())
    notes="This table was generated on " + today + ". It is for testing only!"

    # Update adjoint Hydrogen data
    print bcolors.BOLD + "Updating the adjoint " + atom_name + " native test data ...\n" + bcolors.ENDC

    min_photon_energy = 1e-3
    max_photon_energy = 3.0
    min_electron_energy = 1e-4

    # Set default photon grid tolerances
    photon_grid_convergence_tol = 1e-3
    photon_grid_abs_diff_tol = 1e-42
    photon_grid_dist_tol = 1e-16

    adjoint_pp_energy_dist_norm_const_eval_tol = 1e-3
    adjoint_pp_energy_dist_norm_const_nudge_val = 1e-6
    adjoint_tp_energy_dist_norm_const_eval_tol = 1e-3
    adjoint_tp_energy_dist_norm_const_nudge_val = 1e-6
    adjoint_incoherent_max_energy_nudge_val = 0.2
    adjoint_incoherent_energy_to_max_energy_nudge_val = 1e-5
    adjoint_incoherent_eval_tol = 1e-2
    adjoint_incoherent_grid_convergence_tol = 0.5
    adjoint_incoherent_grid_abs_diff_tol = 1e-42
    adjoint_incoherent_grid_dist_tol = 1e-16

    # Set default electron grid tolerances
    electron_grid_abs_diff_tol = 1e-20
    electron_grid_dist_tol = 1e-16

    cutoff_angle_cosine = 1.0
    num_moment_preserving_angles = 0
    tabular_evaluation_tol = 1e-7
    electron_two_d_interp_policy = "LogLogLog"
    brems_min_energy_nudge_val = 1e-9
    brems_max_energy_nudge_val = 1e-6
    brems_grid_abs_diff_tol = 1e-20
    brems_grid_dist_tol = 1e-16

    electroion_min_energy_nudge_val = 1e-9
    electroion_max_energy_nudge_val = 1e-6
    electroion_abs_diff_tol = 1e-20
    electroion_dist_tol = 1e-16

    # Set the aepr file name
    atomic_number = zaid.atomicNumber()
    aepr_file_name = aepr_directory + "/aepr_native_" + str(atomic_number) + "_v" + str(options.version) + ".xml"


    if options.above_max_mode:
      above_max = "on"
    else:
      above_max = "off"
    print bcolors.BOLD + "Updating file version " + str(options.version) + " with a " + options.grid_policy + " grid policy, " + options.ionization_sampling_mode +" electroionization sampling mode and scatter above max energy mode " + above_max + bcolors.ENDC

    data_container = \
    generateData( epr_file_name,
                  aepr_file_name,
                  True,
                  notes,
                  min_photon_energy,
                  max_photon_energy,
                  min_electron_energy,
                  options.max_energy,
                  photon_grid_convergence_tol,
                  photon_grid_abs_diff_tol,
                  photon_grid_dist_tol,
                  adjoint_pp_energy_dist_norm_const_eval_tol,
                  adjoint_pp_energy_dist_norm_const_nudge_val,
                  adjoint_tp_energy_dist_norm_const_eval_tol,
                  adjoint_tp_energy_dist_norm_const_nudge_val,
                  adjoint_incoherent_max_energy_nudge_val,
                  adjoint_incoherent_energy_to_max_energy_nudge_val,
                  adjoint_incoherent_eval_tol,
                  adjoint_incoherent_grid_convergence_tol,
                  adjoint_incoherent_grid_abs_diff_tol,
                  adjoint_incoherent_grid_dist_tol,
                  options.xs_grid_convergence,
                  electron_grid_abs_diff_tol,
                  electron_grid_dist_tol,
                  cutoff_angle_cosine,
                  num_moment_preserving_angles,
                  tabular_evaluation_tol,
                  options.above_max_mode,
                  electron_two_d_interp_policy,
                  options.grid_policy,
                  brems_min_energy_nudge_val,
                  brems_max_energy_nudge_val,
                  options.brem_eval_tol,
                  options.brem_grid_convergence,
                  brems_grid_abs_diff_tol,
                  brems_grid_dist_tol,
                  options.ionization_sampling_mode,
                  electroion_min_energy_nudge_val,
                  electroion_max_energy_nudge_val,
                  options.ionization_eval_tol,
                  options.ionization_grid_convergence,
                  electroion_abs_diff_tol,
                  electroion_dist_tol )
    # except Exception as e:
    #     print(bcolors.BOLD + bcolors.FAIL + '\nadjoint " + atom_name + " native data FAILED to update: '+ str(e))
    #     sys.exit(1)

    addToDatabase( aepr_file_name,
                   path.dirname( options.db_name ),
                   database,
                   data_container.getAtomicNumber(),
                   data_container.getAtomicWeight(),
                   options.version )

    database.saveToFile( options.db_name, True )

    print bcolors.BOLD + bcolors.OKGREEN + "adjoint " + atom_name + " native data updated successfully!\n" + bcolors.ENDC
    print "New file located at: ", aepr_file_name
