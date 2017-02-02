""" mfdn_v14_b06.py -- MFDn version 14 beta06 wrappers for mcscript

  Created by M. A. Caprio, University of Notre Dame.
  8/2/14 (mac): Originated, based on old wrappers in ncsm.py.
  5/13/15 (mac): Insert "future" statements for Python 2 legacy support.
  Last modified 5/22/15 (mac).

"""
from __future__ import print_function, division

import os
import sys
import glob

import mcscript
import ncsm_config

################################################################
# mfdn-h2 diagonalization run
################################################################

# (05/28/15) V. Constantinou. Adding functionality for full configuration interaction runs..

def call_mfdn_h2 (current_task):
    """ Invoke mfdn for h2 run.

    Takes full task dictionary for h2 run, plus extra parameters:

        Nshell : one-body cutoff (1-based)
        h2_basename : base file name for Hamiltonian file
        obs_basename_list : list of base file names for two-body observables

    """

    # prepare mfdn input file
    mfdn_input_filename = "mfdn.dat"
    # note: emulating missing commas after floats from Pieter's input -- not sure if necessary

    # base parameters
    # set hw value for transition operator oscillator length
    #    value 0 disables transition calculation, e.g., for non-oscillator bases
    hw_for_trans = mcscript.ifelse(current_task["traditional_ho"],current_task["hw"],0)

    # FCI run.. Nmax should be equal to (Z+N)*(Nshell-1) because the many-body basis in a full configuration
    # has all the particles in the Nshell level.. mfdn uses a zero-based system to cound the total number of quanta
    # in the line
    # Nmix, Nmax, Nstep
    # found in mfdn.dat. For a full configuration calculation we must set Nmax equal to the total number of quanta
    # in the configuration with all the particles populating the Nshell level..

    if(current_task.get("fci",False)):
        Nmax = (current_task["nuclide"][0]+current_task["nuclide"][1])*(current_task["Nshell"]-1)
    else:
        Nmax = current_task["Nmax"]

    mfdn_base_parameters = [
        "%d" % ( 0, ), #IFLAGMBSI
        "%d" % ( 0, ), # ndiag (0: no spares, automatic ndiag)
        "%d" % ( 2, ),  # number of classes
        "%d" % ( current_task["nuclide"][0], ),  # protons (class 1 particles)
        "%d" % ( current_task["nuclide"][1], ),  # neutrons (class 2 particles)
        "%d, %d" % ( 1, current_task["Nshell"]),  # min, max # S.P. shells for class 1 particles
        "%d, %d" % ( 1, current_task["Nshell"]),  # min, max # S.P. shells for class 2 particles
        "%d, %d, %d" % ( current_task["Nmin"], Nmax, current_task["Nstep"] ),  # N_min, N_max, delta_N 
        "%d" % ( current_task["2mj"], ), # Total 2 M_j
        "%d, %d, %d, %d" % ( current_task["eigenvectors"], current_task["lanczos"], current_task["initial_vector"], 0 ),  # number of eigenvalues/vectors, max number of its, starting number of its
        "%d, %d" % ( 2, 2),  # rank of input Hamiltonian/interaction
        "%.3f, %s" % (hw_for_trans,"938.92") # h-bar*omega, Nucleon mass (MeV) 
        ]

    # tbo parameters
    mfdn_tbo_parameters = [
        "%s" % ( current_task["h2_basename"],), 
        "%d" % ( 2 + len(current_task["obs_basename_list"]), ) # number of observables (J, T, R2, ...)
        ]
    mfdn_tbo_parameters += current_task["obs_basename_list"]
        
    # obdme parameters
    twice_multipolarity = current_task["obdme_twice_multipolarity"]
    twice_max_delta_J = twice_multipolarity  # hard-coded choice for this script
    if (current_task["obdme_reference_state_list"] == "all2all"):
        num_reference_states = -1
        reference_state_list = []
    else:
        num_reference_states = len(current_task["obdme_reference_state_list"])
        reference_state_list = current_task["obdme_reference_state_list"]
    mfdn_obdme_parameters = [
        "%d, %d" % (1, twice_multipolarity), # static one-body density matrix elements (0: no one-body densities), twice multipolarity
        "%d, %d" % ( num_reference_states, twice_max_delta_J) #number of reference states for transitions (0: no transitions, -1: all2all), max delta2J (?)
        ] + [
        "%d, %d, %d" % reference_state
        for reference_state in reference_state_list
        ]
    mfdn_input = mfdn_base_parameters + mfdn_tbo_parameters + mfdn_obdme_parameters
    
    mcscript.write_input(mfdn_input_filename,mfdn_input)

    # execute mfdn
    # assumes OpenMP parallel environment set up previously
    executable = os.path.join(ncsm_config.exec_dir_mfdn, current_task["mfdn_executable"])
    mcscript.call_parallel([executable], mcscript.run) 

    # test for success
    if (not os.path.exists("mfdn.res")):
        print ("FAILURE: mfdn.res not found.")
        sys.exit(1)


################################################################
# mfdn-tbme conversion run
################################################################

def call_mfdn_tbme_conv (current_task):
    """ Invoke mfdn for tbme conversion run.

    File format is termed 'EJ/PN format' in mfdn.
    
    Conversion is obtained by:
    
    - setting IFLAGMBSI = 1 for construct basis (and convert interaction if rank of input H > 1)

    - recall at least 2 nucleons of each type recommended for conversion, since 2-body

    - Nmax should be N2max (else verified that interaction file is truncated)

    - tbme format has extra line relative to vxx format

      EX: 3, 0, 10, 10 nset, nskip, Ntot1max, Ntot2max (for TBME input file)

      It would appear from the names that these are 0-based cutoffs,
      so we will want Ntot1max = Ntot2max = <IPCUT>.

    - either vxx or tbme has extra line relative to h2 format:

      EX: 1.0, 5.00, 0.0, 1.0	Trel, Hcm, Vcoul, V_NN input strength

    Takes full task dictionary for tbme conversion run, plus extra parameters:
        
        N1bmax : one-body cutoff (0-based)
        N2bmax : two-body cutoff (0-based)

    Files: may require "spectator" tbme input files for the usual Hamiltonian components??

    """

    # prepare mfdn input file
    mfdn_input_filename = "mfdn.dat"
    # note: emulating missing commas after floats from Pieter's input -- not sure if necessary

    Nshell = current_task["N1bmax"] + 1 # one-body cutoff (1-based)

    # set hw value for transition operator oscillator length
    #    value 0 disables transition calculation, e.g., for non-oscillator bases
    mfdn_base_parameters = [
        "%d" % ( 1, ), #IFLAGMBSI
        "%d" % ( 0, ), # ndiag (0: no spares, automatic ndiag)
        "%d" % ( 2, ),  # number of classes
        "%d" % ( 2, ),  # protons (class 1 particles)
        "%d" % ( 2, ),  # neutrons (class 2 particles)
        "%d, %d" % ( 1, Nshell),  # min, max # S.P. shells for class 1 particles
        "%d, %d" % ( 1, Nshell),  # min, max # S.P. shells for class 2 particles
        "%d, %d, %d" % ( 0, current_task["N2bmax"], 2 ),  # N_min, N_max, delta_N 
        "%d" % ( 0, ), # Total 2 M_j
        "%d, %d, %d, %d" % ( 0, 0, -3, 0 ),  # number of eigenvalues/vectors, max number of its, starting number of its
        "%d, %d" % ( 2, 2),  # rank of input Hamiltonian/interaction
        "%.3f, %s" % (current_task["hw"], "938.92"), # h-bar*omega, Nucleon mass (MeV)
        "%d %d %d %d" % (3, 0, current_task["N1bmax"], current_task["N2bmax"]), # nset, nskip, Ntot1max, Ntot2max (for TBME input file)
        "%.3f %.3f %.3f  %.3f" % (0., 0., 0., 1.),  # Trel, Hcm, Vcoul, V_NN input strength
        "%s" % ( current_task["int_basename"],)
        ]
    mfdn_input = mfdn_base_parameters
    
    mcscript.write_input(mfdn_input_filename,mfdn_input)

    # execute mfdn
    # assumes OpenMP parallel environment set up previously
    executable = os.path.join(ncsm_config.mfdn_dir, current_task["mfdn_executable"])
    mcscript.call_parallel([executable], mcscript.run) 

    # test for success
    if (not os.path.exists("mfdn.res")):
        print ("FAILURE: mfdn.res not found.")
        sys.exit(1)
