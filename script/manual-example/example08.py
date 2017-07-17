"""example**.py

    Mark A. Caprio
    University of Notre Dame

    + 11/24/16 (mac).

"""

import os

import generate_input

# h2mixer parameter template
params = {
    # stored radial
    "Nmax_orb" : 2,  # for target space
    "Nmax_orb_int" : 20,  # for overlaps from interaction tbmes
    "orbitals_filename" : "orbitals.dat",  # for target space
    "orbitals_int_filename" : "orbitals-int.dat",  # for overlaps from interaction tbmes
    "radial_me_filename_pattern" : "radial-me-{}{}.dat",  # "{}{}" will be replaced by {"r1","r2","k1","k2"}
    "radial_olap_filename" : "radial-olap.dat",
    # stored input TBMEs
    "VNN_filename" : os.path.join(generate_input.data_dir_h2,"run0164-ob-9/JISP16-ob-9-20.bin"),
    "VC_filename" : os.path.join(generate_input.data_dir_h2,"run0164-ob-9/VC-ob-9-20.bin"),
    # scaling/transformation/basis parameters
    "hw" : 20,
    "target_truncation" : ("tb",2),
    "hw_int" : 20,
    "hw_c" : 20,
    "xform_enabled" : False,
    "xform_truncation" : ("ob",9),
    # Hamiltonian parameters
    "A" : 4,   # atomic mass number of nucleus (for kinematic operators)
    "a" : 20,  # coefficient on Ncm in Lawson term
    "hw_cm" : None,  # hw for Ncm in Lawson term; if None, defaults to "hw"
    "use_coulomb" : True,
    # output file format
    "h2_format" : 15099  # h2 formats: 0, 15099
}

# MFdn parameter template
mfdn_params = {
    "ndiag" : 1,
    "nuclide" : (2,2),
    "Nshell" : params["Nmax_orb"]+1,
    "Nmin" : 0,
    "Nmax" : 2,
    "Nstep" : 2,
    "Mj" : 0,
    "eigenvectors" : 1,
    "lanczos" : 50,
    "initial_vector" : -2,
    "tolerance" : 1e-6,
    "hw_for_trans" : 20,
    "obs_basename_list" : ["tbme-rrel2","tbme-Ncm"],
    "obdme_multipolarity" : 2,
    "obdme_reference_state_list" : [(0,0,1)]
}

with open("run.csh","w") as os:
    os.write(generate_input.run_script(params))
with open("h2mixer.in","w") as os:
    os.write(generate_input.h2mixer_input(params))
with open("mfdn.dat","w") as os:
    os.write(generate_input.mfdn_input(mfdn_params))
