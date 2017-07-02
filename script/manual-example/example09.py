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
    "Nmax_orb" : 13,  # for target space
    "Nmax_orb_int" : 13,  # for overlaps from interaction tbmes
    "orbitals_filename" : "orbitals.dat",  # for target space
    "orbitals_int_filename" : "orbitals-int.dat",  # for overlaps from interaction tbmes
    "radial_me_filename_pattern" : "radial-me-{}{}.dat",  # "{}{}" will be replaced by {"r1","r2","k1","k2"}
    "radial_olap_filename" : "radial-olap.dat",
    # stored input TBMEs
    "VNN_filename" : os.path.join(generate_input.data_dir_h2,"run0164-JISP16-ob-9/JISP16-ob-9-20.bin"),
    "VC_filename" : os.path.join(generate_input.data_dir_h2,"run0164-JISP16-ob-9/VC-ob-9-20.bin"),
    # scaling/transformation/basis parameters
    "hw" : 20,
    "target_truncation" : ("ob",13),
    "hw_int" : 20,
    "hw_c" : 20,
    "xform_enabled" : False,
    "xform_truncation" : ("ob",9),
    # Hamiltonian parameters
    "A" : 2,   # atomic mass number of nucleus (for kinematic operators)
    "a" : 20,  # coefficient on Ncm in Lawson term
    "hw_cm" : None,  # hw for Ncm in Lawson term; if None, defaults to "hw"
    "use_coulomb" : False,
    # output file format
    "h2_format" : 0  # h2 formats: 0, 15099
}

with open("run.csh","w") as os:
    os.write(generate_input.run_script(params))
with open("h2mixer.in","w") as os:
    os.write(generate_input.h2mixer_input(params))
