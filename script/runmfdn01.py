""" runmfdn01.py -- h2mixer+mfdn mscript test run

    Mark A. Caprio
    University of Notre Dame

    - 12/14/16 (mac): Created.
      
"""

import mfdn_h2

# load mcscript
import mcscript
mcscript.init()

##################################################################
# build task list
##################################################################

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
    "h2_format" : 0  # h2 formats: 0, 15099
}

# generate task list

tasks = [
    params
]

##################################################################
# implementation functions for doing a "hello world" task
#
# For a more complicated application, you would separate these out
# into their own module.
##################################################################

def task_descriptor_hello (current_task):
    """ Return task descriptor for hello task.
    """

    return "{world_name}".format(**current_task)

def task_handler_hello (current_task):
    """ Do a hello world task  given current task parameters.

    Expected dictionary keys:

    "world_name" : name of world to greet

    """

    # write greeting message to file

    mcscript.utils.write_input(
        "hello.txt",
        input_lines=[
            "Dear {world_name},".format(**current_task),
            "   Hello!",
            "Your script",
            mcscript.run.name
            ]
        )


##################################################################
# task list entry annotation functions
##################################################################

def task_pool (current_task):
    """ Create task pool identifier.
    """
    
    return "greet"

##################################################################
# master loop
##################################################################

mcscript.task.init(
    tasks,
    task_descriptor=task_descriptor_hello,
    task_pool=task_pool,
    phase_handler_list=[mfdn_h2.task_handler_ho]
    )

################################################################
# termination
################################################################

mcscript.termination()
