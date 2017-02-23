""" jpv2rel.py

    Scripting for jpv2rel conversion runs.

    M. A. Caprio
    Department of Physics
    University of Notre Dame

    2/23/17 (mac): Created.
"""

# generic packages for use below
import math
import os
import sys

# load mcscript
import mcscript
##mcscript.init()

##################################################################
# build task list
##################################################################

# executables

# executable files
projects_root = os.path.join(os.environ["HOME"],"projects")
jpv2rel_executable = os.path.join(projects_root,"shell","libraries","relative","jpv2rel")

# configuration
hw_values = mcscript.utils.value_range(10,40,2.5)

# generate task list

# This should be a list of dictionaries, with keys that your task
# handler function will understand.  Then, mcscript will automatically
# add some extra keys to these dictionaries, e.g., a task identifier
# name.

tasks = [
    {
        "hw" : hw,
        "Jmax" : 4,
        "Nmax" : 20,
        "jpv_filename_template" : "Vrel_JISP16_bare_Jmax{Jmax}.hw{hw}",
        "rel_filename_template" : "jisp16_Nmax{Nmax}_hw{hw:2.1f}_rel.dat"
    }
    for hw in hw_values
]

def convert_jpv2rel(task):
    """ Invoke jpv2rel.
    """
    input_lines = [
        "{Nmax} {Jmax}".format(**task),
        task["jpv_filename_template"].format(**task),
        task["rel_filename_template"].format(**task)
    ]
    mcscript.call(jpv2rel_executable,input_lines=input_lines)
    

for task in tasks:
    convert_jpv2rel(task)
