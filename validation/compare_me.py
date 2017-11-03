""" compare_me.py

    Compare matrix elements in two files, where row order may differ.

    Entries and row order of first file are used.  Entries in second file are "looked up".

    M. A. Caprio
    Department of Physics
    University of Notre Dame

    10/18/17 (mac): Created.
"""

# generic packages for use below
import math
import os
import sys

import numpy as np

# load mcscript
import mcscript
##mcscript.init()

################################################################
# library of format strings
################################################################

# these strings are just for the quantum numbers, not including the RME itself

# LSJT format from basis::lsjt_operator.h:
#
#   T0  N1' l1' N2' l2' L' S' J' T' g'  N1 l1 N2 l2 L S J T g  JT-RME

conversion_list_lsjt = (1+9+9)*[int]
format_string_lsjt = "{:1d} " + (2*9)*"{:2d}"

# JJJT format from basis::jjjt_operator.h:
#
# T0  N1' l1' j1' N2' l2' j2' J' T' g'  N1 l1 j1 N2 l2 j2 J T g  JT-RME

conversion_list_jjjt = [int] + 2*[int,int,float,int,int,float,int,int,int]
format_string_jjjt = "{:1d} " + 2*"{:2d} {:2d} {:4.1f} {:2d} {:2d} {:4.1f} {:2d} {:1d} {:1d} "

# h2 format v0

conversion_list_h2v0 = [int,int,int,int,int,int]
format_string_h2v0 = 6*"{:3d} "


##################################################################
# input routine
##################################################################

def split_and_prune_lines(lines):
    """Split input lines into tokenized lines, suppressing comment
    (beginning with hash token) or empty lines.

    Tokenization is by whitespace, i.e., with split.

    The tokenized lines are returned as tuples rather than lists to
    support downstream processing.  For instance, structured array
    creation with np.array requires the entries to be tuples.

    Arguments:
       (iterable of str): input lines

    Returns:
       (iterator of tuple of str): split and filtered lines

    COPY from mfdnres.tools

    """
    
    def is_active_line(tokens):
        """ Identify nonempty, noncomment line.

        Helper function for file parsing.

        Arguments:
            (list of str): tokenized line

        Returns:
            (bool)
        """
        return bool(tokens) and (tokens[0]!="#")

    tokenized_lines = map(lambda s : tuple(s.split()),lines)
    return filter(is_active_line,tokenized_lines)


def read_me_file(in_stream,conversion_specifier,num_header_lines=0):
    """ Read and store matrix element.

    Arguments:
        in_stream (iterable of lines):
        conversion_specifier (iterable of callables): conversion functions to apply to tokens
        num_header_lines (int,optional): number of header lines to skip

    Returns:
        matrix_element_labels (list of tuble): labels in order of input
        matrix_element_dict (dict): matrix elements indexed by label tuple
             labels (tuple of float): labels all converted to float
             matrix_element (float): matrix element
    """

    # get input tokens
    tokenized_lines = split_and_prune_lines(in_stream)
    tokenized_lines = list(tokenized_lines)[num_header_lines:]  # decapitate

    matrix_element_dict = {}
    matrix_element_labels = []
    for tokenized_line in tokenized_lines:

        # extract data from line
        labels = tuple([conversion(token) for (conversion,token) in zip(conversion_specifier,tokenized_line[:-1])])
        matrix_element = float(tokenized_line[-1])

        # save data
        matrix_element_labels.append(labels)
        matrix_element_dict[labels] = matrix_element

    return (matrix_element_labels,matrix_element_dict)

def print_matrix_elements(out_stream,matrix_element_dict,label_format_specifier,matrix_element_format_specifier):
    """ Print stored matrix elements, in lexicographic order by labels.

    Arguments:
        out_stream (stream): output stream
        matrix_element_dict (dict): matrix elements indexed by label tuple
        label_format_specifier (str): format string for label tuple
        matrix_element_format_specifier (str): format string for matrix element
    """

    sorted_labels = sorted(matrix_element_dict.keys())
    for labels in sorted_labels:
        matrix_element = matrix_element_dict[labels]
        out_line = label_format_specifier.format(*labels) + format(matrix_element,matrix_element_format_specifier)
        print(out_line,file=out_stream)

def compare_matrix_elements(
        out_stream,
        matrix_element_dict1,matrix_element_dict2,
        label_format_specifier,matrix_element_format_specifier,
        threshold,
        matrix_element_labels=None
):
    """ Print stored matrix elements from two files, in lexicographic order by labels, taking labels from first file.

    Arguments:
        out_stream (stream): output stream
        matrix_element_dict (dict): matrix elements indexed by label tuple
        label_format_specifier (str): format string for label tuple
        matrix_element_format_specifier (str): format string for matrix element
        threshold (float): difference threshold for flag
        matrix_element_labels (list of tuple, optional): ordered listing of labels for matrix elements for retrieval
    """

    if (matrix_element_labels is None):
        matrix_element_labels = sorted(matrix_element_dict1.keys())

    for labels in matrix_element_labels:
        matrix_element1 = matrix_element_dict1[labels]
        matrix_element2 = matrix_element_dict2.setdefault(labels,0.)
        difference = matrix_element2-matrix_element1
        try:
            ratio = matrix_element2/matrix_element1
        except:
            ratio = np.nan
        flag = " ***" if abs(difference)>threshold else ""
        out_line = "{:s}  {:s}  {:s}     {:s}     {:s}{:s}".format(
            label_format_specifier.format(*labels),
            format(matrix_element1,matrix_element_format_specifier),
            format(matrix_element2,matrix_element_format_specifier),
            format(difference,matrix_element_format_specifier),
            format(ratio,"+10.3e"),  ## might choose different format
            flag
        )
        print(out_line,file=out_stream)

def simple_test():
    with open("171018-moshinsky-test-files/jisp16_tb_Nmax04.dat") as in_file:
        (matrix_element_labels,matrix_element_dict) = read_me_file(in_file,conversion_list_jjjt,num_header_lines=11)

    with open("test1.dat","w") as out_stream:
        print_matrix_elements(out_stream,matrix_element_dict,format_string_jjjt,"+10.6e")

    with open("test2.dat","w") as out_stream:
        compare_matrix_elements(out_stream,matrix_element_dict,matrix_element_dict,format_string_jjjt,"+16.8e",1e-4,matrix_element_labels)

def compare_jisp16_jjjt():
    """ JISP16: Compare moshinsky with Anna."""
    with open("jisp16_Nmax04_hw20.0_jjjt.dat") as in_stream:
        (matrix_element_labels,matrix_element_dict1) = read_me_file(in_stream,conversion_list_jjjt,num_header_lines=0)
    with open("171018-moshinsky-test-files/jisp16_tb_Nmax04.dat") as in_stream:
        (_,matrix_element_dict2) = read_me_file(in_stream,conversion_list_jjjt,num_header_lines=11)

    with open("jisp16_Nmax04_hw20.0_jjjt-COMPARISON.dat","w") as out_stream:
        compare_matrix_elements(out_stream,matrix_element_dict1,matrix_element_dict2,format_string_jjjt,"+16.8e",1e-4,matrix_element_labels)

def compare_coulomb_jjjt():
    """ Coulomb: Compare moshinsky with Anna."""
    with open("coulomb_Nmax04_jjjt.dat") as in_stream:
        (matrix_element_labels,matrix_element_dict1) = read_me_file(in_stream,conversion_list_jjjt,num_header_lines=0)
    with open("171018-moshinsky-test-files/coulomb_tb_Nmax04.dat") as in_stream:
        (_,matrix_element_dict2) = read_me_file(in_stream,conversion_list_jjjt,num_header_lines=11)

    with open("coulomb_Nmax04_jjjt-COMPARISON.dat","w") as out_stream:
        compare_matrix_elements(out_stream,matrix_element_dict1,matrix_element_dict2,format_string_jjjt,"+16.8e",1e-4,matrix_element_labels)

def compare_nnloopt_jjjt():
    """ NNLOopt: Compare moshinsky with Anna."""
    with open("nnloopt_Nmax04_hw40.0_caveat-nmax30_jjjt.dat") as in_stream:
        (matrix_element_labels,matrix_element_dict1) = read_me_file(in_stream,conversion_list_jjjt,num_header_lines=0)
    with open("171020-moshinsky-test-files/nnlo-opt_tb_Nmax04.dat") as in_stream:
        (_,matrix_element_dict2) = read_me_file(in_stream,conversion_list_jjjt,num_header_lines=11)

    with open("nnloopt_Nmax04_hw40.0_caveat-nmax30_jjjt-COMPARISON.dat","w") as out_stream:
        compare_matrix_elements(out_stream,matrix_element_dict1,matrix_element_dict2,format_string_jjjt,"+16.8e",1e-4,matrix_element_labels)

def compare_quadrupole_lsjt():
    """ Quadrupole: Compare moshinsky with Anna."""
    with open("../work/moshinsky/quadrupole/quadrupole_Nmax04_total_lsjt.dat") as in_stream:
        (matrix_element_labels,matrix_element_dict1) = read_me_file(in_stream,conversion_list_lsjt,num_header_lines=0)
    with open("171101-moshinsky-test-files/quadrupole_Nmax04_total_lsjt_spncci.dat") as in_stream:
        (_,matrix_element_dict2) = read_me_file(in_stream,conversion_list_jjjt,num_header_lines=11)

    with open("quadrupole_Nmax04_total_lsjt-COMPARISON.dat","w") as out_stream:
        compare_matrix_elements(out_stream,matrix_element_dict1,matrix_element_dict2,format_string_jjjt,"+16.8e",1e-4,matrix_element_labels)

def compare_quadrupole_jjjt():
    """ Quadrupole: Compare moshinsky with Anna."""
    with open("../work/moshinsky/quadrupole/quadrupole_Nmax04_total_jjjt.dat") as in_stream:
        (matrix_element_labels,matrix_element_dict1) = read_me_file(in_stream,conversion_list_jjjt,num_header_lines=0)
    with open("171101-moshinsky-test-files/quadrupole_Nmax04_total_jjjt_spncci.dat") as in_stream:
        (_,matrix_element_dict2) = read_me_file(in_stream,conversion_list_jjjt,num_header_lines=11)

    with open("quadrupole_Nmax04_total_jjjt-COMPARISON.dat","w") as out_stream:
        compare_matrix_elements(out_stream,matrix_element_dict1,matrix_element_dict2,format_string_jjjt,"+16.8e",1e-4,matrix_element_labels)

def compare_jisp16_h2v0():
    """ JISP16: Compare moshinsky with h2."""
    with open("jisp16_Nmax06_hw20.0_h2v0.dat") as in_stream:
        (matrix_element_labels,matrix_element_dict1) = read_me_file(in_stream,conversion_list_h2v0,num_header_lines=5)
    with open("run0164-JISP16-coul-tb-06-hw20-dat/JISP16-tb-6-20.dat") as in_stream:
        (_,matrix_element_dict2) = read_me_file(in_stream,conversion_list_h2v0,num_header_lines=5)

    with open("jisp16_Nmax06_hw20.0_h2v0-COMPARISON.dat","w") as out_stream:
        compare_matrix_elements(out_stream,matrix_element_dict1,matrix_element_dict2,format_string_h2v0,"+16.8e",1e-4,matrix_element_labels)

def compare_coulomb_h2v0():
    """ Coulomb: Compare moshinsky with h2."""
    with open("coulomb_Nmax06_h2v0.dat") as in_stream:
        (matrix_element_labels,matrix_element_dict1) = read_me_file(in_stream,conversion_list_h2v0,num_header_lines=5)
    with open("run0164-JISP16-coul-tb-06-hw20-dat/VC-tb-6-20.dat") as in_stream:
        (_,matrix_element_dict2) = read_me_file(in_stream,conversion_list_h2v0,num_header_lines=5)

    with open("coulomb_Nmax06_h2v0-COMPARISON.dat","w") as out_stream:
        compare_matrix_elements(out_stream,matrix_element_dict1,matrix_element_dict2,format_string_h2v0,"+16.8e",1e-4,matrix_element_labels)

def compare_nnloopt_h2v0():
    """ NNLOopt: Compare moshinsky with h2."""
    with open("nnloopt_Nmax06_hw40.0_caveat-nmax30_h2v0.dat") as in_stream:
        (matrix_element_labels,matrix_element_dict1) = read_me_file(in_stream,conversion_list_h2v0,num_header_lines=5)
    with open("run0306-N2LOopt500-tb-06-hw40-dat/N2LOopt500-tb-6-40.dat") as in_stream:
        (_,matrix_element_dict2) = read_me_file(in_stream,conversion_list_h2v0,num_header_lines=5)

    with open("nnloopt_Nmax06_hw40.0_h2v0-COMPARISON.dat","w") as out_stream:
        compare_matrix_elements(out_stream,matrix_element_dict1,matrix_element_dict2,format_string_h2v0,"+16.8e",1e-4,matrix_element_labels)

def main():

    # basic input
    simple_test()

    # jjjt comparisons
    ## compare_jisp16_jjjt()
    ## compare_coulomb_jjjt()
    ## compare_nnloopt_jjjt()
    compare_quadrupole_lsjt()
    compare_quadrupole_jjjt()

    # h2 comparisons
    ## compare_jisp16_h2v0()
    ## compare_coulomb_h2v0()
    ## compare_nnloopt_h2v0()

if (__name__=="__main__"):
    main()
