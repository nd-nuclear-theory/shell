"""mfdnmem.py -- memory estimates for MFDn parallel runs

    Example usage:
        python3 mfdnmem.py < example.in > example.out

    Input format:

        line 1: n_v n_lan

            n_v (int) : number of stored eigenvectors
            n_lan (int) : number of Lanczos iterations (and thus Lanczos vectors)

        line 2: memory_node_GB core_node thread_nums

            memory_node_GB (float) : memory per node (in GB)
            core_node (int): cores per node
            thread_nums (comma-separated ints): OMP depths to consider

        The following lines follow the format of Pieter Maris'z "sparsity" tabulations...

        line n: Z N Nmax Mj D n_nz

            Z (int), N (int), Nmax (int), Mj (float): labels for output
            D (int or float) : dimension of matrix
            n_nz (int or float) : number of nonzero matrix elements in matrix

        Any extra junk on a line is ignored.
        Lines after header beginning with "#" are ignored.

    ----------------------------------------------------------------

    Example input:
        15 500
        64 24 1,6 # Edison configuration
          4      5        8   0.5     63003395        45472165505
          4      5        9   0.5    196163784       181176686001
          4      5       10   0.5    574827349       668289364616
          4      5       11   0.5   1597056996      2304776012970
          4      5       12   0.5   4232122420      7490682147171

    Example output:

        Parameters
            n_v 15 n_lan 500
            memory/node (GB) 64.000
            core/node 24
            memory/core (GB) 2.667
            size range (GB) 0.667..2.667

        Z 4 N 5 Nmax  8 Mj 0.5

          depth 1
             n_d    n_p   n_n  m1 (GB)  m2 (GB)   m (GB)    m (%)
            ---- ------ -----  -------  -------  -------  -------
              25    325    14    1.404    0.803    1.404     52.6
              35    630    27    0.724    0.470    0.724     27.2

          depth 6
             n_d    n_p   n_n  m1 (GB)  m2 (GB)   m (GB)    m (%)
            ---- ------ -----  -------  -------  -------  -------
               9     45    12    1.689    0.758    1.689     63.4
              11     66    17    1.152    0.534    1.152     43.2
              13     91    23    0.835    0.401    0.835     31.3

        Z 4 N 5 Nmax  9 Mj 0.5

          depth 1
             n_d    n_p   n_n  m1 (GB)  m2 (GB)   m (GB)    m (%)
            ---- ------ -----  -------  -------  -------  -------
              45   1035    44    1.657    1.139    1.657     62.1
              55   1540    65    1.114    0.837    1.114     41.8
              65   2145    90    0.800    0.652    0.800     30.0

          depth 6
             n_d    n_p   n_n  m1 (GB)  m2 (GB)   m (GB)    m (%)
            ---- ------ -----  -------  -------  -------  -------
              15    120    30    2.382    1.181    2.382     89.3
              17    153    39    1.868    0.950    1.868     70.1
              19    190    48    1.505    0.784    1.505     56.4
              21    231    58    1.238    0.661    1.238     46.4
              23    276    69    1.036    0.566    1.036     38.8
              25    325    82    0.880    0.492    0.880     33.0
              27    378    95    0.756    0.433    0.756     28.4



    ----------------------------------------------------------------

    Language: Python 3 (compatible with Python 2.7)
    Mark A. Caprio
    University of Notre Dame

    6/30/15 (mac): Initiated.
    11/23/15 (mac): Adapt to take multiple Nmax inputs.
    10/20/17 (pjf): Adapt to node/rank oriented scripting.
    Last modified 11/23/15.

"""

from __future__ import print_function, division

## import configparser
import sys
import math

################################################################
# configuration
################################################################

# constants
GB = 2.**30  # size of GB in bytes

# meshes for run configuration
DIAGONAL_MESH = (
    list(range(1, 15, 2))
    + list(range(15, 155, 10))
    + list(range(155, 1000, 20))
)
DIAGONAL_MESH_OMP = (
    list(range(1, 15, 2))
    + list(range(15, 155, 2))   # make finer for hybrid run
    + list(range(155, 1000, 8))
)


################################################################
# parameter input
################################################################

class MFDnDimensionParameters(object):
    """ Object to store dimension parameters applicable to an MFDn run.

    Attributes:
        n_nz (int or float) : number of nonzero matrix elements in matrix
        D (int or float) : dimension of matrix
        n_lan (int) : number of Lanczos iterations (and thus Lanczos vectors)
        n_v (int) : number of stored eigenvectors
        n_d (int) : number of diagonal processes
        n_p (int) : number of processes (derived from n_d)

    Methods:
        set_diagonal(n_d) : sets both n_d and the derived value n_p
    """

    def set_diagonal(self, n_d):
        """ Sets both n_d and the derived value n_p.

        Args:
            n_d (int) : number of diagonal processes
        """

        self.n_d = n_d
        self.n_p = n_d * (n_d + 1) // 2

################################################################
# memory calculation
################################################################


def memory_per_process(params):
    """ Estimates memory needed for diagonal processes.

    Args:
       params (MFDnDimensionParameters): run dimension parameters

    Returns:
        (tuple of float) : sizes (size_1,size_2) in bytes of estimated
        storage per diagonal process, in the two calculation phases
    """

    # calculate floats stored per MPI process for different purposes
    # number of floats on generic process for matrix storage
    nu_mat = (params.n_nz) / (params.n_p)
    # number of floats on generic process for Lanczos vector storage
    nu_lan = (params.n_lan) * (params.D) / (params.n_p)
    # number of floats on diagonal process for observable calculation
    nu_v = 2 * (params.n_v) * (params.D) / (params.n_d)

    # calculate total memory requirement for diagonal process
    #
    # Note: Naively, a diagonal process would only have nu_mat/2
    # floats to store for the matrix, due to symmetry on the diagonal.
    # However, sparsity is lower on the diagonal, which more than
    # compensates.  We take nu_mat as the estimate, but the number can
    # be even higher.
    #
    # phase 1 (eigensolver)
    #   Diagonal process must store:
    #      nu_mat float32 (for matrix elements)
    #      nu_mat int32 (for matrix indexing)
    #      nu_lan float32 (for Lanczos vectors)
    #
    # phase 2 (observables)
    #   Diagonal process must store:
    #      nu_mat int32 (for matrix indexing)
    #      nu_obs float32 (for eigenvectors)
    #   No matrix elements are stored for TBOs, since those are constructed on the fly.

    float_size = 4
    int_size = 4
    size_1 = nu_mat * float_size + nu_mat * int_size + nu_lan * float_size
    size_2 = nu_mat * int_size + nu_v * float_size

    return (size_1, size_2)

################################################################
# memory calculation
################################################################


def tabulate_usage(params, size_range, threads, core_node, mesh):
    """ Tabulates possible run sizes.

    Args:
       params (MFDnDimensionParameters): run dimension parameters
       size_range (tuple of float) : acceptible range (size_min,size_max) of
           memory usage per core (in bytes)
       threads (int) : OMP threads/rank
       core_node (int) : cores / node
       mesh (list of int) : list of mesh sizes to consider
    """

    (size_min, size_max) = size_range
    template_header = "    {n_d:>4s} {n_p:>6s} {n_n:>5s}  {mem_1:>7s}  {mem_2:>7s}  {mem:>7s}  {mem_perc:>7s}"
    template_line = "    {n_d:4d} {n_p:6d} {n_n:5d}  {mem_1:7.3f}  {mem_2:7.3f}  {mem:7.3f}  {mem_perc:7.1f}"
    print(template_header.format(n_d="n_d", n_p="n_p", n_n="n_n",
                                 mem_1="m1 (GB)", mem_2="m2 (GB)", mem="m (GB)",
                                 mem_perc="m (%)"))
    print(template_header.format(n_d="-" * 4, n_p="-" * 6, n_n="-" * 5,
                                 mem_1="-" * 7, mem_2="-" * 7, mem="-" * 7,
                                 mem_perc="-" * 7))

    # iterate over possible depths
    for n_d in mesh:
        # estimate memory requirement
        params.set_diagonal(n_d)
        (size_1, size_2) = memory_per_process(params)
        (size_core_1, size_core_2) = (size_1 / threads, size_2 / threads)
        size_core = max(size_core_1, size_core_2)  # memory per core
        num_nodes = math.ceil((params.n_p * threads) / core_node)

        # process according to size
        if size_core > size_max:
            # ignore if still have too few cores
            continue
        elif size_core > size_min:
            # memory usage is in acceptable range
            print(template_line.format(
                n_d=params.n_d, n_p=params.n_p, n_n=num_nodes,
                mem_1=size_core_1 / GB, mem_2=size_core_2 / GB,
                mem=size_core / GB, mem_perc=100 * size_core / size_max)
            )
        else:
            # we have too many cores now
            break


################################################################
# main
################################################################

if (__name__ == "__main__"):

    # test case: manually set up dimensions
    #
    # number of processors will be set in tabulation below
    ## params = MFDnDimensionParameters()
    ## params.n_nz = 2.07e11
    ## params.D = 1.87e8
    ## params.n_lan = 1000
    ## params.n_v = 30

    # header input
    # set up structure for input
    params = MFDnDimensionParameters()
    # read line 1: n_v n_lan
    line = sys.stdin.readline()
    tokens = line.split()
    params.n_v = int(tokens[0])
    params.n_lan = int(tokens[1])
    # read line 2: memory_core_GB
    line = sys.stdin.readline()
    tokens = line.split()
    memory_node_GB = float(tokens[0])
    core_node = int(tokens[1])
    thread_nums = [int(val) for val in tokens[2].split(',')]
    memory_node = memory_node_GB * GB
    memory_core_GB = memory_node_GB / core_node
    memory_core = memory_node / core_node

    # set memory range
    (size_min, size_max) = (memory_core / 4, memory_core)

    # head output
    print("Parameters")
    params_line = (
        "    n_v {params.n_v:d} n_lan {params.n_lan:d}\n"
        "    memory/node (GB) {memory_node_GB:.3f}\n"
        "    core/node {core_node:d}\n"
        "    memory/core (GB) {memory_core_GB:.3f}"
    )
    print(params_line.format(params=params, memory_node_GB=memory_node_GB,
                             core_node=core_node, memory_core_GB=memory_core_GB))
    print("    size range (GB) {:.3f}..{:.3f}".format(
        size_min / GB, size_max / GB))
    print()

    # iterate over Nmax cases
    for line in sys.stdin:

        # read line: Z N Nmax Mj D n_nz
        if line.strip() == "":
            continue
        tokens = line.split()
        if tokens[0][0] == "#":
            continue
        params.Z = int(tokens[0])
        params.N = int(tokens[1])
        params.Nmax = int(tokens[2])
        params.Mj = float(tokens[3])
        params.D = float(tokens[4])
        params.n_nz = float(tokens[5])

        # head Nmax output
        params_line = (
            "Z {params.Z:d} N {params.N:d} Nmax {params.Nmax:2d} Mj {params.Mj:.1f} "
        )
        print(params_line.format(params=params))
        print()

        # make tables
        for threads in thread_nums:
            if threads == 1:
                mesh = DIAGONAL_MESH
            else:
                mesh = DIAGONAL_MESH_OMP
            print("  depth {}".format(threads))
            tabulate_usage(params, (size_min, size_max),
                           threads, core_node, mesh)
            print()
