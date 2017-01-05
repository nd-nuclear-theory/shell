"""mfdnmem.py -- memory estimates for MFDn parallel runs

    Example usage:
        python3 mfdnmem.py < example.in > example.out

    Input format:

        line 1: n_v n_lan

            n_lan (int) : number of Lanczos iterations (and thus Lanczos vectors)
            n_v (int) : number of stored eigenvectors

        line 2: memory_core_Gb

            memory_core_Gb (float) : memory per core (in Gb)

        The following lines follow the format of Pieter Maris'z "sparsity" tabulations...

        line n: Z N Nmax Mj D n_nz

            Z (int), N (int), Nmax (int), Mj (float): labels for output
            n_nz (int or float) : number of nonzero matrix elements in matrix
            D (int or float) : dimension of matrix

        Any extra junk on a line is ignored.

    ----------------------------------------------------------------

    Example input:
        10 1000
        1.333
          6	 6	 4	 0	1118926		279405126	5239055438	44134939578
          6	 6	 5	 0	6488004		2320030394	56590581542	620555371150
          6	 6	 6	 0	32598920	15942738578	486945885186	6677973874918
          6	 6	 7	 0	145989820	94374647128	3511904931196	58544810332902
          6	 6	 8	 0	594496743	494946348099	21979112668303	4.3e14
          6	 6	 9	 0	2233681236	2345846247832	1.22e14
          6	 6	10	 0	7830355795	10199749477279
          6	 6	11	 0	25835492536	4.2e13

    Example output:

        Parameters
            n_v 10 n_lan 1000
            memory/core (Gb) 1.333
            size range (Gb) 0.333..1.333
        
        Z 6 N 6 Nmax  4 Mj 0.0 
        
          depth 1
             n_d    n_p  m1 (Gb)  m2 (Gb)   m (Gb)
            ---- ------  -------  -------  -------
               3      6    1.042    0.201    1.042
               5     15    0.417    0.086    0.417
        
          depth 6
             n_d    n_p  m1 (Gb)  m2 (Gb)   m (Gb)
            ---- ------  -------  -------  -------
               1      1    1.042    0.187    1.042
        
        Z 6 N 6 Nmax  5 Mj 0.0 
        
          depth 1
             n_d    n_p  m1 (Gb)  m2 (Gb)   m (Gb)
            ---- ------  -------  -------  -------
               9     45    0.921    0.246    0.921
              11     66    0.628    0.175    0.628
              13     91    0.456    0.132    0.456
              15    120    0.345    0.104    0.345
        
          depth 6
             n_d    n_p  m1 (Gb)  m2 (Gb)   m (Gb)
            ---- ------  -------  -------  -------
               3      6    1.152    0.267    1.152
               5     15    0.461    0.112    0.461
        
        ...
        
             995 495510    0.138    0.375    0.375



    ----------------------------------------------------------------

    Language: Python 3 (compatible with Python 2.7)
    Mark A. Caprio
    University of Notre Dame

    6/30/15 (mac): Initiated.
    11/23/15 (mac): Adapt to take multiple Nmax inputs.
    Last modified 11/23/15.

"""

from __future__ import print_function, division

## import configparser
import sys

################################################################
# configuration
################################################################

# constants
Gb = 2.**30  # size of Gb in bytes

# meshes for run configuration
DEPTH_MESH = [1,6]  # list of OMP depths to consider
DIAGONAL_MESH = (
    list(range(1,15,2)) 
    + list(range(15,155,10))
    + list(range(155,1000,20))
    )
DIAGONAL_MESH_OMP = (
    list(range(1,15,2)) 
    + list(range(15,155,2))   # make finer for hybrid run
    + list(range(155,1000,20))
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

    def set_diagonal(self,n_d):
        """ Sets both n_d and the derived value n_p.
        
        Args:
            n_d (int) : number of diagonal processes
        """
    
        self.n_d = n_d
        self.n_p = n_d*(n_d+1)//2
        
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
    nu_mat = (params.n_nz)/(params.n_p)    # number of floats on generic process for matrix storage 
    nu_lan = (params.n_lan)*(params.D)/(params.n_p)    # number of floats on generic process for Lanczos vector storage
    nu_v = 2*(params.n_v)*(params.D)/(params.n_d) # number of floats on diagonal process for observable calculation

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
    size_1 =  nu_mat*float_size + nu_mat*int_size + nu_lan*float_size
    size_2 =  nu_mat*int_size + nu_v*float_size
    
    return (size_1,size_2)

################################################################
# memory calculation
################################################################

def tabulate_usage(params,size_range,depth,mesh):
    """ Tabulates possible run sizes.

    Args:
       params (MFDnDimensionParameters): run dimension parameters
       size_range (tuple of float) : acceptible range (size_min,size_max) of
           memory usage per core (in bytes)
       depth (int) : OMP depth
       mesh (list of int) : list of mesh sizes to consider
    """
    
    (size_min,size_max) = size_range
    template_header = "    {n_d:>4s} {n_p:>6s}  {mem_1:>7s}  {mem_2:>7s}  {mem:>7s}"
    template_line = "    {n_d:4d} {n_p:6d}  {mem_1:7.3f}  {mem_2:7.3f}  {mem:7.3f}"
    print(template_header.format(n_d="n_d",n_p="n_p",mem_1="m1 (Gb)",mem_2="m2 (Gb)",mem="m (Gb)"))
    print(template_header.format(n_d="-"*4,n_p="-"*6,mem_1="-"*7,mem_2="-"*7,mem="-"*7))

    # iterate over possible depths
    for n_d in mesh:
        # estimate memory requirement
        params.set_diagonal(n_d)
        (size_1,size_2) = memory_per_process(params)
        (size_core_1,size_core_2) = (size_1/depth,size_2/depth)
        size_core = max(size_core_1,size_core_2)  # memory per core
        
        # process according to size
        if (size_core > size_max):
            # ignore if still have too few cores
            continue
        elif (size_core > size_min):
            # memory usage is in acceptable range
            print(template_line.format(
                n_d=params.n_d,n_p=params.n_p,
                mem_1=size_core_1/Gb,mem_2=size_core_2/Gb,mem=size_core/Gb)
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
    # read line 2: memory_core_Gb
    line = sys.stdin.readline()
    tokens = line.split()
    memory_core_Gb = float(tokens[0])
    memory_core = memory_core_Gb * Gb

    # set memory range
    (size_min,size_max) = (memory_core/4,memory_core)

    # head output
    print("Parameters")
    params_line = (
        "    n_v {params.n_v:d} n_lan {params.n_lan:d}\n"
        "    memory/core (Gb) {memory_core_Gb:.3f}"
    )
    print(params_line.format(params=params,memory_core_Gb=memory_core_Gb))
    print("    size range (Gb) {:.3f}..{:.3f}".format(size_min/Gb,size_max/Gb))
    print()

    # iterate over Nmax cases
    for line in sys.stdin:

        # read line: Z N Nmax Mj D n_nz
        if (line.strip() == ""):
            continue
        tokens = line.split()
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
        for depth in DEPTH_MESH:
            if (depth == 1):
                mesh = DIAGONAL_MESH
            else:
                mesh = DIAGONAL_MESH_OMP
            print("  depth {}".format(depth))
            tabulate_usage(params,(size_min,size_max),depth,mesh)
            print()
