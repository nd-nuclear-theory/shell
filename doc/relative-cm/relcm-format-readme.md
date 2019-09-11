# Notre Dame LSJT relative-cm interaction format #

URL: https://github.com/nd-nuclear-theory/shell/doc/relative/lsjt-relative-format-readme.txt

History:
  + 04/07/19 (pjf): Extract descriptions from basis commit e832004.

----------------------------------------------------------------

The LSJT relative-cm matrix element file format is defined in header file comments:

  - The Notre Dame LSJT relative-cm operator file format is defined in

    libraries/basis/lsjt_operator.h

  - To understand the ordering of sectors and matrix elements in the
  operator file, you need to know the ordering of basis subpaces and
  states, which is defined in

    libraries/basis/lsjt_basis.h

The relavant comments from both these files are copied below for your
reference.  Please refer to the original files for further information
on coding (internal data structures, I/O routines, etc.).

Then you can find many example files under the examples subdirectory.

----------------------------------------------------------------

Basis description (lsjt_scheme.h):

  Note: In the following, "~" indicates equality modulo 2.

    ////////////////////////////////////////////////////////////////
    // relative-cm states in LSJT scheme
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    //
    // Labeling
    //
    // subspace labels: (L,S,J,T,g)
    //
    //   L (int): orbital angular momentum
    //   S (int): total spin
    //   J (int): total angular momentum
    //   T (int): isospin
    //   g (int): grade (=0,1) for the parity P
    //
    // state labels within subspace: (Nr,lr,Nc,lc)
    //
    //   Nr (int): oscillator quanta of relative motion
    //   lr (int): orbital angular momentum of relative motion
    //   Nc (int): oscillator quanta of c.m. motion
    //   lc (int): orbital angular momentum of c.m. motion
    //
    ////////////////////////////////////////////////////////////////
    //
    // Subspaces
    //
    // Within the full space, subspaces are ordered by:
    //    -- increasing L (L=0,1,...,Nmax)
    //    -- increasing S (S=0,1)
    //    -- increasing J
    //    -- increasing T (T=0,1)
    //    -- increasing g (g=0,1)
    // subject to:
    //    -- triangularity of (L,S,J)
    //
    // Subspaces are pruned to those of nonzero dimension.
    //
    // Note that ordering of subspaces is lexicographic by (L,S,J,T,g).
    //
    // Truncation of the space is by the two-body Nmax.
    //
    ////////////////////////////////////////////////////////////////
    //
    // States
    //
    // Within a subspace, the states are ordered by:
    //   -- increasing N  (N=Nr+Nc)
    //   -- lexicographically increasing (Nr,lr)
    //   -- lexicographically increasing (Nc,lc)
    // and subject to:
    //   -- triangularity constraint on (lr,lc,L)
    //   -- parity constraint N~g
    //   -- antisymmetry constraint lr+S+T~1 (or, equivalentsly,
    //      Nr+S+T~1)
    //
    // This basis is for *identical* particle states, as enforced by the
    // antisymmetry constraint on Nr.
    //
    ////////////////////////////////////////////////////////////////

Operator file format description (lsjt_operator.h):

    ////////////////////////////////////////////////////////////////
    //
    // general relative-cm two-body operator file format
    //
    ////////////////////////////////////////////////////////////////
    //
    // Header
    //
    //   Header format:
    //
    //     # RELATIVE-CM LSJT
    //     # ...
    //     version                                 # version
    //     J0 g0 T0_min T0_max symmetry_phase_mode # operator tensor properties
    //     Nmax                                    # relative-cm basis truncation
    //
    //   The header may start with one or more contiguous comment lines,
    //   which are designated by a hash character in the first column.
    //
    //   Then, the header contains the following fields:
    //
    //     version : file format version (1 = this version)
    //
    //     J0 (int) : angular momentum of operator
    //
    //     g0 (int) : parity grade of operator [P0=(-)^g0] (g0=0,1)
    //
    //     T0_min T0_max (int) : range of operator isospin components (a
    //       two-body operator may in general have components T0=0,1,2
    //       from the coupling of four isospin-1/2 fermionic operators)
    //
    //     symmetry_phase_mode (int) : RESERVED to describe how to
    //       obtain phase for lower triangle (see "Conjugation symmetry"
    //       below); currently only symmetry_phase_mode=0 is defined
    //
    //       Code note: An enum type basis::SymmetryPhaseMode is defined
    //       for this field, with kHermitian=0.
    //
    //     Nmax (int) : oscillator truncation of two-body space (Nmax>=0)
    //
    // Data
    //
    //   Data lines are of the form:
    //
    //     T0  Nr' lr' Nc' lc' L' S' J' T' g'  Nr lr Nc lc L S J T g  JT-RME
    //
    //   Although the g label is redundant (it can be deduced from lr and
    //   lc), it is included to make the sector structure more easily
    //   apparent to a human reader.
    //
    //   Here JT-RME is the JT-reduced matrix element under group
    //   theory conventions (i.e., no dimension factor in the
    //   Wigner-Eckart theorem):
    //
    //      < Nr' lr' Nc' lc' L' S' J' T' || op || Nr lr Nc lc L S J T >
    //
    //   Iteration follows the usual scheme within the basis module: sectors are
    //   lexicographic by (bra,ket) subspace indices, then matrix elements within
    //   a sector are lexicographic by (bra,ket) state indices.
    //
    //   It is assumed that the RelativeCMSectorsLSJT was constucted with the
    //   direction=kCanonical option.
    //
    //   Note: The concept of AS or NAS matrix elements does not apply in
    //   relative-cm states, so no normalization conversion is needed.
    //
    // For more information...
    //
    //   See the documentation of the *relative* LSJT operator file format for
    //   further general discussion: "Iteration order and symmetry", "Conjugation
    //   symmetry", and "Radial oscillator phase convention".  These discussions
    //   are generic to the LSJT scheme and apply the same to relative or
    //   relative-cm operators.
