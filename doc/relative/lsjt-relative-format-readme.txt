Notre Dame lsjt relative interaction format

URL: https://github.com/nd-nuclear-theory/shell/doc/relative/lsjt-relative-format-readme.txt

History:
  10/28/17 (mac): Extract descriptions from basis commit 9d7b487.
  12/20/17 (mac): Extract all examples to relutils-examples.txt.
----------------------------------------------------------------

The lsjt relative interaction file format is defined in header file comments:

  - The Notre Dame lsjt relative operator file format is defined in

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
  // relative states in LSJT scheme
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (L,S,J,T,g)    P=(-)^g
  //
  //   L (int): orbital angular momentum of *relative* motion (=lr)
  //   S (int): total spin
  //   J (int): total angular momentum of *relative* motion (=Jr)
  //            (i.e., L coupled to S)
  //   T (int): isospin
  //   g (int): grade (=0,1) for the parity P of *relative*
  //            motion (=gr)
  //
  // state labels within subspace: (N)
  //
  //   N (int): oscillator quanta of relative motion (=Nr)
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by:
  //   -- increasing L (L=0,1,...,Nmax)
  //   -- increasing S (S=0,1)
  //   -- increasing J
  //   -- [T forced by L+S+T~1]
  //   -- [g forced by g~L]
  // subject to:
  //   -- triangularity of (L,S,J)
  //   -- parity constraint L~g
  //   -- antisymmetry constraint L+S+T~1
  //
  // Subspaces are asserted to have nonzero dimension (as a sanity
  // check).
  //
  // Note that ordering of subspaces is lexicographic by (L,S,J).
  //
  // Truncation of the space is by the relative Nmax.
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, the states are ordered by:
  //   -- increasing N
  // and subject to:
  //   -- oscillator branching constraints N>=L and N~L (or,
  //      equivalently, parity constraint N~g)
  //
  // This basis is for *identical* particle states, but the
  // antisymmetry constraint is already applied at the level of
  // selecting the subspace labels L+S+T~1.
  //
  ////////////////////////////////////////////////////////////////

Operator file format description (lsjt_operator.h):

  ////////////////////////////////////////////////////////////////
  //
  // general relative two-body operator file format
  //
  ////////////////////////////////////////////////////////////////
  //
  // Header
  //
  //   Header format:
  //
  //     # RELATIVE LSJT
  //     # ...
  //     version                                 # version
  //     J0 g0 T0_min T0_max symmetry_phase_mode # operator tensor properties
  //     Nmax Jmax                               # relative basis truncation
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
  //     Nmax (int) : oscillator truncation of relative space
  //     (Nmax>=0)
  //
  //     Jmax (int) : additional relative angular momentum truncation
  //       of relative space (Jmax<=Nmax+1); Jmax=Nmax+1 includes the
  //       full relative space at the given oscillator truncation, but
  //       operators may often be truncated at lower partial waves
  //
  // Data
  //
  //   Then data lines are of the form:
  //
  //      T0 N' L' S' J' T' N L S J T JT-RME
  //
  //   Here JT-RME is the JT-reduced matrix element under group
  //   theory conventions (i.e., no dimension factor in the
  //   Wigner-Eckart theorem):
  //
  //      < N' L' S' J' T' || op || N L S J T >
  //
  //   For the special case of a rotational scalar (J0=0), isoscalar
  //   (T0=0) operator, this is equivalently the unreduced matrix
  //   element, as is often used to represent Hamiltonians in shell
  //   model codes:
  //
  //      < N' L' S' J' T'; MJ MT | op | N L S J T; MJ MT >
  //
  //   (Note that this matrix element is independent of MJ and MT for a
  //   scalar, isoscalar operator.)
  //
  // Iteration order and symmetry
  //
  //   It is assumed that the RelativeSectorsLSJT was constucted with
  //   the direction=kCanonical option.
  //
  //   Then ordering is...
  //
  //     for each T0 component (T0=T0_min...T0_max)
  //       for each sector
  //         for each matrix element
  //
  //   Sectors are thus ordered lexicographically by
  //   (bra_subspace_index,ket_subspace_index) and subject to these
  //   indices (bra_subspace_index,ket_subspace_index) being in
  //   canonical order.  It is assumed that the matrix elements in the
  //   other direction can be obtained by symmetry.
  //
  //   States are likewise ordered lexicographical by
  //   (bra_state_index,ket_state_index), i.e., row-major ordering.
  //   In diagonal sectors (i.e., when the bra and ket subspaces are
  //   the same), only matrix elements with indices in canonical order
  //   (upper triangle) are read from or written to file.
  //
  // Conjugation symmetry
  //
  //   Since only the upper triangle of the operator is stored, the
  //   user code is responsible for obtaining the lower triangle
  //   (conjugate) matrix elements by "symmetry".  Recall that the
  //   matrix elements are JT-reduced matrix elements, not simple
  //   matrix elements.  Therefore, conjuation will in general involve
  //   phase and dimension factors.
  //
  //   The symmetry_phase_mode field in the header is reserved to
  //   provide information on the correct form to use for this phase
  //   factor.  However, for now, only the placeholder value
  //   0 (=kHermitian) is defined for symmetry_phase_mode, and phase
  //   conventions are only well-defined for a Hamiltonian-like (J0=0,
  //   g0=0) operator.
  //
  //   symmetry_phase_mode=0 (=kHermitian), J0=0, g0=0: For a
  //   Hamiltonian-like operator, we expect Hermiticity, i.e.,
  //   symmetry of (M_J,M_T)-branched matrix elements.  Within a
  //   diagonal sector, this means that the lower triangle is obtained
  //   from the upper triangle by ordinary symmetry.  For off-diagonal
  //   sectors, the appropriate symmetry on the JT-reduced matrix
  //   elements in general includes phase and dimension factors from
  //   the isospin Clebsches:
  //
  //     <a,J,T,g || A_{T0} || a',J,T',g>
  //       = (-)^(T'-T)*Hat(T')/Hat(T)
  //         * <a',J,T',g || A_{T0} || a,J,T,g>
  //
  //   However, note that the symmetry factor only differs from unity
  //   in the isospin-changing <T=0|T=1> sectors.  These sectors
  //   conventionally only contain vanishing matrix elements for
  //   nuclear interactions.
  //
  //   For an operator which transforms under conjugation as a
  //   spherical harmonic, if we look purely at the angular momentum
  //   dependence, we in general obtain the conjugated RME (under
  //   group theory conventions for the RME)
  //
  //          <J||A_{J0}||J'> = (-)^(J'-J)*Hat(J')/Hat(J)*<J'||A_{J0}||J>
  //
  //   This relation applies to both the M1 and E2 operators under
  //   Condon-Shortley phase conventions (Suhonen Ch. 6) for these
  //   operators.  It is derived from the W-E theorem, symmetry of CG
  //   coefficient, and conjugation properties of these operators.
  //
  //   Putting this together, for a JT basis, a spherical
  //   harmonic-like operator conjugates as
  //
  //     <a,J,T,g || A_{J0,T0} || a',J,T',g>
  //       = (-)^(J'-J)*Hat(J')/Hat(J)
  //         * (-)^(T'-T)*Hat(T')/Hat(T)
  //         * <a',J,T',g || A_{J0,T0} || a,J,T,g>
  //
  //   Note that Hermitian conjugation behavior is recovered in the
  //   special case of J0=0.
  //
  // Radial oscillator phase convention
  //
  //   Two phase conventions are in use for radial wave functions and
  //   thus for harmonic oscillator basis functions.  See, e.g.,
  //   footnote 8 of csbasis [PRC 86, 034312 (2012)].  We adopt the
  //   "positive at the origin" convention, which should thus be used
  //   when writing matrix elements to file in the present format.
  //   The "positive at the origin" convention is followed for the
  //   oscillator functions in, e.g., Suhonen (3.42), Moshinsky
  //   (1.1.8) & (1.9.15), and MFDn interaction files.  However, the
  //   "positive at infinity" convention is more natural when
  //   considering oscillator functions as members of SU(3) irreps.
  //   So it may be necessary to convert to this convention for
  //   internal use in SU(3) calculations.
