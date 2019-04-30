# MFDn H2 data format notes (v15200) #

Patrick J. Fasano, Mark A. Caprio

  + 02/12/19 (pjf): Created, duplicated from h2v15099-format-description-180425.md.
  + 02/21/19 (pjf): Updated, with minor revisions from pm.
  + 04/07/19 (pjf): Updated, with additional revisions from pm and mac.
  + 04/25/19 (mac): Typographical and formatting corrections.
----------------------------------------------------------------

## Text vs. binary ##

Both text and binary formats are defined.  The header information and matrix
element ordering are identical.  The text version is meant for human-readability
(for diagnostic purposes) and for quick-and-dirty interchange with other codes.
The binary format is much more compact and provides for faster I/O and thus
should be used for production purposes with large-scale runs.

The I/O routines in the shell project and MFDn assume that a file with extension
`.dat` is in text format, and that a file with extension `.bin` is in binary format.

The data stored in binary format are essentially identical.  We thus first
describe the text format in detail, then, for binary format, we need only
describe the technicalities of binary representation.

## File formatting: Text ##

Here we summarize the syntax for each part of the text-format h2 file.  The
examples are taken from:

    shell/doc/h2/h2v15200/example/identity/identity_Nmax04_h2v15200.dat

### Version number ###

Syntax:

    version_number

Example:

    15200

### Orbital listing header ###

Syntax:

    num_orbitals_p num_orbitals_n

Example:

    15 15

### Orbital listing body ###

Syntax:

     index  n  l  twice_j  twice_tz  weight

  * index (int): 1-based index for orbital
  * n (int): radial n quantum number
  * l (int): orbital l quantum number
  * twice_j (int): twice j quantum number
  * twice_tz (int): twice tz quantum number (+1 for protons, -1 for neutrons)
  * weight (float): weight in weight-max many-body truncation

Order:
  Lexocographic in (weight, -tz, l, j, n)

Example:

       1   0   0   1   1   0.00000000
       2   0   1   1   1   1.00000000
       ...
      15   0   4   9   1   4.00000000
      16   0   0   1  -1   0.00000000
      17   0   1   1  -1   1.00000000
       ...
      30   0   4   9  -1   4.00000000

### Header line 1: Operator properties ###

Syntax:

    J0 g0 Tz0

  * J0 (int): J for operator
  * g0 (int): parity grade for operator (g=0,1; P=(-)^g)
  * Tz0 (int): Tz for operator (HEP "up quark is up" convention, i.e., positive for proton)
      Note: Tz0 MUST be greater than or equal to zero.

Example:

For a scalar, positive parity, Tz-conserving operator...

    0 0 0

### Header line 2: 2-body basis limit ###

These are the maximum weights for truncation on each 2-body subspace (pp, pn, nn).

Syntax:

    wpp wpn wnn

  * wpp (float): weight max for two-body pp state
  * wpn (float): weight max for two-body pn state
  * wnn (float): weight max for two-body nn state

Example:

    4.000000e+00 4.000000e+00 4.000000e+00

### Header line 3: Matrix size ###

Total number of TBMEs within each 2-body *sector* (pp=+1, pn=0, nn=-1).  For
Tz-changing operators, this is specifically the pp/nn/pn label for the
*bra*.

Syntax:

    size_pp size_pn size_nn

Example:

    481 1856 481

### Matrix elements ###

In text format, indexing information is stored along with the matrix element.
This is purely for readability and diagnostic pursposes.  The ordering of matrix
elements is well defined, so the labels for each matrix element can be deduced
purely from the information in the header.

Syntax:

    i1 i2 i3 i4 twice_J_bra twice_J_ket matrix_element

  * i1 i2 i3 i4 (int): single particle orbital indices
  * twice_J_bra twice_J_ket (int): 2*J for bra and ket
  * matrix_element (float): reduced matrix element (RME)

    * This is the RME <i1 i2 J_bra || T || i3 i4 J_ket>.

    * The RME follows the group theory (=Rose) normalization and phase
    convention for the Wigner-Eckart theorem.  Thus, for scalar operators, it is
    simply the M-independent matrix element

          <i1 i2 J || T_0 || i3 i4 J> = <i1 i2 J M || T_0 || i3 i4 J M>

    This matrix element differs by a factor of J_bra-hat from the RME under the
    Edmonds (=Suhonen=Varshalovich) convention for the Wigner-Eckart theorem.

    * The bra and ket in the RME are *normalized* antisymmetrized J-coupled
      two-body states.

Example:

    1   1   1   1   0   0   +1.0000000e+00
    1   1   1   4   0   0   +0.0000000e+00
    1   1   1  11   0   0   +0.0000000e+00
    1   1   2   2   0   0   +0.0000000e+00
    1   1   2   7   0   0   +0.0000000e+00
    1   1   3   3   0   0   +0.0000000e+00
    1   1   3   8   0   0   +0.0000000e+00
    1   1   4   4   0   0   +0.0000000e+00
    1   1   5   5   0   0   +0.0000000e+00
    1   1   6   6   0   0   +0.0000000e+00
    1   4   1   4   0   0   +1.0000000e+00
    1   4   1  11   0   0   +0.0000000e+00
    ...

## File formatting: binary ##

### Data types ###

In binary format, the basic data content is the same as in text mode.  All
integers are stored as 4-byte integers, and all floating point values are stored
as 4-byte IEEE single-precision floating point values.

### FORTRAN record structure ###

The complication is that MFDn uses standard FORTRAN binary I/O, which
is "record" based.  Each record is sandwiched between two delimiters.
These are 4-byte integers, each giving the byte-count of the record.

So, how is the file broken into records?  To ensure forward compatibility, we
write the format code as it own record.  Then the rest of the header is broken
into records roughly corresponding to one record per header line (see below for
details).  Then the rest of the matrix elements are written out in three
records: all the pp matrix elements, all the nn matrix elements, and all the pn
matrix elements.

FORTRAN handles the record structure automatically.  But codes in
other languages (C++, etc.) can add in these record delimiters manually.

### Header data ###

The header data fields are grouped into FORTRAN writes as follows:

  * Version number

  * Orbital listing header

  * A separate array is written for each orbital property:

    * n
    * l
    * twice_j
    * twice_tz
    * weight

  * Header line 1: operator properties

  * Header line 2: 2-body basis limit

  * Header line 3: matrix size


### Matrix elements ###

In binary format, the matrix elements are simply written one after another
without any accompanying indexing information.  All that indexing information is
superfluous, since we have a well-defined ordering for the matrix elements.

The matrix elements are grouped into three FORTRAN writes, for the:

  * pp sector
  * pn sector
  * nn sector

### Example ###

Let consider the binary format file corresponding to the above example.  The
binary file is

    shell/doc/h2/h2v15200/example/identity/identity_Nmax04_h2v15200.bin

We will examine its raw binary contents using the "octal dump" utility od.
First, here is a raw hex dump, for reference

    % od -x identity_Nmax04_h2v15200.bin

    0000000 0004 0000 3b60 0000 0004 0000 0008 0000
    0000020 000f 0000 000f 0000 0008 0000 003c 0000


Let us now try to parse this as 4-byte integers

    % od -Ad -t d4 identity_Nmax04_h2v15200.bin

    0000000           4       15200           4           8
    0000016          15          15           8          60
    ...

Thus, you can see that the "Version number" record consists of the initial byte
count 4, followed by the verion number itself 15200, followed by the terminal
byte count 4.  Then "Orbital listing header" record consists of the byte count
8, followed by the proton and neutron orbital counts 15 and 15, followed by the
terminal byte count 8.  And so on...

To examine the matrix elements, we can re-parse the file as 4-byte floats, using

    % od -Ad -t fF identity_Nmax04_h2v15099.bin

    ...
    0000672               0               1               0               0
    0000688               0               0               0               0
    0000704               0               0               1               0
    0000720               0               0               0               0
    ...

## Matrix elements: ordering ##

### Overview ###

   The ordering of sectors, and then of matrix elements within a sector, is
governed by the indexing scheme for the two-body basis, which is defined in the initial
comments in `basis/jjjpn_scheme.h`.

Then, to translate this indexing for the basis into an ordering for sectors and
matrix elements, the two basic rules to remember are that we always

   (1) run in lexicographic order by (bra,ket), i.e., row-major order,
   and

   (2) restrict attention to the "upper triangle".

### Indexing listing ###

You can also get detailed listings of the indexing for various h2 files by
running h2stat on them, with the "--indexing" option.  For example:

   ~~~~sh
   cd shell/doc/h2/h2v15200/example/zero
   h2stat --indexing zero-e0_Nmax04_h2v15200.dat
   h2stat --indexing zero-m1_Nmax04_h2v15200.dat
   ~~~~

### Sectors ###

   The outermost "loop" is over sectors, i.e., over blocks of the matrix.  These
are pairs of two-body subspaces: ((bra_species,bra_J,bra_g),(ket_species,bra_J,bra_g)

   So, first off, we have to know the ordering of two-body subspaces.  Subspaces
are labeled by

    (Tz,J,g)

where P=(-)^g, and species is in [+1,0,-1].

  Ordering of subspaces is then lexicographic (-Tz,J,g), that is

    # two-body subspace ordering
    for Tz in [+1,0,-1]
      for J=0..Jmax(species)
        for g in [0,1]

where Jmax(species) is the maximum J allowed by the two-body truncation and
antisymmetry.

   Now to traverse the sectors, we must just enumerate pairs of subspaces, which
we traverse in lexicographic order by (bra_subspace,ket_subspace).  We only
include sectors in the "upper triangle" of the operator matrix, that is, where
the subspace indices are in canonical order (bra_subspace_index<=ket_subspace_index).
Matrix elements in the remaining sectors can be deduced by symmetry.  This sector
ordering can be enumerated by a loop of the form:

    # sector ordering
    for each bra_subspace_index
      for each ket_subspace_index >= bra_subspace_index
         if selection rules satisfied
            append sector (bra_subspace_index,ket_subspace_index)

   The selection rules between subspaces depend on the quantum numbers of the
operator (Tz0,J0,g0):

   * Tz selection:

         ket_Tz = bra_Tz - Tz0

     Note that, for a two-body operator, Tz0 is restricted to be 0,+1,+2.

   * A momentum selection:

         triangle(ket_J,J0,bra_J)

   * Parity selection:

         ket_g = g0 + bra_g (modulo 2)

### Matrix elements within a sector ###

   Within each sector (i.e., block of the matrix), we write out the matrix
elements in lexicographic order by (bra_index,ket_index).  Equivalently, we
write them out in row-major order.

   In "diagonal" sectors, i.e., where the bra and ket subspaces are the same, we
only write out matrix elements with indices in canonical order
(bra_index<=ket_index).  Equivalently, we only write out the upper triangle of
that block.

    # matrix element ordering (within a sector)
    for each bra_index
      for each ket_index
        if (bra_sector_index!=ket_sector_index) or (bra_index<=ket_index)
          write (bra_index,ket_index)

