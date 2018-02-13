# MFDn H2 data format notes (v15099) #

Mark A. Caprio

2/12/18 (mac): Created.

----------------------------------------------------------------

## Overview ##

   The ordering of sectors, and then of matrix elements within a sector, is 
governed by the indexing scheme for the two-body basis, which is defined in the initial
comments in

   basis/jjjpn_scheme.h

Then, to translate this indexing for the basis into an ordering for sectors and
matrix elements, the two basic rules to remember are that we

   (1) always run in lexicographic order by (bra,ket), i.e., row-major order,
   and

   (2) for Tz-conserving operators, we restrict attention to the "upper triangle".

## Figure ##

Be sure to see the sketches in h2-sector-ordering_scan180122.pdf as you read the
following.

You can also get detailed listings of the indexing for various h2 files by
running h2stat on them, with the "--indexing" option.  For example:

   ~~~~
   cd shell/doc/h2/h2v15099/example/zero
   h2stat --indexing zero-e0_Nmax04_h2v15099.dat
   h2stat --indexing zero-m1_Nmax04_h2v15099.dat
   ~~~~

## Sectors ##

   The outermost "loop" is over sectors, i.e., over blocks of the matrix.  These
are pairs of two-body subspaces: ((bra_species,bra_J,bra_g),(ket_species,bra_J,bra_g)

   So, first off, we have to know the ordering of two-body subspaces.  Subspaces
are labeled by

    (species,J,g)

where P=(-)^g, and species is in [pp,nn,pn].

Ordering of subspaces is then lexicographic in these labels, that is

    # two-body subspace ordering
    for species in [pp,nn,pn]
      for J=0..Jmax(species)
        for g in [0,1]

   Now to traverse the sectors, we must just enumerate pairs of subspaces, which
we traverse in lexicographic order by (bra_subspace,ket_subspace).  For
Tz-conserving operators, we only include sectors in the "upper triangle" of the
operator matrix, that is, where the subspace indices are in canonical order
(bra_subspace_index<=ket_subspace_index).  Matrix elements in the remaining
sectors can be deduced by symmetry.  This sector ordering can be enumerated by a
loop of the form:

    # Tz0=0 sector ordering
    for each bra_subspace_index
      for each ket_subspace_index >= bra_subspace_index
         if selection rules satisfied
            append sector (bra_subspace_index,ket_subspace_index)

However, for Tz-changing operators, we cannot impose this restriction, as
canonical order would "miss" certain needed sectors (e.g., we would never see
the matrix elements from an nn subspace to a pn subspace, <pn|O|nn>).  We thus
must relax the canonical ordering constraint:

    # Tz0!=0 sector ordering
    for each bra_subspace_index
      for each ket_subspace_index
         if selection rules satisfied
            append sector (bra_subspace_index,ket_subspace_index)

   The selection rules between subspaces depend on the quantum numbers of the
operator (Tz0,J0,g0):

    - Tz selection:

          bra_Tz == Tz0 + ket_Tz

      Note that, for a two-body operator, Tz0 can be -2,-1,0,+1,+2.  To
      interpret this selection rule, it is important to know the sign convention
      for Tz: we follow the particle-physics convention that Tz=+1/2 for proton
      and -1/2 for neutron ("up quark is isospin up").  Thus, e.g., beta-minus
      decay turns a neutron into a proton, so it carries Tz0=+1.

    - angular momentum selection:

          triangle(ket_J,J0,bra_J)

    - parity selection:

          ket_g = g0 + bra_g (modulo 2)

## Matrix elements within a sector ##

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

   The ordering of matrix elements is *identical* in ASCII (.dat) and binary
(.bin) files.  The sector-by-sector I/O routines in shell depend on this.
(Don't believe any comments in Pieter's README files to the contrary.)

## More on header formats ##

   See Pieter's readme file

       shell/doc/h2/h2v15099/README_TBME_2016June20.txt

## More on FORTRAN record delimiters ##

   See my old notes on this in matter in

       shell/doc/h2/h2v15099/h2v15099-format-description-180212.txt
