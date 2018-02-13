
Description of format for input/output of Single-Particle orbitals:

Lines starting with # at the top of the file are not interpreted,
but can be used to put comments at the beginning of the file.

The actual content of the SPorbital file should be:
Version number (the current version number is 15055),
followed by the number of orbitals for two classes (cls) of fermions,
followed by Norb_1 sets of orbital labels for cls = 1 (protons)
followed by Norb_2 sets of orbital labels for cls = 2 (neutrons).

The set of orbital labels are (index, n, l, 2*j, cls, wt),
and their values should be space-delimited.
The index runs from 1 to Norb_i,
the labels (n, l, 2*j) are integers >= 0,
and the weights wt are real numbers >= 0.

Each set of orbital labels should be on a new line.

The order of the (Norb_1 + Norb_2) sets of orbital labels
is irrelevant (at least in the ascii format).

There is no 'end-of-file' symbol; additional lines beyond
the (Norb_1 + Norb_2) sets of orbital labels are ignored.

 Versionnumber
 Norb_1  Norb_2

 i  ni  li  j2i  clsi  wti    (for i = 1, Norb_1) 

 i  ni  li  j2i  clsi  wti    (for i = 1, Norb_2) 

------------------------------------------------------------------------

#
# Below is an example for an SPorbital file 
# with 0s1/2, 0p1/2, and 0p3/2 for protons
# and 0s1/2, 0p1/2, 0p3/2, 1s1/2, 0p3/2, and 0p5/2 for neutrons
# with weights given by (n+l)
#

 15055
 3   6

 1  0  0  1  1  0.0
 2  0  1  3  1  1.0
 3  0  1  5  1  1.0
 
 1  0  0  1  2  0.0
 2  0  1  3  2  1.0
 3  0  1  5  2  1.0
 4  1  0  1  2  1.0
 5  0  2  3  2  2.0
 6  0  2  5  2  2.0

------------------------------------------------------------------------

Notes on the IO and use of Single-Particle orbitals in MFDn
(see src_common/subrts_IO.f for source code -- TO BE IMPLEMENTED/UPDATED)

In MFDn the total number of SPorbitals is Norb = Norb_1 + Norb_2

The orbital indices for the second 'class' (or 'type') of fermions 
start at Norb_1 + 1 and run up to Norb_1 + Norb_2 

Arrays with non-negative INTEGER orbital labels (quantum numbers)

 n_orb(1:Norb)		radial quantum numbers
 l_orb(1:Norb)		orbital quantum numbers
 j2_orb(1:Norb)		twice j-value

 cls_orb(1:Norb)	'class': 1 or 2

Array with non-negative REAL orbital labels
(used for single-particle and/or many-body truncation)

 wt_orb(1:Norb)		weight factor (user-defined, default 2n+l )


NOTE:
If necessary, the orbitals are re-ordered within MFDn, such that
within MFDn orbital weights are monotonically increasing (or more
precisely, not decreasing) with the orbital index within a 'class'.

The output of MFDn one-body (and two-body density matrix elements)
uses the input-indexing of the orbitals.  If the orbitals are
generated in MFDn (instead of read in from file), the order is
increasing in N = 2n+l, starting from N = n = l = 0
for orbitals within the same shell (with same N):
increasing in j, starting with 2*j = 1 up to 2*j = 2*N+1

QUESTION: which indexing of orbitals to use for smwf?
input-indexing or MFDn-internal indexing?

------------------------------------------------------------------------
