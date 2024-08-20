# Conversion from SPS+XPN to h2 #

Mark A. Caprio

  + 08/20/24 (mac): Created.
  
----------------------------------------------------------------

First, please read the initial program documentation in `sps2orbital.cpp` (under
`programs/obutils`) and `xpn2h2.cpp` (under `programs/h2utils`), for a general
understanding of the inputs and output of these programs.  Then, you can follow
the examples below.

In each case, we must first convert the BIGSTICK single-particle space (SPS)
file defining the relevant orbitals.  Then we can process the BIGSTICK explicit
proton-neutron (XPN) interaction file.  The numbering of the orbitals in the SPS
file *must* match the numbering assumed in the XPN file.


## Example: Conversion of shell-model interaction file ##

Here we convert the Wildenthal USD and the Brown & Richter USDB interactions,
for the sd shell.

We must first convert the BIGISTICK SPS file for the sd-shell orbitals
(`sd.sps`):

  ~~~~~~~~~~~~~~~~
  iso
  3
         0.0  2.0  1.5  2
         0.0  2.0  2.5  3
         1.0  0.0  0.5  2
  ~~~~~~~~~~~~~~~~
  
However, there is a complication.  In the last column, this SPS file contain
nonzero "weights" for the orbitals, meant to be used in shell model truncations.
For us, these weights would be a nuisance, as they would need to be taken into
account in defining weight-based truncation parameters when we create the h2
file for the interaction.  So we call `sps2orbital` with the option
`--clear-weights`, to force these weights to zero:

  ~~~~~~~~~~~~~~~~
  sps2orbital --clear-weights sd.sps sd_orbital.dat
  ~~~~~~~~~~~~~~~~

This generates the orbital file `sd_orbital.dat`:

  ~~~~~~~~~~~~~~~~
  # MFDn SPorbital file
  #   version
  #   norb_p norb_n
  #   index n l 2*j species weight
   15099
   3 3
     1   0   2   3   1   0.00000000
     2   0   2   5   1   0.00000000
     3   1   0   1   1   0.00000000
     1   0   2   3   2   0.00000000
     2   0   2   5   2   0.00000000
     3   1   0   1   2   0.00000000
  ~~~~~~~~~~~~~~~~
  
We are also given the BIGSTICK XPN file `w18pn.int`.  [Note: The `w` is for
"Wildenthal.  The 18 indicates that the TBMEs are those for the A=18 nuclei,
i.e., TBMEs for all other A are obtained by a standard (18/A)^(0.3) scaling
relation.  The `pn` indicates "proton-neutron" format.]

If we only wanted the two-body matrix elements (TBMEs) of the residual
interaction, we could simply run:

  ~~~~~~~~~~~~~~~~
  xpn2h2 sd_orbital.dat w18pn.int w-V.dat
  ~~~~~~~~~~~~~~~~

The results are stored in the h2 file (`w-V.dat`).

However, we also need the single-particle energies (SPEs), i.e., the mean-field part
of the Hamiltonian.  Thus, we run:

  ~~~~~~~~~~~~~~~~
  xpn2h2 --obme-filename w_obme.dat sd_orbital.dat w18pn.int w-V.dat
  ~~~~~~~~~~~~~~~~

The mean-field part of the Hamiltonian is stored here as one-body matrix
elements (OBMEs) in the one-body operator file `w_obme.dat`.

This one-body operator can be "upgraded" into an A-dependent two-body operator
[e.g., "intrinsic", JPG 47, 122001 (2020), DOI:10.1088/1361-6471/ab9d38, (9)].
We prepare the following input file for `h2mixer` (`h2mixer_w.in`):

  ~~~~~~~~~~~~~~~~
  # h2mixer_w.in -- upgrade SPEs to TBMEs for two-nucleon system
  
  set-target-indexing sd_orbital.dat 0. 0. 0. 0. 0.
  set-target-multipolarity 0 0 0
  set-output-format 15099
  set-mass 2
  
  # define input OBMEs and TBMEs
  define-ob-source input spes_ob w_obme.dat 0 0 0
  
  # define A-dependent two-body Hamiltonian
  define-tb-source operatorU spes_tb spes_ob
  define-target w-U.dat
    add-source spes_tb 1.0
  ~~~~~~~~~~~~~~~~
  
Then, we run:

  ~~~~~~~~~~~~~~~~
  h2mixer < h2mixer_w.in
  ~~~~~~~~~~~~~~~~
  
The TBMEs stored in `w-U.dat` are appropriate to the system with n=2 valence
particles (A=18).  They must be rescaled by 1/(n-1) for the system with n
valence particles (A=16+n).

We proceed similarly with the USDB interaction:

  ~~~~~~~~~~~~~~~~
  xpn2h2 --obme-filename usdb_obme.dat sd_orbital.dat usdbpn.int usdb-V.dat
  h2mixer < h2mixer_usdb.in
  ~~~~~~~~~~~~~~~~


## Example: Conversion of NCSM two-body interaction file ##

Here we take an XPN interaction file for the JISP16 interaction, in an
oscillator basis (with oscillator parameter hw=20 MeV), on two-body basis
truncated to 6 oscillator quanta (`tb-6` truncation).

We must first convert the BIGISTICK SPS file for the relevant (N<=6) oscillator
orbitals (`ncci-tb-6.sps`):

  ~~~~~~~~~~~~~~~~
  sps2orbital ncci-tb-6.sps ncci-tb-6_orbital.dat
  ~~~~~~~~~~~~~~~~
  
This generates the orbital file `ncci-tb-6_orbital.dat`.  Alternatively, we
could generate an identical list of oscillator orbitals directly, using
`orbital-gen`:

  ~~~~~~~~~~~~~~~~
  orbital-gen --Nmax 6 ncci-tb-6_orbital.dat
  ~~~~~~~~~~~~~~~~

Then, we convert the BIGSTICK XPN interaction file (`JISP16-tb-6-20.int`):

  ~~~~~~~~~~~~~~~~
  xpn2h2 --truncation 6 6 ncci-tb-6_orbital.dat JISP16-tb-6-20.int JISP16-tb-6-20.dat
  ~~~~~~~~~~~~~~~~

The results are stored in `JISP16-tb-6-20.dat`.


