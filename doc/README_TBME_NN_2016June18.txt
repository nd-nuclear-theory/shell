
	CONVENTIONS for TBME_NN input format (coupled J, not coupled T)
	see src_common/subrts_TBME_IO.f for source code  -- TO BE IMPLEMENTED

-----------------------------------------------------------------------------

	two-body nucleon-nucleon coupled J matrix elements 
	for identical or distinguishable particles preserving J and P
	(anticipating future extension to operators 
	that do not preserve J and/or P)

-----------------------------------------------------------------------------

 Norb			Number of single-particle orbitals (input)
 Jmax			maximum J value for 2-body system in MB-basis
 WT2max			maximum Weight for 2-body system in MB-basis

 Ntps			Number of (antisymmetrized) 2-body coupled-J states
 			(depends on Norb, Jmax, and WTmax)
 Ntbme			Number of coupled-J TBME's
 			(depends on Norb, Jmax, and WTmax)
 

Arrays needed for retrieval of desired matrix element

 ntpsJ(0:Jmax+1, 0:1)	number of (antisymmetrized) 2-body coupled-J states
 		 	for each J and IP = 0, 1
			Note: ntpsJ(Jmax+1, 0:1) = 0 
			      Ntps = SUM(ntpsJ)

 ntbmeJ(0:Jmax+1, 0:1)	number of TBME's up to this J and IP = 0, 1
 		  	i.e. the offset in the TBME array for each J and IP
 		  	Note: ntbmeJ(0,0) = 0
			      ntbmeJ(Jmax+1, 0:1) = Ntbme

 tpsJindx_IDN(0:Jmax+1, Norb_a * (Norb_a + 1)/2)
		index for antisymmetrized 2-body coupled-J state A(a, a')_J
		convention: a =< a'

 tpsJindx_DIS(0:Jmax+1, Norb_a, Norb_b)
		index for 2-body coupled-J state (a, b)_J


Array needed for (ascii) reading and/or processing of input TBME's

 TPJstates(2, Ntps)	 orbital indices for 2-body coupled-J states
(more complete would be
 TPJstates(3, Ntps) orbital indices and J values for 2-body coupled-J states)


NOTE: The SPorbitals have to be included in each TBME_NN file !!!

NOTE: If the SPorbitals are NOT YET defined when the TBME_NN file is 
      read in, MFDn will read and define the the orbitals from the TBME_NN file

NOTE: Ideally, for MFDn, the SPorbitals included in each TBME_NN file
      are exactly the ones used in the Many-Body basis;
      More SPorbitals than necessary is allowed the FIRST time, 
      and the FIRST TIME ONLY, that a TBME_NN file is read in ;
      MFDn will stop / crash / produce nonsense if there are less
      (or different) SPorbitals than necessary in the TBME_NN file.

NOTE: If the truncation parameters (Jmax and/or WT2max) are too small
      for the MB-basis, MFDn will stop / crash / produce nonsense ;
      If the truncation parameters (Jmax and/or WT2max) are larger 
      than necessary for the MB-basis, then the FIRST time, 
      and the FIRST TIME ONLY, that a TBME_NN file is read in, 
      Jmax and/or WT2max are increased to match the TBME_NN file.

NOTE: ALL TBME_NN files used in the same MFDn run HAVE to have 
      the same SPorbitals and the same truncations 
      (i.e. the same Jmax and WT2max etc)

NOTE: The parity of the coupled J TBME matrix elements is (currently)
      NOT included in the input TBME_NN files; the parity follows
      directly from the orbital angular momentum l of the SPorbitals.

-----------------------------------------------------------------------------

Supported TBME_NN input formats

-----------------------------------------------------------------------------

	ASCII format

# lines starting with # at the beginning of the file are not interpreted
# so comments can be put at the top of the file, 
# e.g. describing the contents:
# Version number  --   same as in MFDn_SP_Orbitals
# followed by number of orbitals for two classes of fermions,
# followed by Norb_a sets of orbital labels for protons (cls = 1)
# followed by Norb_b sets of orbital labels for protons (cls = 2)
# each set of orbital labels (n, l, 2j, cls, wt) on new line, 
# separated by one or more space(s)
#
# followed by 2*J_operator and P_operator (0 for pos, 1 for neg)
#   (currently only scalar P=0 operators supported/implemented
#   but anticipating future nonscalar parity-changing operators)
#
# followed by WTmax for class 1 (a) and class 2 (b) orbitals
#   NOTE: currently not implemented, 
#         but implicit through orbital input and/or WT2max
#
# followed by WT2max for aa, bb, and ab two-body coupled-J states,
#   NOTE: currently actually only implemented/tested 
#         for WT2max = max(WT2max_pp,  WT2max_nn,  WT2max_pn)
#
# followed by 2*Jmax for aa, bb, and ab two-body coupled-J states,
#   NOTE: currently actually only implemented/tested 
#         for Jmax = max(Jmax_pp,  Jmax_nn,  Jmax_pn)
#
# followed by number of TBMEs for aa, bb, and ab two-body coupled-J states
#
# ordering convention (though not strictly necessary for ASCII input)
#    ia =< ib ; ic =< id ; (ia, ib) =< (ic, id)
#

 15099

<  MFDn_SP_Orbitals.info  >

  2*J_operator   P_operator
  WT1max_p    WT1max_n
  WT2max_pp   WT2max_nn   WT2max_pn
  2*Jmax_pp   2*Jmax_nn   2*Jmax_pn
  Ntbme_pp    Ntbme_nn    Ntbme_pn

  ia   ib   ic   id   2*Jab  2*Jcd   11   matel_pp

  ia   ib   ic   id   2*Jab  2*Jcd   22   matel_nn

  ia   ib   ic   id   2*Jab  2*Jcd   12   matel_pn


NOTE: In the ascii format the order of the TBMEs within a class
      is NOT relevant -- as long as there are first Ntbme_pp pp 
      matrix elements, followed by Ntbme_nn nn matrix elements, 
      followed by Ntbme_pn pn matrix elements.

------------------------------------------------------------------

	BINARY format

  versionnumber = 15099
  Norb_p, Norb_n

  n_orb(1:norb_p)
  l_orb(1:norb_p)
  j2_orb(1:norb_p)
  wt_orb(1:norb_p)

  n_orb(1:norb_n)
  l_orb(1:norb_n)
  j2_orb(1:norb_n)
  wt_orb(1:norb_n)

  2*J_operator   P_operator
  WT1max_p    WT1max_n
  WT2max_pp   WT2max_nn   WT2max_pn
  2*Jmax_pp   2*Jmax_nn   2*Jmax_pn
  Ntbme_pp    Ntbme_nn    Ntbme_pn

  matel_pp(1:Ntbme_pp)
  matel_nn(1:Ntbme_nn)
  matel_pn(1:Ntbme_pn)

-----------------------------------------------------------------------------
