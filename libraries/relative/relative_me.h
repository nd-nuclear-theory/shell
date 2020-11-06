/****************************************************************
  relative_me.h

  Construction of constant and kinematic (r^2 and k^2) relative operators.

  (The identity operator is provided in basis/lsjt_operator.h.)

  Note on oscillator length: Oscillator lengths for two-particle rel-cm
  coordinates are related to oscillator length for single-particle basis as:

    b_rel = 2^(1/2)*b
    b_cm = 2^(-1/2)*b

  The matrix elements generated by these routines are calculated for oscillator
  length parameter b_rel = 2^(1/2) in the relative coordinate, which, after
  Moshinsky transformation, yields matrix elements for two-body states built
  from single-particle wave functions with oscillator length parameter b=1.
  However, in the actual calculations, we first calculate matrix element for
  relative b_rel=1, then, at the end, we analytically rescale the result to
  b_rel = 2^(1/2).  Matrix elements of r^n must scale by (b_rel/b)^n, and those
  for k^n must scale by (b_rel/b)^-n.  Thus, e.g., for Coulomb (n=-1) the scale
  factor is 1/sqrt(2.).

  Note on relative coordinate convention: The relative coordinate and momentum
  are taken in "mechanics" convention

    r_rel = r_1 - r_2

    k_rel = (1/2) (k_1 - k_2)

  This is to be contrasted with the symmetrized "Moshinsky" convention.

  See notes on "internal representation of an operator in JT scheme" in
  lsjt_operator.h for the general principles of how the operators are
  represented.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 7/12/16 (mac): Created (as construct_relative).
  + 3/26/17 (mac): Add ConstructCoulombOperator.
  + 3/28/17 (mac):
    - Rename to relative_me.
    - Remove ConstructDiagonalConstantOperator.
  + 7/1/17 (mac): Add ConstructQuadrupoleOperator.
  + 10/19/17 (mac): Add scale factors for relative oscillator length.
  + 04/28/18 (mac): Add isoscalar L and S operators.
    operators.
  + 05/04/18 (mac):
    - Add isovector variants of operators (L, S, Q, r^2).
    - Rename enum KinematicOperator to CoordinateType.
    - Rename ConstructKinematicOperator to ConstructCoordinateSqr.
  + 05/05/18 (mac):
    - Clean up internal implementation of radial SU(1,1) matrices.
    - Implement dipole operator.
  + 08/08/18 (pjf): Use analytic and spline wrapper functions for radial me.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
  + 06/20/19 (pjf): Add isospin operator.
  + 10/06/19 (pjf): Fix sqrt(5/16pi) prefactor in ConstructQuadrupoleOperator.
  + 11/02/20 (pjf): Add SU(4) Casimir operator.
  + 11/06/20 (pjf): Add diagnostic output to ConstructCoulombOperator.

****************************************************************/

#ifndef RELATIVE_ME_H_
#define RELATIVE_ME_H_

#include <array>
#include "basis/lsjt_operator.h"
#include "basis/proton_neutron.h"

namespace relative {

  ////////////////////////////////////////////////////////////////
  // typedef for operator coordinate/momentum mode
  ////////////////////////////////////////////////////////////////

  enum class CoordinateType {kK,kR};

  ////////////////////////////////////////////////////////////////
  // quadratic operators in coordinates and momenta
  ////////////////////////////////////////////////////////////////

  void ConstructCoordinateSqr(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      relative::CoordinateType coordinate_type,
      int T0
    );
  // Construct simple "kinematic" (r^2 and k^2) operators in relative
  // LSJT basis.
  //
  // These matrix elements are for the relative (r_rel)^2 or (k_rel)^2 operators
  // on the two-body system.  The intrinsic r^2_intr or k^2_intr operators on the
  // A-body system are obtained as the two-body operators
  //
  //   r^2_intr = (1/A) V[(r_rel)^2]
  //   k^2_intr = (4/A) V[(k_rel)^2]
  //
  // and T_intr is obtained as
  //
  //   T_intr = (2/A) (hbar^2/(4m)) V[(k_rel)^2]
  //
  // or, in terms of the two-body relative kinetic energy, as
  //
  //   T_intr = (2/A) V[T_rel]
  //   T_rel = (hbar^2/(4m)) (k_rel)^2
  //
  // These matrix elements are calculated for oscillator length
  // parameter b_rel = 2^(1/2) in the relative coordinate.  (See "Note
  // on oscillator length" at start of this header file.)
  //
  // Makes use of analytic matrix elements for r^2 and k^2 in
  // oscillator basis.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (input) : tensorial properties of operator
  //   relative_space (input) : target space
  //   relative_component_sectors (output) : target sectors
  //   relative_component_matrices (output) : target matrices
  //   coordinate_type (input): whether to construct coordinate
  //     or momentum operator
  //   T0 (input): isoscalar (T0=1) or isovector (T0=1)

  void ConstructQuadrupoleOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      relative::CoordinateType coordinate_type,
      int T0
    );
  // Construct quadrupole operator in relative LSJT basis.
  //
  // These matrix elements are for the relative quadrupole operator Qrel on the
  // two-body system.  The intrinsic quadrupole operator Qintr on the A-body
  // system is obtained as the two-body operator
  //
  //   Qintr = ... V[Qrel]  [TODO fill in coefficient]
  //
  // See "intrinsic transition operators" pencilwork 5/4/18 for derivations and
  // discussion of isoscalar and isovector RMEs.
  //
  // These matrix elements are calculated for oscillator length
  // parameter b_rel = 2^(1/2) in the relative coordinate.  (See "Note
  // on oscillator length" at start of this header file.)
  //
  // Here the operator_type indicates protons/neutrons/both relative
  // to *total* center of mass.
  //
  // Arguments:
  //   operator_labels (input) : tensorial properties of operator
  //   relative_space (input) : target space
  //   relative_component_sectors (output) : target sectors
  //   relative_component_matrices (output) : target matrices
  //   coordinate_type (input): whether to construct coordinate
  //     or momentum operator
  //   T0 (input): isoscalar (T0=1) or isovector (T0=1)

  void ConstructOrbitalAMOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      int T0
    );
  // Construct orbital angular momentum operator in relative LSJT basis.
  //
  // These matrix elements are for the relative orbital angular momentum
  // operator Lrel on the two-body system.  The intrinsic orbital angular
  // momentum operator Lintr on the A-body system is obtained as the two-body
  // operator
  //
  //   Lintr = (2/A) V[Lrel]
  //
  // These matrix elements are independent of oscillator length.  (See "Note on
  // oscillator length" at start of this header file.)
  //
  // See "intrinsic transition operators" pencilwork 5/4/18 for derivations and
  // discussion of isoscalar and isovector RMEs.
  //
  // Arguments:
  //   operator_labels (input) : tensorial properties of operator
  //   relative_space (input) : target space
  //   relative_component_sectors (output) : target sectors
  //   relative_component_matrices (output) : target matrices
  //   T0 (input): isoscalar (T0=1) or isovector (T0=1)

  ////////////////////////////////////////////////////////////////
  // linear operators in coordinates and momenta
  ////////////////////////////////////////////////////////////////

  void ConstructDipoleOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      relative::CoordinateType coordinate_type,
      int T0
    );
  // Construct dipole (E1) operator in relative LSJT basis.
  //
  // These matrix elements are for the relative dipole operator Drel on the
  // two-body system.  The intrinsic dipole operator Dintr on the A-body system
  // is obtained as the two-body operator
  //
  //   Dintr = (1/A) V[Drel]
  //
  // These matrix elements are independent of oscillator length.  (See "Note on
  // oscillator length" at start of this header file.)
  //
  // Note that the isoscalar dipole operator (T0=0) is identically zero, and an
  // attempt to generate it will result in an assertion failure.  However, the
  // T0 parameter is still provided for uniformity of syntax.
  //
  // See "intrinsic transition operators" pencilwork 5/4/18 for derivations and
  // discussion of isoscalar and isovector RMEs.
  //
  // Arguments:
  //   operator_labels (input) : tensorial properties of operator
  //   relative_space (input) : target space
  //   relative_component_sectors (output) : target sectors
  //   relative_component_matrices (output) : target matrices
  //   relative::CoordinateType coordinate_type
  //   T0 (input): isovector (T0=1) ONLY

  ////////////////////////////////////////////////////////////////
  // spin operators
  ////////////////////////////////////////////////////////////////

  void ConstructSpinAMOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      int T0
    );
  // Construct spin angular momentum operator in relative LSJT basis.
  //
  // These matrix elements are for the spin angular momentum operator Srel taken
  // as a relative (Galilean-invariant) operator on the two-body system.  The
  // spin angular momentum operator S on the A-body system is obtained as the
  // two-body operator
  //
  //   S = 1/(A-1) V[Srel]
  //
  // These matrix elements are independent of oscillator length.  (See "Note on
  // oscillator length" at start of this header file.)
  //
  // See "intrinsic transition operators" pencilwork 5/4/18 for derivations and
  // discussion of isoscalar and isovector RMEs.
  //
  // Arguments:
  //   operator_labels (input) : tensorial properties of operator
  //   relative_space (input) : target space
  //   relative_component_sectors (output) : target sectors
  //   relative_component_matrices (output) : target matrices
  //   T0 (input): isoscalar (T0=1) or isovector (T0=1)

  ////////////////////////////////////////////////////////////////
  // isospin operators
  ////////////////////////////////////////////////////////////////

  void ConstructIsospinOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices
    );
  // Construct isospin operator in relative LSJT basis.
  //
  // These matrix elements are for the isospin operator Trel taken as a relative
  // (Galilean-invariant) operator on the two-body system.  The isospin operator
  // T on the A-body system is obtained as the two-body operator
  //
  //   T = 1/(A-1) V[Trel]
  //
  // These matrix elements are independent of oscillator length.  (See "Note on
  // oscillator length" at start of this header file.)
  //
  // Arguments:
  //   operator_labels (input) : tensorial properties of operator
  //   relative_space (input) : target space
  //   relative_component_sectors (output) : target sectors
  //   relative_component_matrices (output) : target matrices

  ////////////////////////////////////////////////////////////////
  // Wigner's SU(4) Casimir
  ////////////////////////////////////////////////////////////////

  void ConstructSU4CasimirOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices
    );
  // Construct SU(4) quadratic Casimir in relative LSJT basis.
  //
  // These matrix elements are for the SU(4) Casimir operator CSU4rel taken as a
  // relative (Galilean-invariant) operator on the two-body system.  The SU(4)
  // Casimir CSU4 on the A-body system is obtained as the two-body operator
  //
  //   CSU4 = V[CSU4rel] - (A-2)/(A-1) V[15/2 I]
  //
  // These matrix elements are independent of oscillator length.  (See "Note on
  // oscillator length" at start of this header file.)
  //
  // Arguments:
  //   operator_labels (input) : tensorial properties of operator
  //   relative_space (input) : target space
  //   relative_component_sectors (output) : target sectors
  //   relative_component_matrices (output) : target matrices

  ////////////////////////////////////////////////////////////////
  // Coulomb
  ////////////////////////////////////////////////////////////////

  void ConstructCoulombOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      basis::OperatorTypePN operator_type,
      int num_steps
    );
  // Construct Coulomb operator in relative LSJT basis.
  //
  // These matrix elements are calculated for oscillator length
  // parameter b_rel = 2^(1/2) in the relative coordinate.  (See "Note
  // on oscillator length" at start of this header file.)
  //
  // Here the operator_type indicates whether this is a single species
  // operator for protons/neutrons or sees both types of nucleons.
  //
  // Makes use of spline integration.  An offset inner integration
  // point x0>0 is required to prevent the integrand from exploding.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (input): tensorial properties of operator
  //   relative_space (input): target space
  //   relative_component_sectors (output): target sectors
  //   relative_component_matrices (output): target matrices
  //   operator_type (input): operator for protons, neutrons, or both
  //   num_steps (input): number of integration steps

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
