  Eigen::MatrixXd 
  KinematicUTSqrMatrixJJJPN(
      const basis::OrbitalSpaceLJPN& radial_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::MatrixVector& radial_matrices,
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int A
    );
  // Populate matrix with two-body matrix elements for the *one-body*
  // operator U_(T^2) obtained from a scalar one-body operator T^2,
  // i.e., r^2 or k^2, on a given JJJPN sector.
  //
  // This one-body operator is related to the two-body operator
  // V_(T^2), also defined in csbasis, by U_(T^2) = 1/(A-1)*V_(T^2).
  // See csbasis (51).
  //
  // Obtained by csbasis (52)-(54).
  //
  // Precondition: It is assumed that the sector is valid for a scalar
  // (and isoscalar and positive parity), i.e., is a diagonal sector.
  //
  // Arguments:
  //   radial_orbital_space, radial_sectors, radial_matrices (...):
  //      definition of T^2 radial matrix elements
  //   sector (basis::TwoBodySectorsJJJPN::SectorType) : The sector to
  //     populate.
  //   A (int): atomic mass number
  //
  // Returns:
  //   (Eigen::MatrixXd) : The matrix for this sector.

  Eigen::MatrixXd 
  KinematicVT1T2MatrixJJJPN(
      const basis::OrbitalSpaceLJPN& radial_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::MatrixVector& radial_matrices,
      bool momentum_space,
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int A
    );
  // Populate matrix with two-body matrix elements for the *two-body*
  // operator V_(T1.T2) obtained from a a dot product of one-body
  // vector operators (i.e., T = r or k), on a given JJJPN sector.
  //
  // Obtained by csbasis (55)-(60).
  //
  // Precondition: It is assumed that the sector is valid for a scalar
  // (and isoscalar and positive parity), i.e., is a diagonal sector.
  //
  // Arguments:
  //   radial_orbital_space, radial_sectors, radial_matrices (...):
  //      definition of T radial matrix elements
  //   momentum_space (bool): if need to include extra momentum space
  //     phase factor from csbasis (59)
  //   sector (basis::TwoBodySectorsJJJPN::SectorType) : The sector to
  //     populate.
  //   A (int): atomic mass number
  //
  // Returns:
  //   (Eigen::MatrixXd) : The matrix for this sector.
