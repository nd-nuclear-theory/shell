${SHELL_DIR}/bin/orbital-gen --oscillator 2 orbitals.dat
${SHELL_DIR}/bin/orbital-gen --oscillator 20 orbitals-int.dat
${SHELL_DIR}/bin/radial-gen --kinematic r 1 oscillator orbitals.dat radial-me-r1.dat
${SHELL_DIR}/bin/radial-gen --kinematic r 2 oscillator orbitals.dat radial-me-r2.dat
${SHELL_DIR}/bin/radial-gen --kinematic k 1 oscillator orbitals.dat radial-me-k1.dat
${SHELL_DIR}/bin/radial-gen --kinematic k 2 oscillator orbitals.dat radial-me-k2.dat
${SHELL_DIR}/bin/h2mixer < h2mixer.in