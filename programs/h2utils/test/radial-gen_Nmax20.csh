time ../radialutils/radial-gen --kinematic r 1 oscillator test/orbitals-ho-Nmax20.dat test/radial-me-ho-r1-Nmax20.dat
time ../radialutils/radial-gen --kinematic r 2 oscillator test/orbitals-ho-Nmax20.dat test/radial-me-ho-r2-Nmax20.dat
time ../radialutils/radial-gen --kinematic k 1 oscillator test/orbitals-ho-Nmax20.dat test/radial-me-ho-k1-Nmax20.dat
time ../radialutils/radial-gen --kinematic k 2 oscillator test/orbitals-ho-Nmax20.dat test/radial-me-ho-k2-Nmax20.dat

time ../radialutils/radial-gen --overlaps 1.000 oscillator test/orbitals-ho-Nmax20.dat test/radial-olap-ho-b1.000-Nmax20.dat
time ../radialutils/radial-gen --overlaps 1.414 oscillator test/orbitals-ho-Nmax20.dat test/radial-olap-ho-b1.414-Nmax20.dat
