# generate radial matrix elements; input from stdin
../radialutils/radial-gen << EOF
set-ket-indexing test/orbitals-ho-Nmax20.dat
set-analytic-basis-type oscillator
define-operator-target r 1 1 test/ob-rme-ho-r1-Nmax20.dat
define-operator-target r 2 0 test/ob-rme-ho-r2-Nmax20.dat
define-operator-target k 1 1 test/ob-rme-ho-k1-Nmax20.dat
define-operator-target k 2 0 test/ob-rme-ho-k2-Nmax20.dat

define-xform-target 1.000 0 test/orbitals-ho-Nmax20.dat test/radial-olap-ho-b1.000-Nmax20.dat
define-xform-target 1.414 0 test/orbitals-ho-Nmax20.dat test/radial-olap-ho-b1.414-Nmax20.dat
EOF
