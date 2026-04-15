
* identity: identity operator (scalar)

    + Generated with obmixer.

* quadrupole: quadrupole operator

    + Generated with obmixer.
    
    + Obtained as a "solid harmonic" (see docstring for
      SolidHarmonicOneBodyOperator, in obutils/obme.h):
    
      Q2 = r^2 * Y2
         = [5/(4*pi)]^(1/2) * r^2 * C2
      
    + For illustration, we also construct it from the product of r.r and Y2,
      constructed separately.
    


