# Path Integral Monte Carlo


Path integrals and elementary lattice QCD using Monte Carlo techniques.

CONTENTS:

UtilityFunctions.hpp contains a handful of functions which help manipulate and print std containers\n
DynamicalSystems.hpp contains definitions of types suited to represent (a limited class of) classical dynamical systems. For the moment, the following classes are defined:\n
\t  DynamicalException (template class, contains a pointer to the DynamicalSystem which raised the exception)\n
\t  DynamicalSystem (template class)\n
\t  SolvableDynamicalSystem (template class)\n
\t  HarmonicOscillator1D\n
\t  FreeParticle1D\n
Methods are implemented to set paths in coordinate space, to fix boundary conditions and to compute the classical motion when an exact solution exists. When an action is defined, a discretized version can be conputed on the given path, for use in the path integral.

The project is still at the early stage...




