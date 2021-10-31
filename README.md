# Path Integral Monte Carlo


Path integrals and elementary lattice QCD using Monte Carlo techniques.

CONTENTS:
UtilityFunctions.hpp contains a handful of functions which help manipulate and print std containers
DynamicalSystems.hpp contains definitions of types suited to represent (a limited class of) classical dynamical systems. For the moment, the following classes are defined:
  DynamicalException (template class, contains a pointer to the DynamicalSystem which raised the exception)
  DynamicalSystem (template class)
  SolvableDynamicalSystem (template class)
  HarmonicOscillator1D
  FreeParticle1D
Methods are implemented to set paths in coordinate space, to fix boundary conditions and to compute the classical motion when an exact solution exists. When an action is defined, a discretized version can be conputed on the given path, for use in the path integral.

The project is still at the early stage...




