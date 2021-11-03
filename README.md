# Path Integral Monte Carlo


Path integrals and elementary lattice QCD using Monte Carlo techniques.

### CONTENTS

`UtilityFunctions.hpp` contains a handful of functions which help handle and print std containers.

`DynamicalSystem.hpp` contains definitions of types suited to represent (a limited class of) classical dynamical systems.

For the moment, the following classes are defined:
- `DynamicalException` (template class)
- `DynamicalSystem` (template class)
- `EuclidFreeParticle1D`
- `EuclidHarmonicOscillator1D`
- `EuclidParticle1D`

Methods are implemented to set paths in coordinate space, to fix boundary conditions and to compute the classical motion (for the moment, only for free particle and harmonic oscillator). Support for handling random paths is included. Functions that compute the discretized action will be defined separately. The exact Feynman amplitude is also provided whenever it is known.

The project is still at the early stage...

**Requires c++17**





