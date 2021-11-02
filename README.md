# Path Integral Monte Carlo


Path integrals and elementary lattice QCD using Monte Carlo techniques.

### CONTENTS

`UtilityFunctions.hpp` contains a handful of functions which help handle and print std containers.

`DynamicalSystems.hpp` contains definitions of types suited to represent (a limited class of) classical dynamical systems.

For the moment, the following classes are defined:
- `DynamicalException` (template class)
- `DynamicalSystem` (template class)
- `FreeParticle1D`
- `HarmonicOscillator1D`
- `Particle1D`

Methods are implemented to set paths in euclidean coordinate space, to fix boundary conditions and to compute the classical motion (for the moment, only for free particle and harmonic oscillator). If the action is defined, a discretized version can be computed on the given path, for use in the path integral. The exact Feynman amplitude is also provided whenever it is known.

The project is still at the early stage...




