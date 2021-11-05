# Path Integral Monte Carlo


Path integrals and elementary lattice QCD using Monte Carlo techniques.

### CONTENTS

`UtilityFunctions.hpp` contains a handful of functions which help handle and print std containers.

`DynamicalSystem.hpp` contains definitions of types suited to represent (a limited class of) classical dynamical systems.

`MonteCarlo.hpp` contains routines for Monte Carlo integration

For the moment, the following classes are defined:
- `DynamicalException` (template class)
- `DynamicalSystem` (template class)
- `EuclidFreeParticle1D`
- `EuclidHarmonicOscillator1D`
- `EuclidParticle1D`

The objective of these classes is to provide a unified and practical interface to set and use paths of dynamical systems in coordinate space. For the moment the simplest 1D case is developed most.

The project is still at the early stage...

**Requires c++17, boost and gsl**





