# GGG\_spin\_alignment

A code for examining spin canting induced by the dipole-dipole interaction in Gadolinium Gallium Garnet (GGG) in the high field, low temperature regime.

Developed as part of a final year project on the University of Warwick's undergraduate Mathematics and Physics programme. This project was supervised by Dr Nick d'Ambrumenil and the core portion of the code is based on a Mathematica notebook of his.

This project was conducted as a follow-up to an earlier publication from Nick in PRL: [Phys. Rev. Lett. **114**, 227203 (2015)](https://doi.org/10.1103/PhysRevLett.114.227203).

## Dependencies

The code makes use of a few standard `C` libraries but is otherwise free of external dependencies. Some of the more computationally demanding functions/subroutines are parallelised using `OpenMP`, but as this is handled using `pragma` statements, a serial build is entirely possible.

Copyright (C) C. D. Woodgate 2018-2025. Released under the GNU Lesser General Public License, version 3.
