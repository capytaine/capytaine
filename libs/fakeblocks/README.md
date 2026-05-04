# Fakeblocks

An implementation of Hierarchical Matrices and block Toeplitz matrices in Capytaine.
This repository contains code defining some additional meshes and matrices classes compatible with Capytaine for advanced resolutions with the Boundary Element Method.

Actually this code used to be part of Capytaine before version 3.0, it has been extracted in another repository for several reasons:
- licensing: the present repository contains GPL licensed code
- lack of documentation and testing,
- lack of usefulness, the performance of the naive Python implementation is not great.

Right now, it is mostly kept as a test of the modularity of Capytaine, allowing external custom mesh and engine classes to be used with its solver.
Since it still works with Capytaine v3, it could be resucitated in the future if necessary.
