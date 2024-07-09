# MRI signal model for configuration state imaging
This repository aims to provide a versatile Bloch solver for simulating MR signal evolution for samples with arbitrary NMR spectra and parameters in the presence of field inhomogeneities and spatial variation.

# Design notes
The Bloch solver components are designed to accomodate scalar-, vector-, and matrix-fields. The scalar/vector/matrix part always occupies the first two dimensions of an array, and any subsequent variable (e.g., spatial dimensions, spectroscopic dimension) vary along the remaining dimensions. For example, if we model a 1D object that has varying T1 along the x-axis, the T1 field would have shape (1,1,Nx) (scalar field) and the position field would have shape (3,1,Nx) (vector field), with the x-component (row 1) reflecting the sample points along x.

# Examples
The code used to generate simulated configuration state imaging data in our paper is available in the `sample` directory. Publication forthcoming in MICCAI 2024.
