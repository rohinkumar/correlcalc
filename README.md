# correlcalc
A Python package to calculate 2-point correlation function(2pCF) from galaxy redshift surveys for any generic model of Cosmology or geometry.

correlcalc calculates two-point correlation function (2pCF) of galaxies/quasars using redshift surveys. It can be used for any assumed geometry or Cosmology model. Using BallTree algorithms to reduce the computational effort for large datasets, it is faster than brute-force methods. It takes redshift (z), Right Ascension (RA) and Declination (DEC) data of galaxies and random catalogs given by redshift survey as inputs. If random catalog is not provided, it generates one of desired size based on the input redshift distribution and a mangle polygon file in .ply format describing the survey geometry. It also calculates anisotropic 2pCF. Optionally it makes healpix maps of the survey providing visualization.

The algorithm and formulae used are presented in the paper entitled **A `Generic' Recipe for Quick Computation of Two-point Correlation function**

It is available on arXiv:1710.01723 at  https://arxiv.org/abs/1710.01723.

Please cite the same if you use this package or the 'recipe' presented herein

Usage of the package is given in jupyter notebook "Using correlcalc example.nb" and in main.py 

To install the package Download this git repositry and in terminal enter the folder that contains setup.py and type "pip install ."
