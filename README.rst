correlcalc
==========

A Python package to calculate 2-point correlation function(2pCF) from
galaxy redshift surveys for any generic model of Cosmology or geometry.

Summary
-------

correlcalc calculates two-point correlation function (2pCF) of
galaxies/quasars using redshift surveys. It can be used for any assumed
geometry or Cosmology model. Using BallTree algorithms to reduce the
computational effort for large datasets, it is faster than brute-force
methods. It takes redshift (z), Right Ascension (RA) and Declination
(DEC) data of galaxies and random catalogs given by redshift survey as
inputs. If random catalog is not provided, it generates one of desired
size based on the input redshift distribution and a mangle polygon file
in .ply format describing the survey geometry. It also calculates
anisotropic 2pCF. Optionally it makes healpix maps of the survey
providing visualization.

Installation
------------

To install this package type "``pip install correlcalc``" in your
terminal. If this method doesn't work

To install the package Download this git repositry and in terminal enter
the folder that contains setup.py and type "``pip install .``" or
"``python setup.py install``"

If you do not have root permission, you can install by adding
"``--user``" at the end of above commands

If you have an older version installed already you can upgrade by
"``pip install correlcalc --upgrade``" command

A note on Dependencies:
~~~~~~~~~~~~~~~~~~~~~~~

All the required dependencies such as sklearn, cython, scipy, numpy etc.
should get automatically installed if installed through pip. In case, if
some of the dependencies do not automatically get installed. The list of
dependencies can be seen in the setup.py file to manually install them.
In case of any problems feel free to raise an issue. "healpix\_util"
package from http://github.com/esheldon/healpix\_util is not available
on pip. So it needs to be manually installed following the commands to
install from git repositry in the above section

Theory
------

The algorithm and formulae used are presented in the paper entitled *A
\`Generic' Recipe for Quick Computation of Two-point Correlation
function*

It is available on arXiv:1710.01723 at https://arxiv.org/abs/1710.01723.

Please cite the same if you use this package or the 'recipe' presented
herein

Usage
-----

Calculation of 2pCF
~~~~~~~~~~~~~~~~~~~

Usage of the package is given in jupyter notebook "Using correlcalc
example.nb" and in ``main.py``

All the methods in correlcalc can be imported using the following
command

``from correlcalc import *``

We first need to define bins (in :math:`c/H_0` units) to calculate 2pCF.
For e.g. to calculate correlation between 0-180Mpc in steps of 6Mpc, we
say

``bins=np.arange(0.002,0.06,0.002)``

To calculate 2pCF using input data file (both ascii and fits files are
supported), use ``tpcf`` method as follows

``correl, poserr=tpcf('/path/to/datfile.dat',bins, randfile='/path/to/randomfile.dat', weights='eq')``

If random file is not available or not provided, we can generate random
catalog by providing the mangle mask file in ``.ply`` format along with
specifying the size of the catalog in multiples of size of data catalog
(default 2x size). To do this

``correl, poserr=tpcf('/path/to/datfile.dat', bins, maskfile='/path/to/maskfile.ply', weights=True, randfact=3)``

This returns ``correl`` and ``poserr`` ``numpy`` arrays corresponding to
Two-point correlation and Poisson error

Keyword Arguments
~~~~~~~~~~~~~~~~~

The following keyword arguments can be included as needed

Data file (Mandatory)
^^^^^^^^^^^^^^^^^^^^^

Data file of galaxy/quasar redshift survey must be passed as the first
argument to both ``tpcf`` and ``atpcf`` methods.

**Supported filetypes**: ascii text files with columns, csv files or
fits files are all supported. Most files provided by SDSS Value added
catalogs should be directly usable.

**To contain**: Any type of file provided must at least have columns
named **Z** (redshift), **RA** (Right Ascension), **DEC** (Declination).
These column names can be in any case.

If one intends to use ``weights=True`` option (must to obtain accurate
results) the data file must also contain radial weights with column
title **radial\_weight** or **WEIGHT\_SYSTOT**

bins (Mandatory)
^^^^^^^^^^^^^^^^

A numpy array with ascending values in :math:`c/H_0` units must be
provided as the second argument to both ``tpcf`` and ``atpcf`` methods.
In case of ``atpcf`` it automatically creates 2D bins as
``bins2d=(bins,bins)`` from provided 1D ``bins``

``randfile=`` Path to random file (semi-Optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If not provided, ``maskfile=`` argument must be given ``.ply`` file.

**Supported filetypes**: ascii text files with columns, csv files or
fits files are all supported. Most files provided by SDSS Value added
catalogs should be directly usable.

**To contain**: Any type of file provided must at least have columns
named **Z** (redshift), **RA** (Right Ascension), **DEC** (Declination).
These column names can be in any case.

If one intends to use ``weights=True`` option (must to obtain accurate
results) the data file must also contain radial weights with column
title **radial\_weight** or **WEIGHT\_SYSTOT**

**Beta Testing:** Beta support for other column titles for weights is
added.

Also added is calculation of weights from n(z) during random catalog
generation.

``mask=`` Path to mangle polygon file (semi-Optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If not provided, ``randfile=`` argument must be provided.

**Supported filetypes**: ``.ply`` file containing Mangle polygons
describing survey geometry in the standard format. Most files provided
by SDSS Value added catalogs should be directly usable.

``randfact=`` (Optional)
^^^^^^^^^^^^^^^^^^^^^^^^

Size of the random catalog in integer multiples of size of data catalog
if random catalog file is not provided. Default value is ``2``

``weights=`` (Optional)
^^^^^^^^^^^^^^^^^^^^^^^

It is highly recommended to use weights argument by providing
``weights=True`` or ``weights='eq'`` to obtain accurate two-point
correlation calculations. This picks up radial weights in the prescribed
format (with column title **radial\_weight** or **WEIGHT\_SYSTOT** )
from the data and random files provided.

``weights=``\ eq'\ ``sets equal weights and hence adds *+1* - This implementation is parallelized and is faster than``\ weights=False\`
implementation on most machines

If ``weights=False``, by default *+1* will be added for each
galaxy/random pair found within the bin instead of adding total weight.
For more details on weights and references, see
http://www.sdss3.org/dr9/tutorials/lss\_galaxy.php

``geometry='flat'`` (Optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Available options**:

``'flat'``\ (default) - for flat geometry of the Universe

``'open'`` - for Open Universe models like Milne

``'close'`` - for Closed Universe

**Customization**

Formulae for calculation of distances between two points (Z1, RA1, DEC1)
and (Z2, RA2, DEC2) is taken from *T. Matsubara, Correlation function in
deep redshift space as a cosmological probe, The Astrophysical Journal
615 (2) (2004) 573*. Using the formulae in this paper, distances squares
(to reduce additional computational time distance squares are calculated
to avoid using expensive ``sqrt`` function every time) are computed in
the ``metrics.pyx`` file for all the above mentioned geometries.
``Cython`` is chosen for implementation to obtain faster results in
building ``BallTree``\ s calculating ``cdist`` and to reduce ``query``
time.

One can customize metric definitions as per one's need by editing this
file. Also **K** (curvature parameter) in the formulae given in this
reference need to be manually changed in the ``metrics.pyx`` for closed
and open cases as per the model. After changing this compile it using
``python metricsetup.py build_ext --inplace``

``cosmology='lcdm'`` (Optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Used to calculate co-moving distances from redshifts.

**Available options**:

``'lcdm'`` (default)- for Lambda CDM model

``'lc'`` - for :math:`R_h=ct` and linear coasting models

**To add**: ``wcdm`` and other popular cosmology models soon

``estimator=`` (Optional)
^^^^^^^^^^^^^^^^^^^^^^^^^

**Available options**:

``'dp'`` - Davis - Peebles estimator (default - fastest)

``'ls'``- Landy - Szalay estimator

``'ph'`` - Peebles- Hauser estimator

``'hew'`` - Hewitt estimator

``'h'`` - Hamilton estimator

For more details on estimator formulae see
https://arxiv.org/pdf/1211.6211.pdf

Calculation of Anisotropic (3D) 2pCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usage of the package is given in jupyter notebook "Using correlcalc
example-anisotropic.nb" and in ``main.py``

All the methods in correlcalc can be imported using the following
command

``from correlcalc import *``

We first need to define bins (in :math:`c/H_0` units) to calculate 2pCF.
For e.g. to calculate correlation between 0-180Mpc in steps of 6Mpc, we
say

``bins=np.arange(0.002,0.06,0.002)``

To calculate anisotropic 2pCF using input data file (both ascii and fits
files are supported), use ``atpcf`` method as follows

``correl3d, poserr=atpcf('/path/to/datfile.dat',binspar, binsper, randfile='/path/to/randomfile.dat', vtype='sigpi', weights=True)``

If random file is not available or not provided, we can generate random
catalog by providing the mangle mask file in ``.ply`` format along with
specifying the size of the catalog in multiples of size of data catalog
(default 2x size). To do this

``correl3d, poserr=atpcf('/path/to/datfile.dat', binspar, binsper, maskfile='/path/to/maskfile.ply', vtype='smu', weights='eq', randfact=3)``

This returns ``correl3d`` and ``poserr`` ``numpy`` arrays corresponding
to anisotropic Two-point correlation and Poisson error

Keyword Arguments
~~~~~~~~~~~~~~~~~

The following keyword arguments can be included as needed

Data file (Mandatory)
^^^^^^^^^^^^^^^^^^^^^

Data file of galaxy/quasar redshift survey must be passed as the first
argument to both ``tpcf`` and ``atpcf`` methods.

**Supported filetypes**: ascii text files with columns, csv files or
fits files are all supported. Most files provided by SDSS Value added
catalogs should be directly usable.

**To contain**: Any type of file provided must at least have columns
named **Z** (redshift), **RA** (Right Ascension), **DEC** (Declination).
These column names can be in any case.

If one intends to use ``weights=True`` option (must to obtain accurate
results) the data file must also contain radial weights with column
title **radial\_weight** or **WEIGHT\_SYSTOT**

binspar (Mandatory)
^^^^^^^^^^^^^^^^^^^

A numpy array with ascending values in :math:`c/H_0` units (for
distances) or :math:`\delta z` as per choice of ``'vtype'`` must be
provided as the second argument to ``atpcf`` method.

binsper (Mandatory)
^^^^^^^^^^^^^^^^^^^

A numpy array with ascending values in :math:`c/H_0` units (for
distances), :math:`z\delta \theta` or :math:`\mu = \cos \alpha` must be
provided as the third argument to ``atpcf`` method.

``randfile=`` Path to random file (semi-Optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If not provided, ``maskfile=`` argument must be given ``.ply`` file.

**Supported filetypes**: ascii text files with columns, csv files or
fits files are all supported. Most files provided by SDSS Value added
catalogs should be directly usable.

**To contain**: Any type of file provided must at least have columns
named **Z** (redshift), **RA** (Right Ascension), **DEC** (Declination).
These column names can be in any case.

If one intends to use ``weights=True`` option the data file must also
contain radial weights with column title **radial\_weight** or
**WEIGHT\_SYSTOT**

**Beta Testing:** Beta support for other column titles for weights is
added.

Also added is calculation of weights from n(z) during random catalog
generation.

``mask=`` Path to mangle polygon file (semi-Optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If not provided, ``randfile=`` argument must be provided.

**Supported filetypes**: ``.ply`` file containing Mangle polygons
describing survey geometry in the standard format. Most files provided
by SDSS Value added catalogs should be directly usable.

``randfact=`` (Optional)
^^^^^^^^^^^^^^^^^^^^^^^^

Size of the random catalog in integer multiples of size of data catalog
if random catalog file is not provided. Default value is ``2``

``weights=`` (Optional)
^^^^^^^^^^^^^^^^^^^^^^^

It is highly recommended to use weights argument by providing
``weights=True`` or ``weights='eq'`` to obtain accurate two-point
correlation calculations. This picks up radial weights in the prescribed
format (with column title **radial\_weight** or **WEIGHT\_SYSTOT** )
from the data and random files provided.

``weights=``\ eq'\ ``sets equal weights and hence adds *+1* - This implementation is parallelized and is faster than``\ weights=False\`
implementation on most machines

If ``weights=False``, by default *+1* will be added for each
galaxy/random pair found within the bin instead of adding total weight.
For more details on weights and references, see
http://www.sdss3.org/dr9/tutorials/lss\_galaxy.php

Metrics in parallel and perpendicular directions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculates anisotropic 2pCF for the following cases.

``vtype=``
^^^^^^^^^^

Valuation method

**Available options**:

``'smu'`` (default)- Calculates 2pCF in s - mu

``'sigpi'`` - Calculates 2pCF using parallel and perpendicular distances

``'ap'`` calculates 2pCF for small :math:`\Delta \theta` and
:math:`z \Delta\theta` . But results can be converted to any cosmology
model of choice (ref: https://arxiv.org/pdf/1312.0003.pdf)

**Customization**

Formulae for calculation of distances in parallel and perpendicular
directions is taken from https://arxiv.org/pdf/1312.0003.pdf. Using the
formulae in this paper, :math:`\Delta z` and :math:`z \Delta \theta` are
computed in the ``metrics.pyx`` file for the above mentioned. ``Cython``
is chosen for implementation to obtain faster results in building
``BallTree``\ s calculating ``cdist`` and to reduce ``query`` time.

One can customize metric definitions as per one's need by editing the
``metrics.pyx`` file. After changing this compile it using
``python metricsetup.py build_ext --inplace``

**To add:**

Direct calculation of distances in LOS and perpendicular to the LOS to
be added to support standard model Cosmology and other popular models.
For now, one needs to manually convert the angular bins to physical
distances to get the approximate results

``cosmology='lcdm'`` (Optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Used to calculate co-moving distances from redshifts.

**Available options**:

``'lcdm'`` (default)- for Lambda CDM model

``'lc'`` - for :math:`R_h=ct` and linear coasting models

**To add**: ``wcdm`` and other popular cosmology models soon

``geometry='flat'`` (Optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Used to calculate co-moving distances between a pair of objects

**Available options**:

``'flat'`` (default)- for Lambda CDM model

``'open'``

``'close'``

``estimator=`` (Optional)
^^^^^^^^^^^^^^^^^^^^^^^^^

**Available options**:

``'dp'`` - Davis - Peebles estimator (default - fastest)

``'ls'``- Landy - Szalay estimator

``'ph'`` - Peebles- Hauser estimator

``'hew'`` - Hewitt estimator

``'h'`` - Hamilton estimator

For more details on estimator formulae see
https://arxiv.org/pdf/1211.6211.pdf
