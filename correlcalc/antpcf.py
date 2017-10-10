__author__ = 'Rohin Kumar Y'


# Calculate anisotropic 2pCF
from tpcf import *
# antpcf(dat,datR,bins,parmetric,permetric) returns numpy 2d array DD, RR, DR correl
# poserr(xi,DD) returns (1.0+xi)/np.sqrt(DD)


def atpcf(datfile, bins, **kwargs):
    """Main function to calculate anisotropic 2pCF. Takes multiple arguments such as randfile, maskfile, calculation estimator etc. for different geometry, cosmology models
    Usage of the package is given in jupyter notebook "Using correlcalc example-anisotropic.nb" and in `main.py`

    All the methods in correlcalc can be imported using the following command

    `from correlcalc import *`

    We first need to define bins (in $c/H_0$ units) to calculate 2pCF. For e.g. to calculate correlation between 0-180Mpc in steps of 6Mpc, we say

    `bins=np.arange(0.002,0.06,0.002)`

    To calculate anisotropic 2pCF using input data file (both ascii and fits files are supported), use `atpcf` method as follows

    `correl2d, poserr=atpcf('/path/to/datfile.dat',bins, randfile='/path/to/randomfile.dat', permetric='apzdth', parmetric='apdz', weights=True)`


    If random file is not available or not provided, we can generate random catalog by providing the mangle mask file in `.ply` format along with specifying the size of the catalog in multiples of size of data catalog (default 2x size). To do this

    `correl2d, poserr=atpcf('/path/to/datfile.dat', bins, maskfile='/path/to/maskfile.ply', permetric='apzdth', parmetric='apdz', weights=True, randfact=3)`

    This returns `correl2d` and `poserr` `numpy` arrays corresponding to anisotropic Two-point correlation and Poisson error

    ### Keyword Arguments
    The following keyword arguments can be included as needed

    #### Data file (Mandatory)

    Data file of galaxy/quasar redshift survey must be passed as the first argument to both `tpcf` and `atpcf` methods.

    **Supported filetypes**: ascii text files with columns, csv files or fits files are all supported. Most files provided by SDSS Value added catalogs should be directly usable.

     **To contain**: Any type of file provided must at least have columns named **Z** (redshift), **RA** (Right Ascension), **DEC** (Declination). These column names can be in any case.

     If one intends to use `weights=True` option (must to obtain accurate results) the data file must also contain radial weights with column title **radial_weight** or **WEIGHT_SYSTOT**

    #### bins (Mandatory)

    A numpy array with ascending values in $c/H_0$ units must be provided as the second argument to both `tpcf` and `atpcf` methods. In case of `atpcf` it automatically creates 2D bins as `bins2d=(bins,bins)` from provided 1D `bins`

    #### `randfile=` Path to random file (semi-Optional)

    If not provided, `maskfile=` argument must be given `.ply` file.

    **Supported filetypes**: ascii text files with columns, csv files or fits files are all supported. Most files provided by SDSS Value added catalogs should be directly usable.

     **To contain**: Any type of file provided must at least have columns named **Z** (redshift), **RA** (Right Ascension), **DEC** (Declination). These column names can be in any case.

     If one intends to use `weights=True` option (must to obtain accurate results) the data file must also contain radial weights with column title **radial_weight** or **WEIGHT_SYSTOT**

     **To add:** In future support for other column titles for weights will be added. To also add calculation of weights from n(z) and for random catalog generation.


    #### `mask=` Path to mangle polygon file (semi-Optional)

    If not provided, `randfile=` argument must be provided.

    **Supported filetypes**: `.ply` file containing Mangle polygons describing survey geometry in the standard format. Most files provided by SDSS Value added catalogs should be directly usable.

    #### `randfact=` (Optional)

    Size of the random catalog in integer multiples of size of data catalog if random catalog file is not provided. Default value is `2`

    #### `weights=` (Optional)

    It is highly recommended to use weights argument by providing `weights=True` to obtain accurate two-point correlation calculations. This picks up radial weights in the prescribed format (with column title **radial_weight** or **WEIGHT_SYSTOT** ) from the data and random files provided.

    If `weights=False`, by default *+1* will be added for each galaxy/random pair found within the bin instead of adding total weight. For more details on weights and references, see http://www.sdss3.org/dr9/tutorials/lss_galaxy.php

    #### Metrics in parallel and perpendicular directions

    Currently calculates anisotropic 2pCF only in angular co-ordinates. Currently only small $\Delta \theta$ and $z \Delta\theta$ are supported. But results can be converted to any cosmology model of choice (ref: https://arxiv.org/pdf/1312.0003.pdf)

    #### `parmetric=` (Mandatory)

    Metric to calculate distance in the line-of-sight direction. Currently only `'apdz'` is supported

    #### `permetric=` (Mandatory)

    Metric to calculate distance perpendicular to the line-of-sight direction. Currently only `'apzdth'` is supported


     **Customization**

     Formulae for calculation of distances in parallel and perpendicular directions is taken from https://arxiv.org/pdf/1312.0003.pdf. Using the formulae in this paper, $\Delta z$ and $z \Delta \theta$ are computed in the `metrics.pyx` file for the above mentioned. `Cython` is chosen for implementation to obtain faster results in building `BallTree`s  calculating `cdist` and to reduce `query` time.

     One can customize metric definitions as per one's need by editing the `metrics.pyx` file. After changing this compile it using `python metricsetup.py build_ext --inplace`

     **To add:**

     Direct calculation of distances in LOS and perpendicular to the LOS to be added to support standard model Cosmology and other popular models. For now, one needs to manually convert the angular bins to  physical distances to get the approximate results


    #### `cosmology='lcdm'` (Optional)

    Used to calculate co-moving distances from redshifts.

    **Available options**:

    `'lcdm'` (default)- for Lambda CDM model

     `'lc'` - for $R_h=ct$ and linear coasting models

     **To add**: `wcdm` and other popular cosmology models soon

    #### `estimator=` (Optional)

    **Available options**:

    `'dp'` - Davis - Peebles estimator (default - fastest)

    `'ls'`- Landy - Szalay estimator

    `'ph'` - Peebles- Hauser estimator

    `'hew'` - Hewitt estimator

    `'h'` - Hamilton estimator

    For more details on estimator formulae see https://arxiv.org/pdf/1211.6211.pdf

    """
    # Default function arguments
    rng = np.array([[min(bins), max(bins)], [min(bins), max(bins)]])
    weightsflag = False
    cosmology = 'lcdm'
    # geometry='flat'
    metric = flatdistsq
    randcatfact = 2
    estimator = 'dp'
    binsq = bins**2
    randfile = None
    maskfile = None

    # Options for correl calculation estimators and cosmology models
    mlist = ['dp', 'ls', 'ph', 'hew', 'h']
    clist = ['lcdm', 'lc'] # to add wcdm

    if kwargs is not None:
        for key, value in kwargs.iteritems():
            # print (key, value)
            if key.lower() == 'randfile':
                randfile = value

            elif key.lower() == 'randfact':
                randcatfact = value

            elif key.lower() == 'parmetric':
                if value.lower() == 'apdz':
                    parmetric = APdz
                    # metric='APdz'
                # elif value.lower()=='spar':
                #     #metric='spar'
                #     parmetric=spar
                # elif value.lower()=='s':
                #     #geometry='close'
                #     parmetric=s
                else:
                    print("Incorrect parallel metric argument provided!")
            elif key.lower() == 'permetric':
                if value.lower() == 'apzdth':
                    permetric = APzdth
                    # geometry='flat'
                # elif value.lower()=='sper':
                #     # geometry='open'
                #     permetric=sper
                # elif value.lower()=='mu':
                #     # geometry='close'
                #     permetric=mu
                else:
                    print("Incorrect perpendicular metric provided!")

            elif key.lower() == 'estimator':
                if value.lower() in mlist:
                    estimator = value.lower()
                else:
                    print("Incorrect estimator provided! Using 'dp' as default")
            elif key.lower() == 'cosmology':
                if value.lower() in clist:
                    cosmology = value.lower()
                else:
                    print("Incorrect Cosmology provided! Using 'lcdm' as default")
            elif key.lower() == 'mask':
                maskfile = value
            elif key.lower() == 'weights':
                if value is True:
                    weightsflag = True
                else:
                    weightsflag = False
            else:
                print ("key argument `%s` not valid" % key)
    else:
        print ("Refer documentation to enter valid keyword arguments")

    print("Calculating Anisotropic Correlation function with the following parameters")
    print ("data file=")
    print(datfile)
    print("random file=")
    print(randfile)
    print("Random catalog size factor(if random file is None)=")
    print(randcatfact)
    print("mask/window file=")
    print(maskfile)
    print ("Cosmology=")
    print(cosmology)
    print("Weights=")
    print(weightsflag)
    print("perpendicular metric=")
    print(permetric)
    print("parallel metric=")
    print(parmetric)
    print("Correl estimator=")
    print(estimator)
    print("-------------------------------")
    # Prepare dat from data file
    dat, weights = datprepz(datfile, 'data', cosmology)
    Nd = len(dat)
    # print (weights)
    # Prepare datR from random file or generate a random catalog
    if randfile is None:
        randcatsize = randcatfact*Nd
        if maskfile is None:
            print ("Mask file compulsory. Please provide mask='maskfilepath.ply'")
        else:
            datR = randcatprepz(datfile, randcatsize, maskfile, cosmology)
            rweights = np.array([])
            # randfile='./randcat.dat'
            # datR, rweights=datprep(randfile,'random',cosmology)
    else:
        datR, rweights = datprepz(randfile, 'random', cosmology)

    print ("Calculating anisotropic 2pCF...")

    # Reference: arXiv: 1211.6211
    if estimator == 'dp':
        if weightsflag is False or len(weights) != Nd:
            # print (weightsflag)
            # print(len(weights))
            # print(len(datR))
            DD = aDDcalc(dat, binsq, parmetric, permetric, rng)
            DR = aDRcalc(dat, datR, binsq, parmetric, permetric, rng)
        else:
            # if len(rweights)!=len(datR):
            DD = aDDwcalc(dat, binsq, parmetric, permetric, rng, weights)
            DR = aRDwcalc(dat, datR, binsq, parmetric, permetric, rng, weights)
            # else:
            #     DD=aDDwcalc(dat,binsq,parmetric,permetric,rng,weights)
            #     DR=aDRwcalc(dat,datR,binsq,parmetric,permetric,rng,weights,rweights)

        print ("Using Davis-Peebles estimator")
        correl = (DD/DR)-1.0

    elif estimator == 'ph':
        if weightsflag is False or len(weights) != Nd or len(rweights) != len(datR):
            DD = aDDcalc(dat, binsq, parmetric, permetric, rng)
            RR = aRRcalc(datR, binsq, parmetric, permetric, rng)
        else:
            DD = aDDwcalc(dat, binsq, parmetric, permetric, rng, weights)
            RR = aRRwcalc(datR, binsq, parmetric, permetric, rng, rweights)
        print ("Using Peebles-Hauser estimator")
        correl = (DD/RR)-1.0
    else:
        if weightsflag is False or len(weights) != Nd or len(rweights) != len(datR):
            DD = aDDcalc(dat, binsq, parmetric, permetric, rng)
            RR = aDDcalc(datR, binsq, parmetric, permetric, rng)
            DR = aDRcalc(dat, datR, binsq, parmetric, permetric, rng)
        else:
            DD = aDDwcalc(dat, binsq, parmetric, permetric, rng, weights)
            RR = aRRwcalc(datR, binsq, parmetric, permetric, rng, rweights)
            DR = aRDwcalc(dat, datR, binsq, parmetric, permetric, rng, weights)
        if estimator == 'ls':
            print ("Using Landy-Szalay estimator")
            correl = (DD-2.0*DR+RR)/RR
        elif estimator == 'hew':
            print ("Using Hewett estimator")
            correl = (DD-DR)/RR
        elif estimator == 'h':
            print ("Using Hamilton estimator")
            correl = (DD*RR)/DR**2 - 1.0
    correlerr = poserr(correl, DD*Nd*(Nd-1.0)*0.5)
    print("Anisotropic Two-point correlation=")
    print (correl, correlerr)
    return correl, correlerr


def aDDcalc(dat, bins, parmetric, permetric, rng):
    print ("Calculating anisotropic DD...\n DD=")
    dd = np.zeros((len(bins)-1, len(bins)-1))
    ddbt = BallTree(dat, metric='pyfunc', func=permetric)
    for i in tqdm(xrange(len(dat))):
        ind = ddbt.query_radius(dat[i].reshape(1, -1), max(bins))
        for j in ind:
            dist0 = dist.cdist([dat[i], ], dat[j], parmetric)[0]
            dist1 = dist.cdist([dat[i], ], dat[j], permetric)[0]
            dd += np.histogram2d(dist0, dist1, range=rng, bins=(bins, bins))[0]
    dd[dd == 0] = 1.0
    Nd = len(dat)
    DD = 2.0*dd/(Nd*(Nd-1.0))
    print (DD)
    return DD


def aRRcalc(datR, bins, parmetric, permetric, rng):
    print ("Calculating anisotropic RR...\n RR=")
    rr = np.zeros((len(bins)-1, len(bins)-1))
    rrbt = BallTree(datR, metric='pyfunc', func=permetric)
    for i in tqdm(xrange(len(datR))):
        ind = rrbt.query_radius(datR[i].reshape(1, -1), max(bins))
        for j in ind:
            dist0 = dist.cdist([datR[i], ], datR[j], parmetric)[0]
            dist1 = dist.cdist([datR[i], ], datR[j], permetric)[0]
            rr += np.histogram2d(dist0, dist1, range=rng, bins=(bins, bins))[0]
    rr[rr == 0] = 1.0
    Nr = len(datR)
    RR = 2.0*rr/(Nr*(Nr-1.0))
    print (RR)
    return RR


def aDRcalc(dat, datR, bins, parmetric, permetric, rng):
    print ("Calculating anisotropic DR...\n DR=")
    dr = np.zeros((len(bins)-1, len(bins)-1))
    rrbt = BallTree(datR, metric='pyfunc', func=permetric)
    for i in tqdm(xrange(len(dat))):
        ind = rrbt.query_radius(dat[i].reshape(1, -1), max(bins))
        for j in ind:
            dist0 = dist.cdist([dat[i], ], datR[j], parmetric)[0]
            dist1 = dist.cdist([dat[i], ], datR[j], permetric)[0]
            dr += np.histogram2d(dist0, dist1, range=rng, bins=(bins, bins))[0]
    dr[dr == 0] = 1.0
    Nd = len(dat)
    Nr = len(datR)
    DR = dr/(Nd*Nr)
    print (DR)
    return DR


def aDDwcalc(dat, bins, parmetric, permetric, rng, weights):
    print ("Calculating anisotropic DD with weights...\n DD=")
    dd = np.zeros((len(bins)-1, len(bins)-1))
    ddbt = BallTree(dat, metric='pyfunc', func=permetric)
    for i in tqdm(xrange(len(dat))):
        ind = ddbt.query_radius(dat[i].reshape(1, -1), max(bins))
        for j in ind:
            dist0 = dist.cdist([dat[i], ], dat[j], parmetric)[0]
            dist1 = dist.cdist([dat[i], ], dat[j], permetric)[0]
            dd += np.histogram2d(dist0, dist1, range=rng, bins=(bins, bins), weights=weights[j])[0]
    dd[dd == 0] = 1.0
    Nd = len(dat)
    DD = 2.0*dd/(Nd*(Nd-1.0))
    print (DD)
    return DD


def aRRwcalc(datR, bins, parmetric, permetric, rng, rweights):
    print ("Calculating anisotropic RR with weights...\n RR=")
    rr = np.zeros((len(bins)-1, len(bins)-1))
    rrbt = BallTree(datR, metric='pyfunc', func=permetric)
    for i in tqdm(xrange(len(datR))):
        ind = rrbt.query_radius(datR[i].reshape(1, -1), max(bins))
        for j in ind:
            dist0 = dist.cdist([datR[i], ], datR[j], parmetric)[0]
            dist1 = dist.cdist([datR[i], ], datR[j], permetric)[0]
            rr += np.histogram2d(dist0, dist1, range=rng, bins=(bins, bins), weights=rweights[j])[0]
    rr[rr == 0] = 1.0
    Nr = len(datR)
    RR = 2.0*rr/(Nr*(Nr-1.0))
    print (RR)
    return RR


def aDRwcalc(dat, datR, bins, parmetric, permetric, rng, weights, rweights):
    print ("Calculating anisotropic DR with weights...\n DR=")
    dr = np.zeros((len(bins)-1, len(bins)-1))
    rrbt = BallTree(datR, metric='pyfunc', func=permetric)
    for i in tqdm(xrange(len(dat))):
        ind = rrbt.query_radius(dat[i].reshape(1, -1), max(bins))
        for j in ind:
            dist0 = dist.cdist([dat[i], ], datR[j], parmetric)[0]
            dist1 = dist.cdist([dat[i], ], datR[j], permetric)[0]
            dr += np.histogram2d(dist0, dist1, range=rng, bins=(bins, bins), weights=rweights[j])[0]
    dr[dr == 0] = 1.0
    Nd = len(dat)
    Nr = len(datR)
    DR = dr/(Nd*Nr)
    print (DR)
    return DR


def aRDwcalc(dat, datR, bins, parmetric, permetric, rng, weights):
    print ("Calculating anisotropic RD with weights...\n DR=")
    dr = np.zeros((len(bins)-1, len(bins)-1))
    bt = BallTree(dat, metric='pyfunc', func=permetric)
    for i in tqdm(xrange(len(datR))):
        ind = bt.query_radius(datR[i].reshape(1, -1), max(bins))
        for j in ind:
            dist0 = dist.cdist([datR[i], ], dat[j], parmetric)[0]
            dist1 = dist.cdist([datR[i], ], dat[j], permetric)[0]
            dr += np.histogram2d(dist0, dist1, range=rng, bins=(bins, bins), weights=weights[j])[0]
    dr[dr == 0] = 1.0
    Nd = len(dat)
    Nr = len(datR)
    DR = dr/(Nd*Nr)
    print (DR)
    return DR
