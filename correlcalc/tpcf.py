from __future__ import division
# Samsiddhir Haritoshanam!
__author__ = 'Rohin Kumar Y'

# tpcf(dat, datR=None, randcatsize=2, bins,**kwargs)
#
# **kwargs for choosing geometry - metric 'flat' 'open' 'close'
# **kwargs for choosing xi estimator - 'simple' 'ls' '...'
# import fileios
from tqdm import *
from datprep import *
import numpy as np
from metrics import *
from sklearn.neighbors import BallTree
from scipy.spatial import distance as dist


def tpcf(datfile, bins, **kwargs):
    """Main function to calculate 2pCF. Takes multiple arguments such as randfile, maskfile, calculation estimator etc. for different geometry, cosmology models
    Usage of the package is given in jupyter notebook "Using correlcalc example.nb" and in `main.py`

    All the methods in correlcalc can be imported using the following command

    `from correlcalc import *`

    We first need to define bins (in $c/H_0$ units) to calculate 2pCF. For e.g. to calculate correlation between 0-180Mpc in steps of 6Mpc, we say

    `bins=np.arange(0.002,0.06,0.002)`

    To calculate 2pCF using input data file (both ascii and fits files are supported), use `tpcf` method as follows

    `correl, poserr=tpcf('/path/to/datfile.dat',bins, randfile='/path/to/randomfile.dat', weights=True)`

    If random file is not available or not provided, we can generate random catalog by providing the mangle mask file in `.ply` format along with specifying the size of the catalog in multiples of size of data catalog (default 2x size). To do this

    `correl, poserr=tpcf('/path/to/datfile.dat', bins, maskfile='/path/to/maskfile.ply', weights=True, randfact=3)`

    This returns `correl` and `poserr` `numpy` arrays corresponding to Two-point correlation and Poisson error

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

    #### `geometry='flat'` (Optional)

    **Available options**:

    `'flat'`(default) - for flat geometry of the Universe

     `'open'` - for Open Universe models like Milne

     `'close'` - for Closed Universe

     **Customization**

     Formulae for calculation of distances between two points (Z1, RA1, DEC1) and (Z2, RA2, DEC2) is taken from *T. Matsubara, Correlation function in deep redshift space as a cosmological probe, The Astrophysical Journal 615 (2) (2004) 573*. Using the formulae in this paper, distances squares (to reduce additional computational time distance squares are calculated to avoid using expensive `sqrt` function every time) are computed in the `metrics.pyx` file for all the above mentioned geometries. `Cython` is chosen for implementation to obtain faster results in building `BallTree`s  calculating `cdist` and to reduce `query` time.

     One can customize metric definitions as per one's need by editing this file. Also **K** (curvature parameter) in the formulae given in this reference need to be manually changed in the `metrics.pyx` for closed and open cases as per the model. After changing this compile it using `python metricsetup.py build_ext --inplace`

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
    # weights = np.array([])
    weightsflag = False
    cosmology = 'lcdm'
    geometry = 'flat'
    metric = flatdistsq
    randcatfact = 2
    estimator = 'dp'
    binsq = bins**2
    randfile = None
    maskfile = None

    # Options for correl calculation estimators and cosmology models
    elist = ['dp', 'ls', 'ph', 'hew', 'h']
    clist = ['lcdm', 'lc']

    if kwargs is not None:
        for key, value in kwargs.iteritems():
            # print (key, value)
            if key.lower() == 'randfile':
                randfile = value
            elif key.lower() == 'randfact':
                randcatfact = value
            elif key.lower() == 'geometry':
                if value.lower() == 'flat':
                    geometry = 'flat'
                    metric = flatdistsq
                elif value.lower() == 'open':
                    geometry = 'open'
                    metric = opendistsq
                elif value.lower() == 'close':
                    geometry = 'close'
                    metric = closedistsq
                else:
                    print("Incorrect geometry argument provided! Using flat geometry")
            elif key.lower() == 'estimator':
                if value.lower() in elist:
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
                    # fdat=readinfile(datfile,ftype='internal')
                    # weights=1.0/(1.0+4.0*np.array(fdat['nz']))
                    # weights=weights/np.mean(weights)
                    # print (weights)
                else:
                    weightsflag = False
            else:
                print ("key argument `%s` not valid" % key)
    else:
        print ("Refer documentation to enter valid keyword arguments")

    print("Calculating Correlation function with the following parameters")
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
    print("Geometry=")
    print(geometry)
    print("Correl estimator=")
    print(estimator)
    print ("Weights=")
    print (weightsflag)
    print("-----------------------------------------")
    # Prepare dat from data file
    dat, weights = datprep(datfile, 'data', cosmology)
    Nd = len(dat)
    # print (weights)
    # Prepare datR from random file or generate a random catalog
    if randfile is None:
        randcatsize = randcatfact*Nd
        if maskfile is None:
            print ("Mask file compulsory. Please provide mask='maskfilepath.ply'")
        else:
            datR = randcatprep(datfile, randcatsize, maskfile, cosmology)
            rweights = np.array([])
            # randfile='./randcat.dat'
            # datR, rweights=datprep(randfile,'random',cosmology)
    else:
        datR, rweights = datprep(randfile, 'random', cosmology)

    # if len(weights)!=0:
        # rfdat=readinfile(randfile,ftype='internal')
        # rweights=1.0/(1.0+4.0*np.array(rfdat['nz']))
        # rweights=rweights/np.mean(rweights)
        # print (rweights)
    # Nr=len(datR)
    print ("Calculating 2pCF...")
    # f=(1.0*Nrd)/N
    # print (weights)
    # Reference: arXiv: 1211.6211
    if estimator == 'dp':
        if weightsflag is False or len(weights) != Nd:
            # print (weightsflag)
            # print(len(weights))
            # print(len(datR))
            DD = DDcalc(dat, binsq, metric)
            DR = DRcalc(dat, datR, binsq, metric)
        else:
            # if len(rweights)!=len(datR):
            DD = DDwcalc(dat, binsq, metric, weights)
            DR = RDwcalc(dat, datR, binsq, metric, weights)
            # else:
            #     DD=DDwcalc(dat,binsq,metric,weights)
            #     DR=DRwcalc(dat,datR,binsq,metric,rweights)
        print ("Using Davis-Peebles estimator")
        correl = (DD/DR)-1.0

    elif estimator == 'ph':
        if weightsflag is False or len(weights) != Nd:
            DD = DDcalc(dat, binsq, metric)
            RR = RRcalc(datR, binsq, metric)
        else:
            DD = DDwcalc(dat, binsq, metric, weights)
            RR = RRwcalc(datR, binsq, metric, rweights)
        print ("Using Peebles-Hauser estimator")
        correl = (DD/RR)-1.0
    else:
        if weightsflag == False or len(weights) != Nd:
            DD = DDcalc(dat, binsq, metric)
            RR = RRcalc(datR, binsq, metric)
            DR = DRcalc(dat, datR, binsq, metric)
        else:
            DD = DDwcalc(dat, binsq, metric, weights)
            RR = RRwcalc(datR, binsq, metric, rweights)
            DR = RDwcalc(dat, datR, binsq, metric, weights)
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
    print("Two-point correlation=")
    print (correl, correlerr)
    return correl, correlerr


def DDcalc(dat, bins, metric):
    print ("Calculating DD...\n DD=")
    DD = autocorr(dat, bins, metric)
    DD [DD == 0] = 1.0
    Nd = len(dat)
    DD = 2.0*DD/(Nd*(Nd-1.0))
    print (DD)
    return DD


def RRcalc(datR, bins, metric):
    print ("Calculating RR...\n RR=")
    RR = autocorr(datR, bins, metric)
    RR [RR == 0] = 1.0
    Nr = len(datR)
    RR = 2.0*RR/(Nr*(Nr-1.0))
    print (RR)
    return RR


def DRcalc(dat, datR, bins, metric):
    print ("Calculating DR...\n DR=")
    DR = crosscorr(dat, datR, bins, metric)
    DR[DR == 0] = 1.0
    Nd = len(dat)
    Nr = len(datR)
    DR = DR/(Nd*Nr)
    print (DR)
    return DR


def autocorr(dat, bins, metric):
    bt = BallTree(dat, metric='pyfunc', func=metric)
    counts_DD = bt.two_point_correlation(dat, bins)
    DD = np.diff(counts_DD)
    return DD


def crosscorr(dat, datR, bins, metric):
    rbt = BallTree(datR, metric='pyfunc', func=metric)
    counts_DR = rbt.two_point_correlation(dat, bins)
    DR = np.diff(counts_DR)
    return DR


def poserr(xi, DD):
    print ("Calculating Poisson error")
    return (1.0+xi)/np.sqrt(DD)
# alternatively
# rbt=BallTree(dat,metric='pyfunc',func=metric)
# counts_RD=rbt.two_point_correlation(dat,bins)


def DDwcalc(dat, bins, metric, weights):
    print ("Calculating DD with weights...\n DD=")
    DD = autocorrw(dat, bins, metric, weights)
    DD[DD == 0] = 1.0
    Nd = len(dat)
    DD = 2.0*DD/(Nd*(Nd-1.0))
    print (DD)
    return DD


def RRwcalc(datR, bins, metric, weights):
    print ("Calculating RR with weights...\n RR=")
    RR = autocorrw(datR, bins, metric, weights)
    RR[RR == 0] = 1.0
    Nr = len(datR)
    RR = 2.0*RR/(Nr*(Nr-1.0))
    print (RR)
    return RR


def DRwcalc(dat, datR, bins, metric, rweights):
    print ("Calculating DR with weights...\n DR=")
    DR = crosscorrw(dat, datR, bins, metric, rweights)
    DR[DR == 0] = 1.0
    Nd = len(dat)
    Nr = len(datR)
    DR = DR/(Nd*Nr)
    print (DR)
    return DR


def RDwcalc(dat, datR, bins, metric, weights):
    print ("Calculating RD with weights...\n RD=")
    DR = crosscorrwrd(dat, datR, bins, metric, weights)
    DR[DR == 0] = 1.0
    Nd = len(dat)
    Nr = len(datR)
    DR = DR/(Nd*Nr)
    print (DR)
    return DR


def autocorrw(dat, bins, metric, weights):
    bt = BallTree(dat, metric='pyfunc', func=metric)
    DD = np.zeros(len(bins)-1)
    for i in tqdm(xrange(len(dat))):
        ind = bt.query_radius(dat[i].reshape(1, -1), max(bins))
        # wts=np.array([])
        for j in ind:
            dist0 = dist.cdist([dat[i], ], dat[j], metric)[0]
            DD += np.histogram(dist0, bins=bins, weights=weights[j])[0]
            # print (dist0,weights[j])
    return DD


def crosscorrw(dat, datR, bins, metric, rweights):
    rbt = BallTree(datR, metric='pyfunc', func=metric)
    DR = np.zeros(len(bins)-1)
    for i in tqdm(xrange(len(dat))):
        ind = rbt.query_radius(dat[i].reshape(1, -1), max(bins))
        # wts=np.array([])
        for j in ind:
            dist0 = dist.cdist([dat[i], ], datR[j], metric)[0]
            DR += np.histogram(dist0, bins=bins, weights=rweights[j])[0]
            # print (dist0,weights[j])
    return DR


def crosscorrwrd(dat, datR, bins, metric, weights):
    bt = BallTree(dat, metric='pyfunc', func=metric)
    RD = np.zeros(len(bins)-1)
    for i in tqdm(xrange(len(datR))):
        ind = bt.query_radius(datR[i].reshape(1, -1), max(bins))
        #  wts=np.array([])
        for j in ind:
            dist0 = dist.cdist([datR[i], ], dat[j], metric)[0]
            RD += np.histogram(dist0, bins=bins, weights=weights[j])[0]
            # print (dist0,weights[j])
    return RD
