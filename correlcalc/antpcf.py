__author__ = 'Rohin Kumar Y'


# Calculate anisotropic 2pCF
from tpcf import *
import scipy as sp
# antpcf(dat,datR,bins,parmetric,permetric) returns numpy 2d array DD, RR, DR correl
# poserr(xi,DD) returns (1.0+xi)/np.sqrt(DD)


def atpcf(datfile, binspar, binsper, **kwargs):
    """Main function to calculate anisotropic 2pCF. Takes multiple arguments such as randfile, maskfile, calculation estimator etc. for different geometry, cosmology models
    Usage of the package is given in jupyter notebook "Using correlcalc example-anisotropic.nb" and in `main.py`

    All the methods in correlcalc can be imported using the following command

    `from correlcalc import *`

    We first need to define bins (in $c/H_0$ units) to calculate 2pCF. For e.g. to calculate correlation between 0-180Mpc in steps of 6Mpc, we say

    `bins=np.arange(0.002,0.06,0.002)`

    To calculate anisotropic 2pCF using input data file (both ascii and fits files are supported), use `atpcf` method as follows

    `correl3d, poserr=atpcf('/path/to/datfile.dat',binspar, binsper, randfile='/path/to/randomfile.dat', vtype='sigpi', weights=True)`


    If random file is not available or not provided, we can generate random catalog by providing the mangle mask file in `.ply` format along with specifying the size of the catalog in multiples of size of data catalog (default 2x size). To do this

    `correl3d, poserr=atpcf('/path/to/datfile.dat', binspar, binsper, maskfile='/path/to/maskfile.ply', vtype='smu', weights='eq', randfact=3)`

    This returns `correl3d` and `poserr` `numpy` arrays corresponding to anisotropic Two-point correlation and Poisson error

    ### Keyword Arguments
    The following keyword arguments can be included as needed

    #### Data file (Mandatory)

    Data file of galaxy/quasar redshift survey must be passed as the first argument to both `tpcf` and `atpcf` methods.

    **Supported filetypes**: ascii text files with columns, csv files or fits files are all supported. Most files provided by SDSS Value added catalogs should be directly usable.

     **To contain**: Any type of file provided must at least have columns named **Z** (redshift), **RA** (Right Ascension), **DEC** (Declination). These column names can be in any case.

     If one intends to use `weights=True` option (must to obtain accurate results) the data file must also contain radial weights with column title **radial_weight** or **WEIGHT_SYSTOT**

    #### binspar (Mandatory)

    A numpy array with ascending values in $c/H_0$ units (for distances) or $\delta z$ as per choice of `'vtype'` must be provided as the second argument to  `atpcf` method.

    #### binsper (Mandatory)

    A numpy array with ascending values in $c/H_0$ units (for distances), $z\delta \theta$ or $\mu = \cos \alpha$ must be provided as the third argument to  `atpcf` method.


    #### `randfile=` Path to random file (semi-Optional)

    If not provided, `maskfile=` argument must be given `.ply` file.

    **Supported filetypes**: ascii text files with columns, csv files or fits files are all supported. Most files provided by SDSS Value added catalogs should be directly usable.

     **To contain**: Any type of file provided must at least have columns named **Z** (redshift), **RA** (Right Ascension), **DEC** (Declination). These column names can be in any case.

     If one intends to use `weights=True` option the data file must also contain radial weights with column title **radial_weight** or **WEIGHT_SYSTOT**

    **Beta Testing:** Beta support for other column titles for weights is added.

     Also added is calculation of weights from n(z) during random catalog generation.

    #### `mask=` Path to mangle polygon file (semi-Optional)

    If not provided, `randfile=` argument must be provided.

    **Supported filetypes**: `.ply` file containing Mangle polygons describing survey geometry in the standard format. Most files provided by SDSS Value added catalogs should be directly usable.

    #### `randfact=` (Optional)

    Size of the random catalog in integer multiples of size of data catalog if random catalog file is not provided. Default value is `2`

    #### `weights=` (Optional)

    It is highly recommended to use weights argument by providing `weights=True` or `weights='eq'` to obtain accurate two-point correlation calculations. This picks up radial weights in the prescribed format (with column title **radial_weight** or **WEIGHT_SYSTOT** ) from the data and random files provided.

    `weights=`eq'` sets equal weights and hence adds *+1* - This implementation is parallelized and is faster than `weights=False` implementation on most machines

    If `weights=False`, by default *+1* will be added for each galaxy/random pair found within the bin instead of adding total weight. For more details on weights and references, see http://www.sdss3.org/dr9/tutorials/lss_galaxy.php

    #### Metrics in parallel and perpendicular directions

    Calculates anisotropic 2pCF for the following cases.

    #### `vtype=`

    Valuation method

    **Available options**:

    `'smu'` (default)- Calculates 2pCF in s - mu

     `'sigpi'` - Calculates 2pCF using parallel and perpendicular distances

    `'ap'` calculates 2pCF for small $\Delta \theta$ and $z \Delta\theta$ . But results can be converted to any cosmology model of choice (ref: https://arxiv.org/pdf/1312.0003.pdf)

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

    #### `geometry='flat'` (Optional)

    Used to calculate co-moving distances between a pair of objects

    **Available options**:

    `'flat'` (default)- for Lambda CDM model

     `'open'`

     `'close'`


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
    global binsparv
    global binsperv
    global maxrad
    global dat
    global datR
    global Nd
    global Nr
    DD = DR = RD = RR = np.zeros((len(binspar)-1, len(binsper)-1))
    weightsflag = True
    useones = True
    cosmology = 'lcdm'
    sflag = True
    geometry='flat'
    filtermetric = flatdistsq
    permetric = musqlcdmf
    parmetric = flatdistsq
    vtype = 'smu'
    randcatfact = 2
    estimator = 'dp'
    binsparv = binspar**2
    binsperv = binsper**2
    randfile = None
    maskfile = None

    # Options for correl calculation estimators and cosmology models
    mlist = ['dp', 'ls', 'ph', 'hew', 'h']
    clist = ['lcdm', 'lc'] # to add wcdm
    glist = ['flat', 'open', 'close']
    parper = ['ap', 'sigpi', 'smu']

    if kwargs is not None:
        for key, value in kwargs.items():
            # print (key, value)
            if key.lower() == 'randfile':
                randfile = value

            elif key.lower() == 'randfact':
                randcatfact = value

            elif key.lower() == 'geometry':
                if value.lower() in glist:
                    geometry = value.lower()
                #     geometry = 'flat'
                #     filtermetric = flatdistsq
                # elif value.lower() == 'open':
                #     geometry = 'open'
                #     filtermetric = opendistsq
                # elif value.lower() == 'close':
                #     geometry = 'close'
                #     filtermetric = closedistsq

            elif key.lower() == 'cosmology':
                if value.lower() in clist:
                    cosmology = value.lower()
                else:
                    print("Incorrect Cosmology provided! Using 'lcdm' as default")

            elif key.lower() == 'vtype':
                if value.lower() in parper:
                    vtype = value.lower()

            elif key.lower() == 'estimator':
                if value.lower() in mlist:
                    estimator = value.lower()
                else:
                    print("Incorrect estimator provided! Using 'dp' as default")

            elif key.lower() == 'mask':
                maskfile = value
            elif key.lower() == 'weights':
                if value is True:
                    weightsflag = True
                    useones = False
                elif isinstance(value, str):
                    if value.lower() == 'eq':
                        weightsflag = True
                        useones = True
                else:
                    weightsflag = False
            else:
                print ("key argument `%s` not valid" % key)
    else:
        print ("Refer documentation to enter valid keyword arguments")

    if vtype == 'ap':
        parmetric = APdz
        binsparv = binspar
        binsperv = binsper
        sflag = False
        filtermetric = APzdth
        permetric = APzdth
        maxrad = max(np.sqrt(binsparv**2 + binsperv**2))

    elif vtype == 'smu':
        # binsparv = binspar**2
        # binsperv = binsper**2
        maxrad = max(binsparv)
        if geometry == 'open':
            parmetric = opendistsq
            filtermetric = opendistsq
            if cosmology == 'lc':
                permetric = musqlco
            else:
                permetric = musqlcdmo

        elif geometry == 'close':
            parmetric = closedistsq
            filtermetric = closedistsq
            if cosmology == 'lc':
                permetric = musqlcc
            else:
                permetric = musqlcdmc
        else:
            parmetric = flatdistsq
            filtermetric = flatdistsq
            if cosmology == 'lc':
                permetric = musqlcf
            else:
                permetric = musqlcdmf

    elif vtype == 'sigpi':
        # binsparv = binspar**2
        # binsperv = binsper**2
        maxrad = max(binsparv+binsperv)
        if geometry == 'open':
            filtermetric = opendistsq
            if cosmology == 'lc':
                parmetric = sparsqlc
                permetric = spersqlco
            else:
                parmetric = sparsqlcdm
                permetric = spersqlcdmo

        elif geometry == 'close':
            filtermetric = closedistsq
            if cosmology == 'lc':
                parmetric = sparsqlc
                permetric = spersqlcc
            else:
                parmetric = sparsqlcdm
                permetric = spersqlcdmc
        else:
            filtermetric = flatdistsq
            if cosmology == 'lc':
                parmetric = sparsqlc
                permetric = spersqlcf
            else:
                parmetric = sparsqlcdm
                permetric = spersqlcdmf

    else:
        print ("No valid valuation method provided. Using 'smu' as default")


    print("Calculating Anisotropic Correlation function with the following parameters")
    print("data file=")
    print(datfile)
    print("random file=")
    print(randfile)
    print("Random catalog size factor(if random file is None)=")
    print(randcatfact)
    print("mask/window file=")
    print(maskfile)
    print ("Cosmology=")
    print(cosmology)
    print ("Geometry=")
    print (geometry)
    print("Weights=")
    print(weightsflag)
    print ("Using ones as weights?=")
    print (useones)
    print("perpendicular metric=")
    print(permetric)
    print("parallel metric=")
    print(parmetric)
    print("Correl estimator=")
    print(estimator)
    print("Valuation type=")
    print(vtype)
    print ("binsparv=")
    print (binsparv)
    print ("binsperv=")
    print (binsperv)
    print("---------------------------------------------------------------------------")

    if sflag is False:
        # Prepare dat from data file
        dat, weights = datprepz(datfile, 'data', cosmology)
        Nd = len(dat)
        # Prepare datR from random file or generate a random catalog
        if randfile is None:
            randcatsize = randcatfact*Nd
            if maskfile is None:
                print ("Mask file compulsory. Please provide mask='maskfilepath.ply'")
            else:
                datR, rweights = randcatprepz(datfile, randcatsize, maskfile, cosmology)
        else:
            datR, rweights = datprepz(randfile, 'random', cosmology)

    else:
        # Prepare dat from data file
        dat, weights = datprep(datfile, 'data', cosmology)

        Nd = len(dat)
        # Prepare datR from random file or generate a random catalog
        if randfile is None:
            randcatsize = randcatfact*Nd
            if maskfile is None:
                print ("Mask file compulsory. Please provide mask='maskfilepath.ply'")
            else:
                datR, rweights = randcatprep(datfile, randcatsize, maskfile, cosmology)
        else:
            datR, rweights = datprep(randfile, 'random', cosmology)

    Nr = len(datR)
    fact = (1.0*Nr)/Nd
    global adbt
    global arbt

    print ("Creating BallTree for data points using ...")
    print (filtermetric)
    adbt = BallTree(dat, metric='pyfunc', func=filtermetric)

    print ("Creating BallTree for random points using ...")
    print (filtermetric)
    arbt = BallTree(datR, metric='pyfunc', func=filtermetric)

    rng = np.array([[min(binsparv), max(binsparv)], [min(binsperv), max(binsperv)]])
    print ("Calculating anisotropic 2pCF...")

    # Reference: arXiv: 1211.6211
    if estimator == 'dp':
        if weightsflag is False: # or len(weights) != Nd
            # print (weightsflag)
            # print(len(weights))
            # print(len(datR))
            DD = aDDcalc(dat, binsparv, binsperv, parmetric, permetric, rng)
            DR = aDRcalc(dat, datR, binsparv, binsperv, parmetric, permetric, rng)
            RD = aRDcalc(dat, datR, binsparv, binsperv, parmetric, permetric, rng)
        else:
            # if len(rweights)!=len(datR):
            # DD = aDDwcalc(dat, binsq, parmetric, permetric, rng, weights)
            if useones is True or len(weights) != Nd:
                weights = np.ones(Nd)
                rweights = np.ones(Nr)
            print ("Calculating anisotropic DD with weights (parallelized)...\n DD=")
            DD = amulti_autocp(dat, binsparv, binsperv, parmetric, permetric, rng, weights, Nd, pcpus)
            # DR = aRDwcalc(dat, datR, binsq, parmetric, permetric, rng, weights)
            print ("Calculating anisotropic DR with weights (parallelized)...\n DR=")
            DR = amulti_crosscp(dat, datR, binsparv, binsperv, parmetric, permetric, rng, weights, Nr, pcpus)
            print ("Calculating anisotropic RD with weights (parallelized)...\n RD=")
            RD = amulti_crosscpr(dat, datR, binsparv, binsperv, parmetric, permetric, rng, rweights, Nd, pcpus)
            # else:
            #     DD=aDDwcalc(dat,binsq,parmetric,permetric,rng,weights)
            #     DR=aDRwcalc(dat,datR,binsq,parmetric,permetric,rng,weights,rweights)

        print ("Using Davis-Peebles estimator")
        correl = fact*(DD*2.0/(DR+RD))-1.0

    elif estimator == 'ph':
        if weightsflag is False: # or len(weights) != Nd or len(rweights) != len(datR):
            DD = aDDcalc(dat, binsparv, binsperv, parmetric, permetric, rng)
            RR = aRRcalc(datR, binsparv, binsperv, parmetric, permetric, rng)
        else:
            if useones is True or len(weights) != Nd:
                weights = np.ones(Nd)
                rweights = np.ones(Nr)
            print ("Calculating anisotropic DD with weights (parallelized)...\n DD=")
            # DD = aDDwcalc(dat, binsq, parmetric, permetric, rng, weights)
            DD = amulti_autocp(dat, binsparv, binsperv, parmetric, permetric, rng, weights, Nd, pcpus)
            # if len(rweights) != Nr:
            #     RR = aRRcalc(datR, binsparv, binsperv, parmetric, permetric, rng)
            # else:
            print ("Calculating anisotropic RR with weights (parallelized)...\n RR=")
            RR = amulti_autocpr(datR, binsparv, binsperv, parmetric, permetric, rng, rweights, Nr, pcpus)
        print ("Using Peebles-Hauser estimator")
        correl = fact**2*(DD/RR)-1.0
    else:
        if weightsflag is False: # or len(weights) != Nd or len(rweights) != len(datR):
            DD = aDDcalc(dat, binsparv, binsperv, parmetric, permetric, rng)
            RR = aRRcalc(datR, binsparv, binsperv, parmetric, permetric, rng)
            DR = aDRcalc(dat, datR, binsparv, binsperv, parmetric, permetric, rng)
            RD = aRDcalc(dat, datR, binsparv, binsperv, parmetric, permetric, rng)
        else:
            if useones is True or len(weights) != Nd:
                weights = np.ones(Nd)
                rweights = np.ones(Nr)
            print ("Calculating anisotropic DD with weights (parallelized)...\n DD=")
            # DD = aDDwcalc(dat, binsq, parmetric, permetric, rng, weights)
            DD = amulti_autocp(dat, binsparv, binsperv, parmetric, permetric, rng, weights, Nd, pcpus)
            # print ("Calculating anisotropic RR with weights (parallelized)...\n RR=")
            # RR = aRRwcalc(datR, binsq, parmetric, permetric, rng, rweights)
            # RR = amulti_autocpr(datR, binsq, parmetric, permetric, rng, rweights, Nr, pcpus)
            # DR = aRDwcalc(dat, datR, binsq, parmetric, permetric, rng, weights)
            print ("Calculating anisotropic DR with weights (parallelized)...\n DR=")
            DR = amulti_crosscp(dat, datR, binsparv, binsperv, parmetric, permetric, rng, weights, Nr, pcpus)
            print ("Calculating anisotropic RD with weights (parallelized)...\n RD=")
            RD = amulti_crosscpr(dat, datR, binsparv, binsperv, parmetric, permetric, rng, rweights, Nd, pcpus)
            # if len(rweights) != Nr:
            #     RR = aRRcalc(datR, binsparv, binsperv, parmetric, permetric, rng)
            # else:
            print ("Calculating anisotropic RR with weights (parallelized)...\n RR=")
            RR = amulti_autocpr(datR, binsparv, binsperv, parmetric, permetric, rng, rweights, Nr, pcpus)
        if estimator == 'ls':
            print ("Using Landy-Szalay estimator")
            correl = fact**2*(DD/RR)-fact*(DR+RD)/RR+1.0
            # correl = fact**2*(DD/RR)-2.0*fact*(DR/RR)+1.0
        elif estimator == 'hew':
            print ("Using Hewett estimator")
            correl = fact**2*(DD/RR)-fact*0.5*(DR+RD)/RR
            # correl = fact**2*(DD/RR)-fact*(DR/RR)
        elif estimator == 'h':
            print ("Using Hamilton estimator")
            correl = (4.0*DD*RR)/(DR+RD)**2 - 1.0
            # correl = (DD*RR)/DR**2 - 1.0
    correlerr = poserr(correl, DD)
    print("Anisotropic Two-point correlation=")
    np.savetxt("aDD_"+str(cosmology)+"_"+str(geometry)+"_"+str(vtype)+"_"+str(estimator)+".txt", DD)
    np.savetxt("aDR_"+str(cosmology)+"_"+str(geometry)+"_"+str(vtype)+"_"+str(estimator)+".txt", DR)
    np.savetxt("aRD_"+str(cosmology)+"_"+str(geometry)+"_"+str(vtype)+"_"+str(estimator)+".txt", RD)
    np.savetxt("aRR_"+str(cosmology)+"_"+str(geometry)+"_"+str(vtype)+"_"+str(estimator)+".txt", RR)
    np.savetxt("abinspar_"+str(cosmology)+"_"+str(geometry)+"_"+str(vtype)+"_"+str(estimator)+".txt", binspar)
    np.savetxt("abinsper_"+str(cosmology)+"_"+str(geometry)+"_"+str(vtype)+"_"+str(estimator)+".txt", binsper)
    np.savetxt("atpcf_"+str(cosmology)+"_"+str(geometry)+"_"+str(vtype)+"_"+str(estimator)+".txt", correl)
    np.savetxt("atpcferr_"+str(cosmology)+"_"+str(geometry)+"_"+str(vtype)+"_"+str(estimator)+".txt", correlerr)
    print (correl, correlerr)
    return correl, correlerr


def aDDcalc(dat, binspar, binsper, parmetric, permetric, rng):
    print ("Calculating anisotropic DD...\n DD=")
    dd = np.zeros((len(binspar)-1, len(binsper)-1))
    for i in tqdm(range(len(dat))):
        ind = adbt.query_radius(dat[i].reshape(1, -1), maxrad)
        for j in ind:
            dist0 = dist.cdist([dat[i], ], dat[j[j>i]], parmetric)[0]
            # print("dist0")
            # print dist0
            dist1 = dist.cdist([dat[i], ], dat[j[j>i]], permetric)[0]
            # print("dist1")
            # print dist1
            # print np.histogram2d(dist0, dist1, range=rng, bins=(binspar, binsper))[0]
            dd += np.histogram2d(dist0, dist1, range=rng, bins=(binspar, binsper))[0]
            # print ("rng")
            # print rng
            # print("binspar")
            # print binspar
            # print("binsper")
            # print binsper
            # print dd
    dd[dd == 0] = 1.0
    # Nd = len(dat)
    # DD = dd/(Nd*(Nd-1.0))
    print (dd)
    return dd


def aRRcalc(datR, binspar, binsper, parmetric, permetric, rng):
    print ("Calculating anisotropic RR...\n RR=")
    rr = np.zeros((len(binspar)-1, len(binsper)-1))
    # rrbt = BallTree(datR, metric='pyfunc', func=permetric)
    for i in tqdm(range(len(datR))):
        ind = arbt.query_radius(datR[i].reshape(1, -1), maxrad)
        for j in ind:
            dist0 = dist.cdist([datR[i], ], datR[j[j>i]], parmetric)[0]
            dist1 = dist.cdist([datR[i], ], datR[j[j>i]], permetric)[0]
            rr += np.histogram2d(dist0, dist1, range=rng, bins=(binspar, binsper))[0]
    rr[rr == 0] = 1.0
    # Nr = len(datR)
    # RR = rr/(Nr*(Nr-1.0))
    print (rr)
    return rr


def aDRcalc(dat, datR, binspar, binsper, parmetric, permetric, rng):
    print ("Calculating anisotropic DR...\n DR=")
    dr = np.zeros((len(binspar)-1, len(binsper)-1))
    # rrbt = BallTree(datR, metric='pyfunc', func=permetric)
    for i in tqdm(range(len(dat))):
        ind = arbt.query_radius(dat[i].reshape(1, -1), maxrad)
        for j in ind:
            dist0 = dist.cdist([dat[i], ], datR[j[j>i]], parmetric)[0]
            dist1 = dist.cdist([dat[i], ], datR[j[j>i]], permetric)[0]
            dr += np.histogram2d(dist0, dist1, range=rng, bins=(binspar, binsper))[0]
    dr[dr == 0] = 1.0
    # Nd = len(dat)
    # Nr = len(datR)
    # DR = dr/(Nd*Nr)
    print (dr)
    return dr


def aRDcalc(dat, datR, binspar, binsper, parmetric, permetric, rng):
    print ("Calculating anisotropic RD...\n RD=")
    rd = np.zeros((len(binspar)-1, len(binsper)-1))
    # rrbt = BallTree(datR, metric='pyfunc', func=permetric)
    for i in tqdm(range(len(datR))):
        ind = arbt.query_radius(datR[i].reshape(1, -1), maxrad)
        for j in ind:
            dist0 = dist.cdist([datR[i], ], dat[j[j>i]], parmetric)[0]
            dist1 = dist.cdist([datR[i], ], dat[j[j>i]], permetric)[0]
            rd += np.histogram2d(dist0, dist1, range=rng, bins=(binspar, binsper))[0]
    rd[rd == 0] = 1.0
    # Nd = len(dat)
    # Nr = len(datR)
    # DR = dr/(Nd*Nr)
    print (rd)
    return rd


def aDDwcalc(dat, binspar, binsper, parmetric, permetric, rng, weights):
    print ("Calculating anisotropic DD with weights...\n DD=")
    dd = np.zeros((len(binspar)-1, len(binsper)-1))
    # ddbt = BallTree(dat, metric='pyfunc', func=permetric)
    for i in tqdm(range(len(dat))):
        ind = adbt.query_radius(dat[i].reshape(1, -1), maxrad)
        for j in ind:
            dist0 = dist.cdist([dat[i], ], dat[j[j>i]], parmetric)[0]
            dist1 = dist.cdist([dat[i], ], dat[j[j>i]], permetric)[0]
            dd += np.histogram2d(dist0, dist1, range=rng, bins=(binspar, binsper), weights=weights[j[j>i]])[0]
    dd[dd == 0] = 1.0
    # Nd = len(dat)
    # DD = dd/(Nd*(Nd-1.0)) # factor of 2 cancels with 1/2 that needs to be done to remove double counting of pairs
    # print (dd)
    return dd


def aRRwcalc(datR, binspar, binsper, parmetric, permetric, rng, rweights):
    print ("Calculating anisotropic RR with weights...\n RR=")
    rr = np.zeros((len(binspar)-1, len(binsper)-1))
    for i in tqdm(range(len(datR))):
        ind = arbt.query_radius(datR[i].reshape(1, -1), maxrad)
        for j in ind:
            dist0 = dist.cdist([datR[i], ], datR[j[j>i]], parmetric)[0]
            dist1 = dist.cdist([datR[i], ], datR[j[j>i]], permetric)[0]
            rr += np.histogram2d(dist0, dist1, range=rng, bins=(binspar, binsper), weights=rweights[j[j>i]])[0]
    rr[rr == 0] = 1.0
    # Nr = len(datR)
    # RR = rr/(Nr*(Nr-1.0)) # factor of 2 cancels with 1/2 that needs to be done to remove double counting of pairs
    # print (rr)
    return rr


def aDRwcalc(dat, datR, binspar, binsper, parmetric, permetric, rng, rweights):
    print ("Calculating anisotropic DR with weights...\n DR=")
    dr = np.zeros((len(binspar)-1, len(binsper)-1))
    # rrbt = BallTree(datR, metric='pyfunc', func=permetric)
    for i in tqdm(range(len(dat))):
        ind = arbt.query_radius(dat[i].reshape(1, -1), maxrad)
        for j in ind:
            dist0 = dist.cdist([dat[i], ], datR[j], parmetric)[0]
            dist1 = dist.cdist([dat[i], ], datR[j], permetric)[0]
            dr += np.histogram2d(dist0, dist1, range=rng, bins=(binspar, binsper), weights=rweights[j])[0]
    dr[dr == 0] = 1.0
    # Nd = len(dat)
    # Nr = len(datR)
    # DR = dr/(Nd*Nr)
    # print (dr/2.0)
    return dr/2.0


def aRDwcalc(dat, datR, binspar, binsper, parmetric, permetric, rng, weights):
    print ("Calculating anisotropic RD with weights...\n DR=")
    dr = np.zeros((len(binspar)-1, len(binsper)-1))
    # bt = BallTree(dat, metric='pyfunc', func=permetric)
    for i in tqdm(range(len(datR))):
        ind = arbt.query_radius(datR[i].reshape(1, -1), maxrad)
        for j in ind:
            dist0 = dist.cdist([datR[i], ], dat[j], parmetric)[0]
            dist1 = dist.cdist([datR[i], ], dat[j], permetric)[0]
            dr += np.histogram2d(dist0, dist1, range=rng, bins=(binspar, binsper), weights=weights[j])[0]
    dr[dr == 0] = 1.0
    # DR = dr/(Nd*Nr)
    # print (dr/2.0)
    return dr/2.0


def aDDwcalcp(dat, binspar, binsper, parmetric, permetric, rng, weights, rNd, multi=False, queue=0):
    dd = np.zeros((len(binspar)-1, len(binsper)-1))
    # ddbt = BallTree(dat, metric='pyfunc', func=permetric)
    for i in tqdm(rNd):
        ind = adbt.query_radius(dat[i].reshape(1, -1), maxrad)
        for j in ind:
            dist0 = dist.cdist([dat[i], ], dat[j[j>i]], parmetric)[0]
            dist1 = dist.cdist([dat[i], ], dat[j[j>i]], permetric)[0]
            dd += np.histogram2d(dist0, dist1, range=rng, bins=(binspar, binsper), weights=weights[j[j>i]])[0]
    if multi:
        queue.put(dd)
    else:
        return dd
    # print (DD)
    return dd


def aRRwcalcp(datR, binspar, binsper, parmetric, permetric, rng, rweights, rNr, multi=False, queue=0):
    rr = np.zeros((len(binspar)-1, len(binsper)-1))
    # rrbt = BallTree(datR, metric='pyfunc', func=permetric)
    for i in tqdm(rNr):
        ind = arbt.query_radius(datR[i].reshape(1, -1), maxrad)
        for j in ind:
            dist0 = dist.cdist([datR[i], ], datR[j[j>i]], parmetric)[0]
            dist1 = dist.cdist([datR[i], ], datR[j[j>i]], permetric)[0]
            rr += np.histogram2d(dist0, dist1, range=rng, bins=(binspar, binsper), weights=rweights[j[j>i]])[0]
    if multi:
        queue.put(rr)
    else:
        return rr
    # rr[rr == 0] = 1.0
    # Nr = len(datR)
    # RR = rr/(Nr*(Nr-1.0)) # factor of 2 cancels with 1/2 that needs to be done to remove double counting of pairs
    # print (RR)
    return rr


def aDRwcalcp(dat, datR, binspar, binsper, parmetric, permetric, rng, rweights, rNd, multi=False, queue=0):
    # print ("Calculating anisotropic DR with weights (parallelized)...\n DR=")
    dr = np.zeros((len(binspar)-1, len(binsper)-1))
    # rrbt = BallTree(datR, metric='pyfunc', func=permetric)
    for i in tqdm(rNd):
        ind = arbt.query_radius(dat[i].reshape(1, -1), maxrad)
        for j in ind:
            dist0 = dist.cdist([dat[i], ], datR[j], parmetric)[0]
            dist1 = dist.cdist([dat[i], ], datR[j], permetric)[0]
            dr += np.histogram2d(dist0, dist1, range=rng, bins=(binspar, binsper), weights=rweights[j])[0]
    if multi:
        queue.put(dr)
    else:
        return dr
    # print (DR)
    return dr


def aRDwcalcp(dat, datR, binspar, binsper, parmetric, permetric, rng, weights, rNr, multi=False, queue=0):
    # print ("Calculating anisotropic RD with weights (parallelized)...\n DR=")
    dr = np.zeros((len(binspar)-1, len(binsper)-1))
    # bt = BallTree(dat, metric='pyfunc', func=permetric)
    for i in tqdm(rNr):
        ind = adbt.query_radius(datR[i].reshape(1, -1), maxrad)
        for j in ind:
            dist0 = dist.cdist([datR[i], ], dat[j], parmetric)[0]
            dist1 = dist.cdist([datR[i], ], dat[j], permetric)[0]
            dr += np.histogram2d(dist0, dist1, range=rng, bins=(binspar, binsper), weights=weights[j])[0]
    if multi:
        queue.put(dr)
    else:
        return dr
    return dr


def amulti_autocp(dat, binspar, binsper, parmetric, permetric, rng, weights, Nd, CORES=pcpus):

    DD = np.zeros((len(binspar)-1, len(binsper)-1))
    queues = [RetryQueue() for i in range(CORES)]
    args = [(dat, binspar, binsper, parmetric, permetric, rng, weights, range(int(Nd*i/CORES),int(Nd*(i+1)/CORES)), True, queues[i]) for i in range(CORES)]
    jobs = [Process(target=aDDwcalcp, args=(a)) for a in args]
    for j in jobs: j.start()
    for q in queues: DD += q.get()
    for j in jobs: j.join()
    DD[DD == 0] = 1.0
    # DD = DD/(Nd*(Nd-1.0)) # factor of 2 cancels with 1/2 that needs to be done to remove double counting of pairs
    print (DD)
    return DD


def amulti_autocpr(datR, binspar, binsper, parmetric, permetric, rng, rweights, Nr, CORES=pcpus):

    RR = np.zeros((len(binspar)-1, len(binsper)-1))
    queues = [RetryQueue() for i in range(CORES)]
    args = [(datR, binspar, binsper, parmetric, permetric, rng, rweights, range(int(Nr*i/CORES),int(Nr*(i+1)/CORES)), True, queues[i]) for i in range(CORES)]
    jobs = [Process(target=aRRwcalcp, args=(a)) for a in args]
    for j in jobs: j.start()
    for q in queues: RR += q.get()
    for j in jobs: j.join()
    RR[RR == 0] = 1.0
    # RR = RR/(Nr*(Nr-1.0)) # factor of 2 cancels with 1/2 that needs to be done to remove double counting of pairs
    print (RR)
    return RR


def amulti_crosscp(dat, datR, binspar, binsper, parmetric, permetric, rng, weights, Nr, CORES=pcpus):

    RD = np.zeros((len(binspar)-1, len(binsper)-1))
    queues = [RetryQueue() for i in range(CORES)]
    args = [(dat, datR, binspar, binsper, parmetric, permetric, rng, weights, range(int(Nr*i/CORES), int(Nr*(i+1)/CORES)), True, queues[i]) for i in range(CORES)]
    jobs = [Process(target=aRDwcalcp, args=(a)) for a in args]
    for j in jobs: j.start()
    for q in queues: RD += q.get()
    for j in jobs: j.join()
    RD[RD == 0] = 1.0
    # Nd=len(dat)
    # DR = DR/(Nd*Nr)
    print (RD/2.0)
    return RD/2.0


def amulti_crosscpr(dat, datR, binspar, binsper, parmetric, permetric, rng, rweights, Nd, CORES=pcpus):

    DR = np.zeros((len(binspar)-1, len(binsper)-1))
    queues = [RetryQueue() for i in range(CORES)]
    args = [(dat, datR, binspar, binsper, parmetric, permetric, rng, rweights, range(int(Nd*i/CORES), int(Nd*(i+1)/CORES)), True, queues[i]) for i in range(CORES)]
    jobs = [Process(target=aDRwcalcp, args=(a)) for a in args]
    for j in jobs: j.start()
    for q in queues: DR += q.get()
    for j in jobs: j.join()
    DR[DR == 0] = 1.0
    # Nd=len(dat)
    # DR = DR/(Nd*Nr)
    print (DR/2.0)
    return DR/2.0


def ximonopole(correlsmu, mu):
    xi0 = np.sum(correlsmu*sp.special.legendre(0)(mu),axis=1)/len(mu)
    np.savetxt("xi0.txt",xi0)
    return xi0


def xidipole(correlsmu, mu):
    xi2 = np.sum(5.0*correlsmu*sp.special.legendre(2)(mu),axis=1)/len(mu)
    np.savetxt("xi2.txt",xi2)
    return xi2


def xiquadpole(correlsmu, mu):
    xi4 = np.sum(9.0*correlsmu*sp.special.legendre(4)(mu),axis=1)/len(mu)
    np.savetxt("xi4.txt",xi4)
    return xi4

def beta(correlsmu, mu):
    xis0 = ximonopole(correlsmu,mu)
    xis2 = xidipole(correlsmu,mu)
    xis4 = xiquadpole(correlsmu,mu)
    xir = xis0*sp.special.legendre(0)(mu) + xis2*sp.special.legendre(2)(mu) + xis4*sp.special.legendre(4)(mu)
    r = xir/xis0
    return 5.0/3.0*(np.sqrt(1.8*r-0.8)-1.0)

# def beta(correlsmu, mu):
#     betav =
