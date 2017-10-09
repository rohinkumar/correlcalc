__author__ = 'Rohin Kumar Y'
#Calculate anisotropic 2pCF
from tpcf import *
#antpcf(dat,datR,bins,parmetric,permetric) returns numpy 2d array DD, RR, DR correl
#poserr(xi,DD) returns (1.0+xi)/np.sqrt(DD)

def atpcf(datfile, bins, **kwargs):
    """Main function to calculate anisotropic 2pCF. Takes multiple arguments such as randfile, maskfile, calculation estimator etc. for different geometry, cosmology models"""
    #Default function arguments
    rng=np.array([[min(bins), max(bins)], [min(bins), max(bins)]])
    weightsflag=False
    cosmology='lcdm'
    #geometry='flat'
    metric=flatdistsq
    randcatfact = 2
    estimator='dp'
    binsq=bins**2
    randfile=None
    maskfile=None

    #Options for correl calculation estimators and cosmology models
    mlist=['dp','ls','ph','hew','h']
    clist=['lcdm','lc']#to add wcdm

    if kwargs is not None:
        for key, value in kwargs.iteritems():
            #print (key, value)
            if key.lower()=='randfile':
                randfile=value

            elif key.lower()=='randfact':
                randcatfact=value

            elif key.lower()=='parmetric':
                if value.lower()=='apdz':
                    parmetric=APdz
                    #metric='APdz'
                # elif value.lower()=='spar':
                #     #metric='spar'
                #     parmetric=spar
                # elif value.lower()=='s':
                #     #geometry='close'
                #     parmetric=s
                else:
                    print("Incorrect parallel metric argument provided!")
            elif key.lower()=='permetric':
                if value.lower()=='apzdth':
                    permetric=APzdth
                    #geometry='flat'
                # elif value.lower()=='sper':
                #     #geometry='open'
                #     permetric=sper
                # elif value.lower()=='mu':
                #     #geometry='close'
                #     permetric=mu
                else:
                    print("Incorrect perpendicular metric provided!")

            elif key.lower()=='estimator':
                if value.lower() in mlist:
                    estimator=value.lower()
                else:
                    print("Incorrect estimator provided! Using 'dp' as default")
            elif key.lower()=='cosmology':
                if value.lower() in clist:
                    cosmology=value.lower()
                else:
                    print("Incorrect Cosmology provided! Using 'lcdm' as default")
            elif key.lower()=='mask':
                maskfile=value
            elif key.lower()=='weights':
                if value==True:
                    weightsflag=True
                else:
                    weightsflag=False
            else:
                print ("key argument `%s` not valid"%key)
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
    #Prepare dat from data file
    dat, weights=datprep(datfile, 'data', cosmology)
    Nd=len(dat)
    #print (weights)
    #Prepare datR from random file or generate a random catalog
    if randfile==None:
        randcatsize=randcatfact*Nd
        if maskfile==None:
            print ("Mask file compulsory. Please provide mask='maskfilepath.ply'")
        else:
            datR=randcatprep(datfile,randcatsize,maskfile,cosmology)
            rweights=np.array([])
            #randfile='./randcat.dat'
            #datR, rweights=datprep(randfile,'random',cosmology)
    else:
        datR, rweights=datprep(randfile,'random',cosmology)

    #if len(weights)!=0:
        #rfdat=readinfile(randfile,ftype='internal')
        #rweights=1.0/(1.0+4.0*np.array(rfdat['nz']))
        #rweights=rweights/np.mean(rweights)
        #print (rweights)
    #Nr=len(datR)

    #Prepare dat from data file
    # dat=datprep(datfile, 'data', cosmology)
    # Nd=len(dat)
    #
    # #Prepare datR from random file or generate a random catalog
    # if randfile==None:
    #     randcatsize=randcatfact*Nd
    #     if maskfile==None:
    #         print ("Mask file compulsory. Please provide mask='maskfilepath.ply'")
    #     else:
    #         datR=randcatprep(datfile,randcatsize,maskfile,cosmology)
    # else:
    #     datR=datprep(randfile,'random',cosmology)

    #Nr=len(datR)
    print ("Calculating anisotropic 2pCF...")

    #f=(1.0*Nrd)/N

    #Reference: arXiv: 1211.6211
    if estimator=='dp':
        if weightsflag==False or len(weights)!=Nd:
            # print (weightsflag)
            # print(len(weights))
            # print(len(datR))
            DD=aDDcalc(dat,binsq,parmetric,permetric,rng)
            DR=aDRcalc(dat,datR,binsq,parmetric,permetric,rng)
        else:
            # if len(rweights)!=len(datR):
            DD=aDDwcalc(dat,binsq,parmetric,permetric,rng,weights)
            DR=aRDwcalc(dat,datR,binsq,parmetric,permetric,rng,weights)
            # else:
            #     DD=aDDwcalc(dat,binsq,parmetric,permetric,rng,weights)
            #     DR=aDRwcalc(dat,datR,binsq,parmetric,permetric,rng,weights,rweights)

        print ("Using Davis-Peebles estimator")
        correl=(DD/DR)-1.0

    elif estimator=='ph':
        if weightsflag==False or len(weights)!=Nd or len(rweights)!=len(datR):
            DD=aDDcalc(dat,binsq,parmetric,permetric,rng)
            RR=aRRcalc(datR,binsq,parmetric,permetric,rng)
        else:
            DD=aDDwcalc(dat,binsq,parmetric,permetric,rng,weights)
            RR=aRRwcalc(datR,binsq,parmetric,permetric,rng,rweights)
        print ("Using Peebles-Hauser estimator")
        correl=(DD/RR)-1.0
    else:
        if weightsflag==False or len(weights)!=Nd or len(rweights)!=len(datR):
            DD=aDDcalc(dat,binsq,parmetric,permetric,rng)
            RR=aDDcalc(datR,binsq,parmetric,permetric,rng)
            DR=aDRcalc(dat,datR,binsq,parmetric,permetric,rng)
        else:
            DD=aDDwcalc(dat,binsq,parmetric,permetric,rng,weights)
            RR=aRRwcalc(datR,binsq,parmetric,permetric,rng,rweights)
            DR=aRDwcalc(dat,datR,binsq,parmetric,permetric,rng,weights)
        if estimator=='ls':
            print ("Using Landy-Szalay estimator")
            correl=(DD-2.0*DR+RR)/RR
        elif estimator=='hew':
            print ("Using Hewett estimator")
            correl=(DD-DR)/RR
        elif estimator=='h':
            print ("Using Hamilton estimator")
            correl=(DD*RR)/DR**2 - 1.0
    correlerr = poserr(correl,DD*Nd*(Nd-1.0)*0.5)
    print("Anisotropic Two-point correlation=")
    print (correl, correlerr)
    return correl, correlerr

def aDDcalc(dat,bins,parmetric,permetric,rng):
    print ("Calculating anisotropic DD...\n DD=")
    dd=np.zeros((len(bins)-1,len(bins)-1))
    ddbt=BallTree(dat,metric='pyfunc',func=permetric)
    for i in tqdm(xrange(len(dat))):
        ind=ddbt.query_radius(dat[i].reshape(1,-1),max(bins))
        for j in ind:
            dist0=dist.cdist([dat[i],],dat[j],parmetric)[0]
            dist1=dist.cdist([dat[i],],dat[j],permetric)[0]
            dd+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins))[0]
    dd[dd==0]=1.0
    Nd=len(dat)
    DD=2.0*dd/(Nd*(Nd-1.0))
    print (DD)
    return DD

def aRRcalc(datR,bins,parmetric,permetric,rng):
    print ("Calculating anisotropic RR...\n RR=")
    rr=np.zeros((len(bins)-1,len(bins)-1))
    rrbt=BallTree(datR,metric='pyfunc',func=permetric)
    for i in tqdm(xrange(len(datR))):
        ind=rrbt.query_radius(datR[i].reshape(1,-1),max(bins))
        for j in ind:
            dist0=dist.cdist([datR[i],],datR[j],parmetric)[0]
            dist1=dist.cdist([datR[i],],datR[j],permetric)[0]
            rr+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins))[0]
    rr[rr==0]=1.0
    Nr=len(datR)
    RR=2.0*rr/(Nr*(Nr-1.0))
    print (RR)
    return RR

def aDRcalc(dat,datR,bins,parmetric,permetric,rng):
    print ("Calculating anisotropic DR...\n DR=")
    dr=np.zeros((len(bins)-1,len(bins)-1))
    rrbt=BallTree(datR,metric='pyfunc',func=permetric)
    for i in tqdm(xrange(len(dat))):
        ind=rrbt.query_radius(dat[i].reshape(1,-1),max(bins))
        for j in ind:
            dist0=dist.cdist([dat[i],],datR[j],parmetric)[0]
            dist1=dist.cdist([dat[i],],datR[j],permetric)[0]
            dr+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins))[0]
    dr[dr==0]=1.0
    Nd=len(dat)
    Nr=len(datR)
    DR=dr/(Nd*Nr)
    print (DR)
    return DR

def aDDwcalc(dat,bins,parmetric,permetric,rng,weights):
    print ("Calculating anisotropic DD with weights...\n DD=")
    dd=np.zeros((len(bins)-1,len(bins)-1))
    ddbt=BallTree(dat,metric='pyfunc',func=permetric)
    for i in tqdm(xrange(len(dat))):
        ind=ddbt.query_radius(dat[i].reshape(1,-1),max(bins))
        for j in ind:
            dist0=dist.cdist([dat[i],],dat[j],parmetric)[0]
            dist1=dist.cdist([dat[i],],dat[j],permetric)[0]
            dd+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins),weights=weights[j])[0]
    dd[dd==0]=1.0
    Nd=len(dat)
    DD=2.0*dd/(Nd*(Nd-1.0))
    print (DD)
    return DD

def aRRwcalc(datR,bins,parmetric,permetric,rng,rweights):
    print ("Calculating anisotropic RR with weights...\n RR=")
    rr=np.zeros((len(bins)-1,len(bins)-1))
    rrbt=BallTree(datR,metric='pyfunc',func=permetric)
    for i in tqdm(xrange(len(datR))):
        ind=rrbt.query_radius(datR[i].reshape(1,-1),max(bins))
        for j in ind:
            dist0=dist.cdist([datR[i],],datR[j],parmetric)[0]
            dist1=dist.cdist([datR[i],],datR[j],permetric)[0]
            rr+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins),weights=rweights[j])[0]
    rr[rr==0]=1.0
    Nr=len(datR)
    RR=2.0*rr/(Nr*(Nr-1.0))
    print (RR)
    return RR

def aDRwcalc(dat,datR,bins,parmetric,permetric,rng,weights,rweights):
    print ("Calculating anisotropic DR with weights...\n DR=")
    dr=np.zeros((len(bins)-1,len(bins)-1))
    rrbt=BallTree(datR,metric='pyfunc',func=permetric)
    for i in tqdm(xrange(len(dat))):
        ind=rrbt.query_radius(dat[i].reshape(1,-1),max(bins))
        for j in ind:
            dist0=dist.cdist([dat[i],],datR[j],parmetric)[0]
            dist1=dist.cdist([dat[i],],datR[j],permetric)[0]
            dr+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins),weights=rweights[j])[0]
    dr[dr==0]=1.0
    Nd=len(dat)
    Nr=len(datR)
    DR=dr/(Nd*Nr)
    print (DR)
    return DR

def aRDwcalc(dat,datR,bins,parmetric,permetric,rng,weights):
    print ("Calculating anisotropic RD with weights...\n DR=")
    dr=np.zeros((len(bins)-1,len(bins)-1))
    bt=BallTree(dat,metric='pyfunc',func=permetric)
    for i in tqdm(xrange(len(datR))):
        ind=bt.query_radius(datR[i].reshape(1,-1),max(bins))
        for j in ind:
            dist0=dist.cdist([datR[i],],dat[j],parmetric)[0]
            dist1=dist.cdist([datR[i],],dat[j],permetric)[0]
            dr+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins),weights=weights[j])[0]
    dr[dr==0]=1.0
    Nd=len(dat)
    Nr=len(datR)
    DR=dr/(Nd*Nr)
    print (DR)
    return DR
# def antpcf(dat,datR,bins,parmetric,permetric,estimator,**kwargs):
#     rng=np.array([[min(bins), max(bins)], [min(bins), max(bins)]])
#     print "Calculating anisotropic DD..."
#     dd=np.zeros((len(bins),len(bins))
#     ddbt=BallTree(dat,metric='pyfunc',func=permetric)
#     for i in tqdm(xrange(len(dat))):
#         ind=ddbt.query_radius(dat[i].reshape(1,-1),max(bins))
#         for j in ind:
#             dist0=dist.cdist([dat[i],],dat[j],parmetric)[0]
#             dist1=dist.cdist([dat[i],],dat[j],permetric)[0]
#             dd+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins))[0]
#     print dd
#
#     print "Calculating anisotropic RR..."
#     rr=np.zeros((len(bins),len(bins)))
#     rrbt=BallTree(datR,metric='pyfunc',func=permetric)
#     for i in tqdm(xrange(len(datR))):
#         ind=rrbt.query_radius(datR[i].reshape(1,-1),max(bins))
#         for j in ind:
#             dist0=dist.cdist([datR[i],],datR[j],parmetric)[0]
#             dist1=dist.cdist([datR[i],],datR[j],permetric)[0]
#             rr+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins))[0]
#     print rr
#
#     print "Calculating anisotropic DR..."
#     dr=np.zeros((len(bins),len(bins)))
#     for i in tqdm(xrange(len(dat))):
#         ind=rrbt.query_radius(dat[i].reshape(1,-1),max(bins))
#         for j in ind:
#             dist0=dist.cdist([dat[i],],datR[j],parmetric)[0]
#             dist1=dist.cdist([dat[i],],datR[j],permetric)[0]
#             dr+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins))[0]
#     print dr
#
#     print "Calculating anisotropic 2pCF with Poisson error"
#     rr[rr==0]=1.0
#     dd[dd==0]=1.0
#     Nrd=len(datR)
#     N=len(dat)
#     f=(1.0*Nrd)/N
#     if estimator=='ls':
#         correl=1.0+f**2*dd/rr-2.0*f*dr/rr
#     elif estimator=='simple':
#         correl=f**2*dd/rr-1.0
#
#     correlerr = poserr(correl,dd)
#     print(correl, correlerr)
#     return correl, correlerr


#def poserr2d(xi,dd2d):
#    return (1.0+xi)/np.sqrt(dd2d)
