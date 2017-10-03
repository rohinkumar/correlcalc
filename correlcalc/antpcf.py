__author__ = 'Rohin Kumar Y'
#Calculate anisotropic 2pCF
#import tpcf
#import numpy as np
from tpcf import *
from scipy.spatial import distance as dist
#antpcf(dat,datR,bins,parmetric,permetric) returns numpy 2d array DD, RR, DR correl
#poserr(xi,DD) returns (1.0+xi)/np.sqrt(DD)

def antpcf(dat,datR,bins,parmetric,permetric,method,**kwargs):
    rng=np.array([[min(bins), max(bins)], [min(bins), max(bins)]])
    print "Calculating anisotropic DD..."
    dd=np.zeros((20,20))
    ddbt=BallTree(dat,metric='pyfunc',func=permetric)
    for i in tqdm(xrange(len(dat))):
        ind=ddbt.query_radius(dat[i].reshape(1,-1),max(bins))
        for j in ind:
            dist0=dist.cdist([dat[i],],dat[j],parmetric)[0]
            dist1=dist.cdist([dat[i],],dat[j],permetric)[0]
            dd+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins))[0]
    print dd

    print "Calculating anisotropic RR..."
    rr=np.zeros((20,20))
    rrbt=BallTree(datR,metric='pyfunc',func=permetric)
    for i in tqdm(xrange(len(datR))):
        ind=rrbt.query_radius(datR[i].reshape(1,-1),max(bins))
        for j in ind:
            dist0=dist.cdist([datR[i],],datR[j],parmetric)[0]
            dist1=dist.cdist([datR[i],],datR[j],permetric)[0]
            rr+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins))[0]
    print rr

    print "Calculating anisotropic DR..."
    dr=np.zeros((20,20))
    for i in tqdm(xrange(len(dat))):
        ind=rrbt.query_radius(dat[i].reshape(1,-1),max(bins))
        for j in ind:
            dist0=dist.cdist([dat[i],],datR[j],parmetric)[0]
            dist1=dist.cdist([dat[i],],datR[j],permetric)[0]
            dr+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins))[0]
    print dr

    print "Calculating anisotropic 2pCF with Poisson error"
    rr[rr==0]=1.0
    dd[dd==0]=1.0
    Nrd=len(datR)
    N=len(dat)
    f=(1.0*Nrd)/N
    if method=='ls':
        correl=1.0+f**2*dd/rr-2.0*f*dr/rr
    elif method=='simple':
        correl=f**2*dd/rr-1.0

    correlerr = poserr(correl,dd)
    print(correl, correlerr)
    return correl, correlerr


#def poserr2d(xi,dd2d):
#    return (1.0+xi)/np.sqrt(dd2d)
