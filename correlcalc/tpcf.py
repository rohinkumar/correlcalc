#Samsiddhir Haritoshanam!
__author__ = 'Rohin Kumar Y'
#tpcf(dat, datR=None, randcatsize=2, bins,**kwargs)
#
#**kwargs for choosing geometry - metric 'flat' 'open' 'close'
#**kwargs for choosing xi method - 'simple' 'ls' '...'
#import fileios
from tqdm import *
from dataprep import *
import numpy as np
from sklearn.neighbors import BallTree
from __future__ import divison

def tpcf(dat, datR, randcatfact=2, bins, **kwargs):
    print "Calculating DD..."
    DD=autocorr(dat,bins,metric)
    print DD
    if datR==None:
        randcatsize=randcatfact*len(dat)
        datR=randcatprep(z,randcatsize,model)
    print "Calculating DR..."
    DR=crosscorr(dat,datR,bins,metric)
    print DR
    print "Calculating RR..."
    RR=autocorr(datR,bins,metric)
    print RR

    print "Calculating 2pCF with Poisson error"

    Nrd=len(datR)
    N=len(dat)
    f=(1.0*Nrd)/N
    if method='ls':
        correl=1.0+f**2*DD/RR-2.0*f*DR/RR
    else if method='simple':
        correl=f**2*DD/RR-1.0

    correlerr = poserr(correl,DD)

    return correl, correlerr

def autocorr(dat,bins,metric):
    bt=BallTree(dat,metric='pyfunc',func=metric)
    counts_DD=bt.two_point_correlation(dat,bins)
    DD = np.diff(counts_DD)
    return DD

def crosscorr(dat,datR,bins,metric,**kwargs):
    bt=BallTree(dat,metric='pyfunc',func=metric)
    counts_DR=bt.two_point_correlation(datR,bins)
    DR=np.diff(counts_DR)
    return DR

def poserr(xi,DD):
    return (1.0+xi)/np.sqrt(DD)
#alternatively
#rbt=BallTree(dat,metric='pyfunc',func=metric)
#counts_RD=rbt.two_point_correlation(dat,bins)
