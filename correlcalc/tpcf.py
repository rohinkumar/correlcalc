from __future__ import division
#Samsiddhir Haritoshanam!
__author__ = 'Rohin Kumar Y'

#tpcf(dat, datR=None, randcatsize=2, bins,**kwargs)
#
#**kwargs for choosing geometry - metric 'flat' 'open' 'close'
#**kwargs for choosing xi method - 'simple' 'ls' '...'
#import fileios
from tqdm import *
from datprep import *
import numpy as np
from metrics import *
from sklearn.neighbors import BallTree


def tpcf(datfile, bins, **kwargs):
    """Main function to calculate 2pCF"""

    #Default function arguments
    cosmology='lcdm'
    geometry='flat'
    metric=flatdistsq
    randcatfact = 2
    method='dp'
    binsq=bins**2
    randfile=None
    maskfile=None

    #Options for correl calculation methods and cosmology models
    mlist=['dp','ls','ph','hew','h']
    clist=['lcdm','lc']

    if kwargs is not None:
        for key, value in kwargs.iteritems():
            #print (key, value)
            if key.lower()=='randfile':
                randfile=value

            elif key.lower()=='randfact':
                randcatfact=value

            elif key.lower()=='geometry':
                if value.lower()=='flat':
                    geometry='flat'
                    metric=flatdistsq
                elif value.lower()=='open':
                    geometry='open'
                    metric=opendistsq
                elif value.lower()=='close':
                    geometry='close'
                    metric=closedistsq
                else:
                    print("Incorrect geometry argument provided! Using flat geometry")

            elif key.lower()=='method':
                if value.lower() in mlist:
                    method=value.lower()
                else:
                    print("Incorrect method provided! Using 'dp' as default")
            elif key.lower()=='cosmology':
                if value.lower() in clist:
                    cosmology=value.lower()
                else:
                    print("Incorrect Cosmology provided! Using 'lcdm' as default")
            elif key.lower()=='mask':
                maskfile=value
            else:
                print ("key argument not valid")
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
    print("Correl method=")
    print(method)
    print("---------------")
    #Prepare dat from data file
    dat=datprep(datfile, 'data', cosmology)
    Nd=len(dat)

    #Prepare datR from random file or generate a random catalog
    if randfile==None:
        randcatsize=randcatfact*Nd
        if maskfile==None:
            print ("Mask file compulsory. Please provide mask='maskfilepath.ply'")
        else:
            datR=randcatprep(datfile,randcatsize,maskfile,cosmology)
    else:
        datR=datprep(randfile,'random',cosmology)

    #Nr=len(datR)


    print ("Calculating 2pCF...")

    #f=(1.0*Nrd)/N

    #Reference: arXiv: 1211.6211

    if method=='ls':
        print ("Using Landy-Szalay method")
        DD=DDcalc(dat,binsq,metric)
        RR=RRcalc(datR,binsq,metric)
        DR=DRcalc(dat,datR,binsq,metric)
        correl=1.0+(DD-2.0*DR)/RR

    elif method=='ph':
        print ("Using Peebles-Hauser method")
        DD=DDcalc(dat,binsq,metric)
        RR=RRcalc(datR,binsq,metric)
        correl=(DD/RR)-1.0

    elif method=='hew':
        print ("Using Hewett method")
        DD=DDcalc(dat,binsq,metric)
        RR=RRcalc(datR,binsq,metric)
        DR=DRcalc(dat,datR,binsq,metric)
        correl=(DD-DR)/RR

    elif method=='dp':
        print ("Using Davis-Peebles method")
        DD=DDcalc(dat,binsq,metric)
        DR=DRcalc(dat,datR,binsq,metric)
        correl=(DD/DR)-1.0

    elif method=='h':
        print ("Using Hamilton method")
        DD=DDcalc(dat,binsq,metric)
        RR=RRcalc(datR,binsq,metric)
        correl=(DD*RR)/DR**2 - 1.0

    correlerr = poserr(correl,DD)
    print("Two-point correlation=")
    print (correl, correlerr)
    return correl, correlerr

def DDcalc(dat,bins,metric):
    print ("Calculating DD...\n DD=")
    DD=autocorr(dat,bins,metric)
    DD[DD==0]=1.0
    Nd=len(dat)
    DD=2.0*DD/(Nd*(Nd-1.0))
    print (DD)
    return DD

def RRcalc(datR,bins,metric):
    print ("Calculating RR...\n RR=")
    RR=autocorr(datR,bins,metric)
    RR[RR==0]=1.0
    Nr=len(datR)
    RR=2.0*RR/(Nr*(Nr-1.0))
    print (RR)
    return RR

def DRcalc(dat,datR,bins,metric):
    print ("Calculating DR...\n DR=")
    DR=crosscorr(dat,datR,bins,metric)
    DR[DR==0]=1.0
    Nd=len(dat)
    Nr=len(datR)
    DR=DR/(Nd*Nr)
    print (DR)
    return DR

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
    print ("Calculating Poisson error")
    return (1.0+xi)/np.sqrt(DD)
#alternatively
#rbt=BallTree(dat,metric='pyfunc',func=metric)
#counts_RD=rbt.two_point_correlation(dat,bins)
