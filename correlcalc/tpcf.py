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
    """Main function to calculate 2pCF. Takes multiple arguments such as randfile, maskfile, calculation estimator etc. for different geometry, cosmology models"""
    # Default function arguments
    weights=np.array([])
    weightsflag=False
    cosmology='lcdm'
    geometry='flat'
    metric=flatdistsq
    randcatfact = 2
    estimator='dp'
    binsq=bins**2
    randfile=None
    maskfile=None

    #Options for correl calculation estimators and cosmology models
    elist=['dp','ls','ph','hew','h']
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
            elif key.lower()=='estimator':
                if value.lower() in elist:
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
                if value is True:
                    weightsflag=True
                    # fdat=readinfile(datfile,ftype='internal')
                    # weights=1.0/(1.0+4.0*np.array(fdat['nz']))
                    # weights=weights/np.mean(weights)
                    #print (weights)
                else:
                    weightsflag=False
            else:
                print ("key argument `%s` not valid"%key)
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


    print ("Calculating 2pCF...")

    #f=(1.0*Nrd)/N
    #print (weights)
    #Reference: arXiv: 1211.6211
    if estimator=='dp':
        if weightsflag==False or len(weights)!=Nd:
            # print (weightsflag)
            # print(len(weights))
            # print(len(datR))
            DD=DDcalc(dat,binsq,metric)
            DR=DRcalc(dat,datR,binsq,metric)
        else:
            # if len(rweights)!=len(datR):
            DD=DDwcalc(dat,binsq,metric,weights)
            DR=RDwcalc(dat,datR,binsq,metric,weights)
            # else:
            #     DD=DDwcalc(dat,binsq,metric,weights)
            #     DR=DRwcalc(dat,datR,binsq,metric,rweights)
        print ("Using Davis-Peebles estimator")
        correl=(DD/DR)-1.0

    elif estimator=='ph':
        if weightsflag==False or len(weights)!=Nd:
            DD=DDcalc(dat,binsq,metric)
            RR=RRcalc(datR,binsq,metric)
        else:
            DD=DDwcalc(dat,binsq,metric,weights)
            RR=RRwcalc(datR,binsq,metric,rweights)
        print ("Using Peebles-Hauser estimator")
        correl=(DD/RR)-1.0
    else:
        if weightsflag==False or len(weights)!=Nd:
            DD=DDcalc(dat,binsq,metric)
            RR=RRcalc(datR,binsq,metric)
            DR=DRcalc(dat,datR,binsq,metric)
        else:
            DD=DDwcalc(dat,binsq,metric,weights)
            RR=RRwcalc(datR,binsq,metric,rweights)
            DR=RDwcalc(dat,datR,binsq,metric,weights)
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

def crosscorr(dat,datR,bins,metric):
    rbt=BallTree(datR,metric='pyfunc',func=metric)
    counts_DR=rbt.two_point_correlation(dat,bins)
    DR=np.diff(counts_DR)
    return DR

def poserr(xi,DD):
    print ("Calculating Poisson error")
    return (1.0+xi)/np.sqrt(DD)
#alternatively
#rbt=BallTree(dat,metric='pyfunc',func=metric)
#counts_RD=rbt.two_point_correlation(dat,bins)
def DDwcalc(dat,bins,metric,weights):
    print ("Calculating DD with weights...\n DD=")
    DD=autocorrw(dat,bins,metric,weights)
    DD[DD==0]=1.0
    Nd=len(dat)
    DD=2.0*DD/(Nd*(Nd-1.0))
    print (DD)
    return DD

def RRwcalc(datR,bins,metric,weights):
    print ("Calculating RR with weights...\n RR=")
    RR=autocorrw(datR,bins,metric,weights)
    RR[RR==0]=1.0
    Nr=len(datR)
    RR=2.0*RR/(Nr*(Nr-1.0))
    print (RR)
    return RR

def DRwcalc(dat,datR,bins,metric,rweights):
    print ("Calculating DR with weights...\n DR=")
    DR=crosscorrw(dat,datR,bins,metric,rweights)
    DR[DR==0]=1.0
    Nd=len(dat)
    Nr=len(datR)
    DR=DR/(Nd*Nr)
    print (DR)
    return DR

def RDwcalc(dat,datR,bins,metric,weights):
    print ("Calculating RD with weights...\n RD=")
    DR=crosscorrwrd(dat,datR,bins,metric,weights)
    DR[DR==0]=1.0
    Nd=len(dat)
    Nr=len(datR)
    DR=DR/(Nd*Nr)
    print (DR)
    return DR

def autocorrw(dat,bins,metric,weights):
    bt=BallTree(dat,metric='pyfunc',func=metric)
    DD=np.zeros(len(bins)-1)
    for i in tqdm(xrange(len(dat))):
        ind=bt.query_radius(dat[i].reshape(1,-1),max(bins))
        #wts=np.array([])
        for j in ind:
            dist0=dist.cdist([dat[i],],dat[j],metric)[0]
            DD+=np.histogram(dist0,bins=bins,weights=weights[j])[0]
            #print (dist0,weights[j])
    return DD

def crosscorrw(dat,datR,bins,metric,rweights):
    rbt=BallTree(datR,metric='pyfunc',func=metric)
    DR=np.zeros(len(bins)-1)
    for i in tqdm(xrange(len(dat))):
        ind=rbt.query_radius(dat[i].reshape(1,-1),max(bins))
        #wts=np.array([])
        for j in ind:
            dist0=dist.cdist([dat[i],],datR[j],metric)[0]
            DR+=np.histogram(dist0,bins=bins,weights=rweights[j])[0]
            #print (dist0,weights[j])
    return DR

def crosscorrwrd(dat,datR,bins,metric,weights):
    bt=BallTree(dat,metric='pyfunc',func=metric)
    RD=np.zeros(len(bins)-1)
    for i in tqdm(xrange(len(datR))):
        ind=bt.query_radius(datR[i].reshape(1,-1),max(bins))
        #wts=np.array([])
        for j in ind:
            dist0=dist.cdist([datR[i],],dat[j],metric)[0]
            RD+=np.histogram(dist0,bins=bins,weights=weights[j])[0]
            #print (dist0,weights[j])
    return RD
