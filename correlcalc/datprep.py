__author__ = 'Rohin Kumar Y'
import numpy as np
from fileios import *
from genrand import *
from comovdist import *
from math import pi
#calculate comoving distances and convert to 3xN np array
#comov(z,**kwargs)
#**kwargs 'lcdm', 'lc' model

def datprep(fname,ftype,model):
    """Method to convert input ascii files into 3xN matrices needed for metric distance calculations. Calculates comoving distances as per input model."""
    data=readinfile(fname,ftype)
    for x in data.colnames:
        if x.lower()=='z':
            z=np.array(data[x])
        elif x.lower()=='ra':
            ra=np.array(data[x])
        elif x.lower()=='dec':
            dec=np.array(data[x])
        else:
            print("Some problem with your file?")
    s=comov(z,model)
    rar=ra*pi/180.0
    decr=dec*pi/180.0
    dat=np.array([s,rar,decr])
    dat.reshape(3,len(data))
    dat=dat.transpose()
    #print ("Printing 3xN array")
    #print(dat)
    return dat

def randcatprep(z,randcatsize,model):
    """Method to generate random catalog from mangle mask and input redshift distribution"""
    z=randz(z,randcatsize)
    rar,decr=randang(maskfile,randcatsize)
    s=comov(z,model)
    datR=np.array([s,rar,decr])
    datR.reshape(3,len(data))
    datR=datR.transpose()
    print(datR)
    return datR
