__author__ = 'Rohin Kumar Y'
import numpy as np
from fileios import *
from genrand import *
from comovdist import *
from math import pi
#calculate comoving distances and convert to 3xN np array
#comov(z,**kwargs)
#**kwargs 'lcdm', 'lc' model

def datprep (fname,ftype,model):
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
            pass
    s=comov(z,model)
    rar=ra*pi/180.0
    decr=dec*pi/180.0
    dat=np.array([s,rar,decr])
    dat.reshape(3,len(data))
    dat=dat.transpose()
    #print ("Printing 3xN array")
    #print(dat)
    return dat


def randcatprep(datfname,randcatsize,maskfile,model):
    """Method to generate random catalog from mangle mask and input redshift distribution"""
    print("Generating random catalog of %d size in file randcat.dat... "%randcatsize)
    data=readinfile(datfname,'internal')
    for x in data.colnames:
        if x.lower()=='z':
            z=np.array(data[x])
    zr=randz(z,randcatsize)
    ra,dec=randang(maskfile,randcatsize)
    rar=ra*pi/180.0
    decr=dec*pi/180.0
    rcatfname="randcat.dat"#%(datfname)
    storerandcat(z,rar,decr,rcatfname)
    s=comov(zr,model)
    datR=np.array([s,rar,decr])
    datR.reshape(3,len(zr))
    datR=datR.transpose()
    #print(datR)
    return datR
