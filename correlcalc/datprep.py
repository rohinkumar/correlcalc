__author__ = 'Rohin Kumar Y'
import numpy as np
from fileios import *
from genrand import *
from comovdist import *
from math import pi
#calculate comoving distances and convert to 3xN np array
#comov(z,**kwargs)
#**kwargs 'lcdm', 'lc' model

def datprep(ftype,model):
    filename=inputfiles(ftype)
    data=ascii.read(filename)
    z=np.array(data['z'])
    ra=np.array(data['ra'])
    dec=np.array(data['dec'])
    s=comov(z,model)
    rar=ra*pi/180.0
    decr=dec*pi/180.0
    dat=np.array([s,rar,decr])
    dat.reshape(3,len(data))
    dat=dat.transpose()
    return dat

def randcatprep(z,randcatsize,model):
    maskfile=inputmaskfile()
    z=randz(z,randcatsize)
    rar,decr=randang(maskfile,randcatsize)
    s=comov(z,model)
    datR=np.array([s,rar,decr])
    datR.reshape(3,len(data))
    datR=datR.transpose()
    return datR
