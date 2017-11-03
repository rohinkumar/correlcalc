__author__ = 'Rohin Kumar Y'
import numpy as np
from fileios import *
from genrand import *
from comovdist import *
from math import pi
# calculate comoving distances and convert to 3xN np array
# comov(z,**kwargs)
# **kwargs 'lcdm', 'lc' model


def datprep(fname, ftype, model):
    """Method to convert input ascii files into 3xN matrices needed for metric distance calculations. Calculates comoving distances as per input model."""
    if fname.lower().endswith('.fits'):
        data = readfitsfile(fname, ftype)
    else:
        data = readinfile(fname, ftype)
    weights = np.array([])
    print ("Reading column data from input file...")
    for x in data.colnames:
        if x.lower() == 'z':
            z = np.array(data[x])
        elif x.lower() == 'ra':
            ra = np.array(data[x])
        elif x.lower() == 'dec':
            dec = np.array(data[x])
        elif x.lower() == 'radial_weight':
            weights = np.array(data[x])
        elif x.lower() == 'weight_fkp':
            weights = np.array(data[x])
        elif x.lower() == 'weight_systot':
            weights = np.array(data[x])
        # elif x.lower()=='nz':
        #     weights=1.0/(1.0+4.0*np.array(data[x]))
        #     weights=weights/np.mean(weights)
            # print (weights)
        else:
            pass
    # print ("Calcuating comoving distances and converting angles to radians...")
    s = comovp(z, model)
    rar = ra*pi/180.0
    decr = dec*pi/180.0
    print ("Preparing %s into 3xN matrices in [s,rar,decr] format..." % ftype)
    dat = np.array([s, rar, decr, z])
    dat.reshape(4, len(data))
    dat = dat.transpose()
    # print ("Printing 3xN array")
    # print(dat)
    return dat, weights


def datprepz(fname, ftype, model):
    """Method to convert input ascii files into 3xN matrices needed for metric distance calculations. passes z instead of s in matrices."""
    if fname.lower().endswith('.fits'):
        data = readfitsfile(fname, ftype)
    else:
        data = readinfile(fname, ftype)
    weights = np.array([])
    print ("Reading column data from input file...")
    for x in data.colnames:
        if x.lower() == 'z':
            z = np.array(data[x])
        elif x.lower() == 'ra':
            ra = np.array(data[x])
        elif x.lower() == 'dec':
            dec = np.array(data[x])
        elif x.lower() == 'radial_weight':
            weights = np.array(data[x])
        elif x.lower() == 'weight_fkp':
            weights = np.array(data[x])
        elif x.lower() == 'weight_systot':
            weights = np.array(data[x])
        # elif x.lower()=='nz':
        #     weights=1.0/(1.0+4.0*np.array(data[x]))
        #     weights=weights/np.mean(weights)
            # print (weights)
        else:
            pass
    # print ("Converting angles to radians...")
    rar = ra*pi/180.0
    decr = dec*pi/180.0
    print ("Preparing %s into 3xN matrices in [z,rar,decr] format..." % ftype )
    dat = np.array([z, rar, decr])
    dat.reshape(3, len(data))
    dat = dat.transpose()
    # print ("Printing 3xN array")
    # print(dat)
    return dat, weights


def randcatprep(datfname, randcatsize, maskfile, model):
    """Method to generate random catalog from mangle mask and input redshift distribution for given cosmology"""
    print("Generating random catalog of %d size in file randcat.dat... "% randcatsize)
    if datfname.lower().endswith('.fits'):
        data = readfitsfile(datfname, 'internal')
    else:
        data = readinfile(datfname, 'internal')

    for x in data.colnames:
        if x.lower() == 'z':
            z = np.array(data[x])
        # elif x.lower()=='nz':
        #     rweights=1.0/(1.0+4.0*np.array(data[x]))
        #     rweights=weights/np.mean(rweights)
        # else:
        #
    zr, rweights = randz(z, randcatsize)
    ra, dec = randang(maskfile, randcatsize)
    print ("Generated random catalog...")
    print ("zr = ")
    print (zr)
    print ("random ra =")
    print (ra)
    print ("random dec =")
    print (dec)
    print ("radial weights =")
    print (rweights)
    rar = ra*pi/180.0
    decr = dec*pi/180.0
    rcatfname = "randcat.dat"# %(datfname)
    storerandcat(zr, ra, dec, rweights, rcatfname)
    s = comovp(zr, model)
    print ("Preparing random catalog into 3xN matrices in [s,rar,decr] format" )
    datR = np.array([s, rar, decr, zr])
    datR.reshape(4, len(zr))
    datR = datR.transpose()
    # print(datR)
    return datR, rweights


def randcatprepz(datfname, randcatsize, maskfile, model):
    """Method to generate random catalog from mangle mask and input redshift distribution for given cosmology"""
    print("Generating random catalog of %d size in file randcat.dat... "% randcatsize)
    if datfname.lower().endswith('.fits'):
        data = readfitsfile(datfname, 'internal')
    else:
        data = readinfile(datfname, 'internal')

    for x in data.colnames:
        if x.lower() == 'z':
            z = np.array(data[x])
        # elif x.lower()=='nz':
        #     rweights=1.0/(1.0+4.0*np.array(data[x]))
        #     rweights=weights/np.mean(rweights)
        # else:
        #
    zr, rweights = randz(z, randcatsize)
    ra, dec = randang(maskfile, randcatsize)
    rar = ra*pi/180.0
    decr = dec*pi/180.0
    rcatfname = "randcat.dat"# %(datfname)
    storerandcat(zr, ra, dec, rweights, rcatfname)
    print ("Preparing random catalog into 3xN matrices in [zr,rar,decr] format" )
    datR = np.array([zr, rar, decr])
    datR.reshape(3, len(zr))
    datR = datR.transpose()
    # print(datR)
    return datR, rweights

