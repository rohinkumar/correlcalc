__author__ = 'Rohin Kumar Y'
from scipy import integrate
import numpy as np
import math as m
from metrics import *
#from . import *
from param import *
# from correlcalc import *
from multiprocessing import Pool
from multiprocessing import cpu_count
# from __future__ import division
# Om = param.Om
# Ol = param.Ol

pcpus = cpu_count() - 1


def Ezs(zv):
    """Hubble parameter in H0 units"""
    return 1.0/m.sqrt(Om*(1.0+zv)**3+(1.0-Om-Ol)*(1.0+zv)**2+Ol)


Ez = np.vectorize(Ezs)


def DC_LCDMs(z):
    """Method to calculate comoving distance for LCDM model"""
    return integrate.quad(Ez, 0, z)[0]


DC_LCDM = np.vectorize(DC_LCDMs)


def DC_LC(z):
    """Method for comoving distance in Linear Coasting model"""
    return np.log(1.0+z)


def comov(z, model):
    """Method to calculate comoving distance of given redshifts for input model. Units in c/H0"""
    # More models such as wcdm to be added
    print ("Calculating comoving distances...")
    if model == 'lcdm':
        return DC_LCDM(z)
    elif model == 'lc':
        return DC_LC(z)
    else:
        print("Only 'lcdm' and 'lc' models supported for now")
        return None


def comovp(z, model):
    """Method to calculate comoving distance of given redshifts for input model. Units in c/H0"""
    # More models such as wcdm to be added
    print ("Calculating comoving distances (parallelized)...")
    pool = Pool(processes=pcpus)
    if model == 'lcdm':
        res = pool.map(DC_LCDM, z)
        pool.close()
        pool.join()
        return res
    elif model == 'lc':
        res = pool.map(DC_LC, z)
        pool.close()
        pool.join()
        return res
    else:
        print("Only 'lcdm' and 'lc' models supported for now")
        return None


# def sparsq(dat1, dat2, cosmology):
#     if cosmology == 'lc':
#         Omv = 0
#         Olv = 0
#     else:
#         Omv = Om
#         Olv = Ol
#     return csparsq(dat1, dat2, Omv, Olv)


