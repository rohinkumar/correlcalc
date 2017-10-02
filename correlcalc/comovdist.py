__author__ 'Rohin Kumar Y'
from scipy import integrate
import numpy as np 
import math as m
#from . import *
from .param import *
#from correlcalc import *

#Om = param.Om
#Ol = param.Ol

def Ez(zv):
    return 1.0/m.sqrt(Om*(1.0+zv)**3+(1.0-Om-Ol)*(1.0+zv)**2+Ol)

def DC_LCDM(z):
    np.vectorize(Ez)
    return integrate.quad(Ez, 0.0, z)[0]

def DC_LC(z):
    return np.log(1.0+z)

def comov(z,model='lcdm'):
    if model='lcdm':
        return DC_LCDM(z)
    else if model='lc':
        return DC_LC(z)
    else:
        return None
