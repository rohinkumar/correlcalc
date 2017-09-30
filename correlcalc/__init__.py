import numpy as np
import math as m
from scipy import integrate
from cosmology import *
try:
    from correlcalc.metrics import *
except ImportError:
    print "correlcalc.metrics ImportError"
from metrics import *
from hw import hw
def version():
    print "correlcalc version no. 0.917"
