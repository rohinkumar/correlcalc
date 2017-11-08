__author__ = 'Rohin Kumar Y'
from scipy.stats import gaussian_kde
import numpy as np
from fileios import *

# Need to parallelize these... Very slow indeed! :(


def kde(x, x_grid, bandwidth=0.2):
    """Kernel Density Estimation with Scipy"""
    kdev = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1))
    return kdev.evaluate(x_grid)


def generate_rand_from_pdf(pdf, x_grid, N):
    """Method to generate 'N' no. of random numbers from input probability distribution function (pdf) in form of kernel density"""
    cdf = np.cumsum(pdf)
    cdf = cdf / cdf[-1]
    values = np.random.rand(N)
    value_bins = np.searchsorted(cdf, values)
    random_from_cdf, nz = x_grid[value_bins], pdf[value_bins]
    print (nz)
    return random_from_cdf, nz


def randang(maskfile, randcatsize):
    """Method to calculate RA and DEC from mangle .ply file"""
    mangle = readmaskfile(maskfile)
    rar, decr = mangle.genrand(randcatsize)
    rar = np.asfarray(rar)
    decr = np.asfarray(decr)
    return rar, decr


def randz(z, randcatsize):
    """Method to calculate random redshift values from input redshift distribution"""
    x_grid = np.linspace(min(z), max(z), num=randcatsize)
    pdf = kde(z, x_grid, bandwidth=1e-3)
    randzv, nz = generate_rand_from_pdf(pdf, x_grid, randcatsize)
    # randzv = np.array(randzv, np.double)
    rweights = 1.0/(1.0+10000*nz)
    rweights = rweights/np.mean(rweights)
    return randzv, rweights
