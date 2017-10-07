__author__ = 'Rohin Kumar Y'
from scipy.stats import gaussian_kde
import numpy as np
from fileios import *

def kde(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scipy"""
    kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)

def generate_rand_from_pdf(pdf, x_grid, N):
    """Method to generate 'N' no. of random numbers from input probability distribution function (pdf) in form of kernel density"""
    cdf = np.cumsum(pdf)
    cdf = cdf / cdf[-1]
    values = np.random.rand(N)
    value_bins = np.searchsorted(cdf, values)
    random_from_cdf = x_grid[value_bins]
    return random_from_cdf

def randang(maskfile,randcatsize):
    """Method to calculate RA and DEC from mangle .ply file"""
    mangle=readmaskfile(maskfile)
    rar,decr=mangle.genrand(randcatsize)
    return rar,decr

def randz(z,randcatsize):
    """Method to calculate random redshift values from input redshift distribution"""
    x_grid=np.linspace(min(z),max(z),randcatsize)
    pdf=kde(z,x_grid,bandwidth=1e-3)
    return generate_rand_from_pdf(pdf, x_grid, randcatsize)
