import healpix_util as hu
import numpy as np
# import matplotlib.pyplot as plt


def healpixmap(ra, dec, weights):
    """Method to view/visualize angular distribution of galaxies in terms of healpix maps"""
    NSIDE = 512
    hpix = hu.HealPix("ring", NSIDE)
    pix = hpix.eq2pix(ra, dec)
    hpixdata = np.array(np.zeros(hu.nside2npix(NSIDE)))
    for j in range(len(pix)):
        hpixdata[pix[j]] += weights[j]
    hu.mollview(hpixdata, rot=180)
    return hpixdata
