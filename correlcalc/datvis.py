import healpix_util as hu
import numpy as np
def healpixmap(ra,dec):
    NSIDE=512
    hpix=hu.HealPix("ring",NSIDE)
    pix=hpix.eq2pix(ra,dec)
    hpixdata=np.array(np.zeros(hu.nside2npix(NSIDE)))
    for j in range(len(pix)):
        hpixdata[pix[j]]+=1
    return hpixdata

