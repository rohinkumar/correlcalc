__author__ = 'Rohin Kumar Y'
import os
import astropy.io.ascii as ascii
import astropy.io.fits as fits
from astropy.table import Table
import pymangle


def readinfile(filename, ftype):
    """Method to run basic checks on input data/random files and return contents"""
    if os.path.isfile(filename):
        dat = ascii.read(filename)
        cols = ['z', 'ra', 'dec']
        colnames = []
        for x in dat.colnames:
            colnames.append(x.lower())
        if all(x in colnames for x in cols):
            if ftype.lower() == 'data':
                print ("Entered ascii data file")
                print (dat)
                return dat
            elif ftype.lower() == 'random':
                print ("Entered ascii random file")
                print (dat)
                return dat
            elif ftype.lower() == 'internal':
                return dat
            else:
                print ("Please provide the file as 'data' or 'random'")
                return None
        else:
            print ("File should at least contain columns named 'z', 'ra' and 'dec'. Column names in your ascii file are")
            print(colnames)
            return None
    else:
        print ("Invalid File Path, File Doesn't exist")
        return None


def readfitsfile(fname, ftype):
    """Basic checks for fits file on input data/random files and returns table data from fits"""
    if os.path.isfile(fname):
        hdu_list = fits.open(fname, memmap=True)
        print("Entered fits file containing following data")
        hdu_list.info()
        print(hdu_list[1].columns)
        dat = Table(hdu_list[1].data)
        cols = ['z', 'ra', 'dec']
        colnames = []
        for x in dat.colnames:
            colnames.append(x.lower())
        if all(x in colnames for x in cols):
            if ftype.lower() == 'data':
                print ("Entered fits data file")
                print (dat)
                return dat
            elif ftype.lower() == 'random':
                print ("Entered fits random file")
                print (dat)
                return dat
            elif ftype.lower() == 'internal':
                return dat
            else:
                print ("Please provide the file as 'data' or 'random'")
                return None
        else:
            print ("File should at least contain columns named 'z', 'ra' and 'dec'. Column names in your fits file are")
            print(colnames)
            return None
    else:
        print ("Invalid File Path, File Doesn't exist")
        return None
    # else:
    #         print ("Please provide fits file with .fits extension")
    #         return None


def readmaskfile(fname):
    """Basic checks for mangle ply file. Returns mangle object"""
    if os.path.isfile(fname):
        if fname.lower().endswith('.ply'):
            mangle = pymangle.Mangle(fname)
            print(mangle)
            return mangle
        else:
            print ("Please provide mangle polygon file with .ply extension")
            return None
    else:
        print ("Invalid File Path, File Doesn't exist")
        return None


def storerandcat(z, ra, dec, rweights, rcatname):
    """Create random catalog file to store generated random catalog from randcatprep method"""
    fobj = open(rcatname, 'w')
    fobj.write("z\t ra\t dec\t radial_weight\n ")
    for i in range(0, len(z)):
        fobj.write("%f\t " % z[i])
        fobj.write("%f\t %f\t %f\n" % (ra[i], dec[i], rweights[i]))
    fobj.close()
