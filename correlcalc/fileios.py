__author__ 'Rohin Kumar Y'
import os
import astropy.io.ascii as ascii

def checkfile(filename):
    if os.path.isfile(filename):
        dat=ascii.read(filename)
        cols=['z','ra','dec']
        if all(x in dat.colnames for x in cols):
            return filename
        else:
            print "File should contain columns named 'z', 'ra' and 'dec'"
            return None
    else:
        print "Invalid File Path, File Doesn't exist"
        return None

def infiles(ftype='data'):
    if ftype='data':
        filename=raw_input("Enter ascii data file path:")
        filename=checkfile(filename)
    else if ftype='random':
        filename=raw_input("Enter random ascii file path:")
        filename=checkfile(filename)
    else:
        print "Some problem with your file?"
        filename=None
    return filename
