__author__ = 'Rohin Kumar Y'
import os
import astropy.io.ascii as ascii

def checkdatfile(filename):
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

def inputfiles(ftype):
    if ftype=='data':
        filename=input("Enter ascii data file path:")
        print filename
        filename=checkdatfile(filename)
    if ftype=='random':
        filename=input("Enter random ascii file path:")
        filename=checkdatfile(filename)
    else:
        print "Some problem with your file?"
        filename=None
    return filename

def inputmaskfile():
    maskfile=raw_input("Please provide mangle mask file:")
    if os.path.isfile(maskfile):
        if filename.lower().endswith('.ply'):
            return maskfile
        else:
            print "Provide mask file with .ply extension"
            return None
    else:
        print "Invalid File Path, File Doesn't exist"
        return None
