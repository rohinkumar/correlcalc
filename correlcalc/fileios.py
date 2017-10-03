__author__ = 'Rohin Kumar Y'
import os
import astropy.io.ascii as ascii
import pymangle

def readinfile(filename,ftype):
    """Method to run basic checks on input data/random files and return contents"""
    if os.path.isfile(filename):
        dat=ascii.read(filename)
        cols=['z','ra','dec']
        colnames=[]
        for x in dat.colnames:
            colnames.append(x.lower())
        if all(x in colnames for x in cols):
            if ftype=='data':
                print ("Entered data ascii file")
                print (dat)
                return dat
            elif ftype== 'random':
                print ("Entered random ascii file")
                print (dat)
                return dat
            elif ftype== 'internal':
                return dat
            else:
                print ("Please provide the file as 'data' or 'random'")
                return None
        else:
            print "File should contain columns named 'z', 'ra' and 'dec'. Column names in your ascii file are"
            print(colnames)
            return None
    else:
        print "Invalid File Path, File Doesn't exist"
        return None

def readmaskfile(fname):
    """Basic checks for mangle ply file. Returns mangle object"""
    if os.path.isfile(fname):
        if fname.lower().endswith('.ply'):
            mangle=pymangle.Mangle(fname)
            print(mangle)
            return mangle
        else:
            print "Please provide mangle polygon file with .ply extension"
            return None
    else:
        print "Invalid File Path, File Doesn't exist"
        return None

def storerandcat(z,ra,dec,rcatname):
    """Create random catalog file to store generated random catalog from randcatprep method"""
    fobj = open(rcatname,'w')
    fobj.write("z\t ra\t dec\n ")
    for i in range(0,len(z)):
        fobj.write("%f\t " %z[i])
        fobj.write("%f\t %f\n " %(ra[i],dec[i]))
    fobj.close()
