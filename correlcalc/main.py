__author__ = 'Rohin Kumar Y'
#from fileios import *
#msg = 'Enter Absolute Path to file: '
#f_name = raw_input(msg).strip()

#path = file_data_and_path(f_name)
#if path != None:
#       print 'Path:',path
#import tpcf
#from fileios import *
#from comovdist import *
from datprep import *
from tpcf import *
import matplotlib.pyplot as plt
from genrand import *
from runtimeupdate import *
#from metrics import *
# t1 = checkdatfile('./testfile.dat')
# print t1
#
# t2=inputfiles('data')
# print t2
#%matplotlib osx
from datvis import *

#dat=checkinfile('./testfile.dat','data')
#print dat
from metrics import *
from tpcf import *
#dat=readinfile('./testfile.dat','data')
#dat=datprep('./testfile.dat','data','lcdm')
from antpcf import *
bins=np.arange(0.002,0.082,0.002)
#tpcf('./testfile.dat',bins,randfile='./testfile.dat',method='ls')
tpcf('./testfile.dat',bins,mask='/Users/rohin/Documents/ipy_notebooks/galsurveystudy/masks/boss_geometry_2011_06_10.ply',cosmology='lc',method='ls')
#bins=np.arange(0,0.201,0.01)
#tpcf(dat,dat,3,bins,flatdistsq,'ls')
#antpcf(dat,dat,bins,flatdistsq,flatdistsq,'ls')
#z=np.array(dat['Z'])
#ra=np.array(dat['ra'])
#dec=np.array(dat['Dec'])
#hpix=healpixmap(ra,dec)
#print(hpix)
#hu.mollview(hpix,rot=180)
# zr=randz(z,2*len(dat['Z']))
# print(zr)
# spinner.text='Generating Random RA, DEC...'
# spinner.start()
# ra,dec=randang('/Users/rohin/Downloads/mask.dr72safe0.ply',2*len(dat['Z']))
# spinner.stop()
# print(ra,dec)
# plt.hist(z, 10, normed=True, alpha=0.5, label='hist')
# x_grid=np.linspace(min(z),max(z),2*len(dat['Z']))#check  nb file
# kdepdf=kde(z,x_grid,bandwidth=1e-3)#check  nb file#,bandwidth=0.01
# plt.plot(x_grid, kdepdf, color='r', alpha=0.5, lw=3, label='kde')
#plt.show()
#bins=np.arange(0,0.1,0.01)
#print(autocorr(dat,bins,flatdistsq))
#print(crosscorr(dat,dat,bins,flatdistsq))
#print(closedistsq(dat[1],dat[2]))
#mask=readmaskfile('/Users/rohin/Downloads/mask.dr72safe0.ply')
#
#z=dat['Z']
#print(z)
#np.vectorize(Ez)
#ez=Ez(z)
#print ez
#np.vectorize(DC_LCDM)
#s=DC_LCDM(z)
#s=comov(z,'lcdm')
#print s
