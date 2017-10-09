__author__ = 'Rohin Kumar Y'
# from fileios import *
# msg = 'Enter Absolute Path to file: '
# f_name = raw_input(msg).strip()
# Ran multiple tests for most methods in the package... need to clean up
# path = file_data_and_path(f_name)
# if path != None:
#       print 'Path:',path
# import tpcf
# from fileios import *
# from comovdist import *
from datprep import *
from tpcf import *
import matplotlib.pyplot as plt
from genrand import *
from runtimeupdate import *
# from metrics import *
# t1 = checkdatfile('./testfile.dat')
# print t1
#
# t2=inputfiles('data')
# print t2
# %matplotlib osx
from datvis import *

# dat=checkinfile('./testfile.dat','data')
# print dat
from metrics import *
from tpcf import *
# dat=readinfile('./testfile.dat','data')
# dat=datprep('./testfile.dat','data','lcdm')
from antpcf import *
# bins=np.arange(0.002,0.082,0.002)
# tpcf('./testfile.dat',bins,randfile='./testfile.dat',method='ls')
# tpcf('./testfile.dat',bins,mask='/Users/rohin/Documents/ipy_notebooks/galsurveystudy/masks/boss_geometry_2011_06_10.ply',cosmology='lc',method='ls')
# bins=np.arange(0,0.201,0.01)
# tpcf(dat,dat,3,bins,flatdistsq,'ls')
# antpcf(dat,dat,bins,flatdistsq,flatdistsq,'ls')
# z=np.array(dat['Z'])
# ra=np.array(dat['ra'])
# dec=np.array(dat['Dec'])
# hpix=healpixmap(ra,dec)
# print(hpix)
# hu.mollview(hpix,rot=180)
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
# plt.show()
# bins=np.arange(0,0.1,0.01)
# print(autocorr(dat,bins,flatdistsq))
# print(crosscorr(dat,dat,bins,flatdistsq))
# print(closedistsq(dat[1],dat[2]))
# mask=readmaskfile('/Users/rohin/Downloads/mask.dr72safe0.ply')
#
# z=dat['Z']
# print(z)
# np.vectorize(Ez)
# ez=Ez(z)
# print ez
# np.vectorize(DC_LCDM)
# s=DC_LCDM(z)
# s=comov(z,'lcdm')
# print s
# For 2pCF for LCDM model
from correlcalc import *
bins = np.arange(0.002, 0.082, 0.002)
correldr72 = tpcf('/Users/rohin/Downloads/DR7-Full.ascii', bins, mask='/Users/rohin/Documents/ipy_notebooks/galsurveystudy/masks/window.dr72safe0.ply', randfact=2)
import matplotlib.pyplot as plt
binMpc = 3000*bins
plt.plot(binMpc[1:], correldr72[0], 'ro-')
plt.show()
plt.yscale('log')
plt.plot(binMpc[4:], correldr72[0][3:], 'ro-')
plt.show()
plt.yscale('log')
plt.xscale('log')
plt.plot(binMpc[3:], correldr72[0][2:], 'ro-')
plt.show()
# For anisotropic 2pCF using delta Z and Z delta theta as in arXiv: 1312.0003
from antpcf import *
bins=np.arange(0.01, 0.201, 0.01)
atpcf('./testfile.dat', bins, randfile='./testfile.dat', estimator='ls', permetric='apzdth', parmetric='apdz')
