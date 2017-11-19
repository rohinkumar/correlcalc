from correlcalc import *
bins = np.arange(0.0005,0.0505,0.0005)
acorrdr122xap=atpcf('/usr3/vstr/yrohin/Downloads/galaxy_DR12v5_CMASS_North.fits',bins,bins,randfile='/usr3/vstr/yrohin/randcat_dr12cmn_2x_pdf10k.dat',vtype='ap',estimator='ls',weights='eq')

