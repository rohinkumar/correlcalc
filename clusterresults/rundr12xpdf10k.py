from correlcalc import *
bins = np.arange(0.002,0.062,0.002)
#corrdr12flcdmls=tpcf('/usr3/vstr/yrohin/Downloads/galaxy_DR12v5_CMASS_North.fits',bins,randfile='/usr3/vstr/yrohin/randcat_dr12cmn_2x_pdf10k.dat',estimator='ls',cosmology='lcdm',weights='eq')
print("--------------------------------------------")
corrdr12flcls=tpcf('/usr3/vstr/yrohin/Downloads/galaxy_DR12v5_CMASS_North.fits',bins,randfile='/usr3/vstr/yrohin/randcat_dr12cmn_2x_pdf10k.dat',estimator='ls',cosmology='lcdm',weights=True)
print("--------------------------------------------")
#corrdr12flcls=tpcf('/usr3/vstr/yrohin/Downloads/galaxy_DR12v5_CMASS_North.fits',bins,randfile='/usr3/vstr/yrohin/randcat_dr12cmn_2x_pdf10k.dat',estimator='ls',cosmology='lc',weights='eq')
#print("--------------------------------------------")
#corrdr12olcls=tpcf('/usr3/vstr/yrohin/Downloads/galaxy_DR12v5_CMASS_North.fits',bins,randfile='/usr3/vstr/yrohin/randcat_dr12cmn_2x_pdf10k.dat',estimator='ls',cosmology='lc',weights='eq',geometry='open')
print("--------------------------------------------")
corrdr12flclsw=tpcf('/usr3/vstr/yrohin/Downloads/galaxy_DR12v5_CMASS_North.fits',bins,randfile='/usr3/vstr/yrohin/randcat_dr12cmn_2x_pdf10k.dat',estimator='ls',cosmology='lc',weights=True)
print("--------------------------------------------")
corrdr12flolsw=tpcf('/usr3/vstr/yrohin/Downloads/galaxy_DR12v5_CMASS_North.fits',bins,randfile='/usr3/vstr/yrohin/randcat_dr12cmn_2x_pdf10k.dat',estimator='ls',cosmology='lc',weights=True,geometry='open')
