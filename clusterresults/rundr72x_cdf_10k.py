from correlcalc import *
bins = np.arange(0.002,0.062,0.002)
corrdr7flcdmlsw=tpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,randfile='/usr3/vstr/yrohin/randcat_dr7_cdf_10k.dat',estimator='ls',weights=True)

