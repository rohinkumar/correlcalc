from correlcalc import *
bins = np.arange(0.002,0.062,0.002)
corrdr7flcols=tpcf('/Users/rohin/Downloads/DR7-Full.ascii',bins,randfile='/Users/rohin/Downloads/random-DR7-Full.ascii',estimator='ls',cosmology='lc',geometry='open',weights='eq')
