from correlcalc import *
bins = np.arange(0.002,0.062,0.002)
corrdr7flcdmls=tpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,randfile='/usr3/vstr/yrohin/Downloads/random-DR7-Full.ascii',estimator='ls',weights='eq')

