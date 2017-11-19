from correlcalc import *
bins = np.arange(0.0005,0.0205,0.0005)
acorrdr7ap=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,bins,randfile='/usr3/vstr/yrohin/Downloads/random-DR7-Full.ascii',vtype='ap',estimator='ls',weights='eq')

