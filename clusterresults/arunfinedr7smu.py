from correlcalc import *
print("Calulating for sflat-mu")
binspar = np.arange(0.001,0.051,0.001)
binsper = np.arange(0.01,1.01,0.01)
acorrdr7flcdmsmu=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',binspar,binsper,randfile='/usr3/vstr/yrohin/Downloads/random-DR7-Full.ascii',vtype='smu',estimator='ls',weights='eq')
