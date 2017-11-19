from correlcalc import *
print("Calulating for sflat-mu")
binspar = np.arange(0.002,0.052,0.002)
binsper = np.arange(0.025,1.025,0.025)
acorrdr122xflcdmsmu=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',binspar,binsper,randfile='/usr3/vstr/yrohin/Downloads/random-DR7-Full.ascii',vtype='smu',estimator='ls',weights='eq')
