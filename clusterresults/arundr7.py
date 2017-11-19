from correlcalc import *
bins = np.arange(0.002,0.052,0.002)
acorrdr7flcdmls=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,bins,randfile='/usr3/vstr/yrohin/Downloads/randcat_DR7_2xw.dat',parmetric='sparf',permetric='sperf',estimator='ls',weights='eq')
print("-------------------------------------")
acorrdr7flcdmlsw=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,bins,randfile='/usr3/vstr/yrohin/Downloads/randcat_DR7_2xw.dat',parmetric='sparf',permetric='sperf',estimator='ls',weights=True)
print("-------------------------------------")
print("Calulating for sflat-mu")
binsper = np.arange(0.002,0.052,0.002)
binspar = np.arange(0.952,1.00,0.002)
acorrdr7flcdmsmu=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',binspar,binsper,randfile='/usr3/vstr/yrohin/Downloads/randcat_DR7_2xw.dat',parmetric='mu',permetric='sflat',estimator='ls',weights='eq')
print("-------------------------------------")
acorrdr7flcdmsmu=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',binspar,binsper,randfile='/usr3/vstr/yrohin/Downloads/randcat_DR7_2xw.dat',parmetric='mu',permetric='sflat',estimator='ls',weights=True)
print("-------------------------------------")
bins = np.arange(0.002,0.022,0.002)
acorrdr7ap=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,bins,randfile='/usr3/vstr/yrohin/Downloads/randcat_DR7_2xw.dat',parmetric='apdz',permetric='apzdth',estimator='ls',weights='eq')
print("-------------------------------------")
acorrdr7apw=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,bins,randfile='/usr3/vstr/yrohin/Downloads/randcat_DR7_2xw.dat',parmetric='apdz',permetric='apzdth',estimator='ls',weights=True)