from correlcalc import *
bins = np.arange(0.001,0.021,0.001)
#acorrdr7flcdmls=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,bins,randfile='/usr3/vstr/yrohin/Downloads/randcat_DR7_2xw.dat',vtype='sigpi',estimator='ls',weights='eq')
#print("-------------------------------------")
#acorrdr7flcdmlsw=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,bins,randfile='/usr3/vstr/yrohin/Downloads/randcat_DR7_2xw.dat',vtype='sigpi',estimator='ls',weights=True)
#print("-------------------------------------")
#print("Calulating for sflat-mu")
#binspar = np.arange(0.001,0.02,0.001)
#binsper = np.arange(0.05,1.05,0.05) 
#acorrdr7flcdmsmu=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',binspar,binsper,randfile='/usr3/vstr/yrohin/Downloads/randcat_DR7_2xw.dat',vtype='smu',estimator='ls',weights='eq')
#print("-------------------------------------")
#acorrdr7flcdmsmu=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',binspar,binsper,randfile='/usr3/vstr/yrohin/Downloads/randcat_DR7_2xw.dat',vtype='smu',estimator='ls',weights=True)
#print("-------------------------------------")
#bins = np.arange(0.002,0.022,0.002)
#acorrdr7ap=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,bins,randfile='/usr3/vstr/yrohin/Downloads/randcat_DR7_2xw.dat',vtype='ap',estimator='ls',weights='eq')
#print("-------------------------------------")
#acorrdr7apw=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,bins,randfile='/usr3/vstr/yrohin/Downloads/randcat_DR7_2xw.dat',vtype='ap',estimator='ls',weights=True)

#print("-------------------------------------")
print("Calulating for Rh=ct flat lc sigpi")
#binspar = np.arange(0.001,0.02,0.001)
#binsper = np.arange(0.05,1.05,0.05)
acorrdr7flcsmu=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,bins,randfile='/usr3/vstr/yrohin/Downloads/randcat_DR7_2xw.dat',vtype='sigpi',estimator='ls',weights='eq',cosmology='lc',geometry='flat')
print("-------------------------------------")
acorrdr7flcsmu=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,bins,randfile='/usr3/vstr/yrohin/Downloads/random-DR7-Full.ascii',vtype='sigpi',estimator='ls',weights='eq',cosmology='lc',geometry='flat')

print("-------------------------------------")
print("Calulating for Milne open lc sigpi")
#binspar = np.arange(0.001,0.02,0.001)
#binsper = np.arange(0.05,1.05,0.05)
acorrdr7olcsmu=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,bins,randfile='/usr3/vstr/yrohin/Downloads/randcat_DR7_2xw.dat',vtype='sigpi',estimator='ls',weights='eq',cosmology='lc',geometry='open')
print("-------------------------------------")
acorrdr7olcsmu=atpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,bins,randfile='/usr3/vstr/yrohin/Downloads/random-DR7-Full.ascii',vtype='sigpi',estimator='ls',weights='eq',cosmology='lc',geometry='open')
