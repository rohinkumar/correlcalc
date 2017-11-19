from correlcalc import *
bins = np.arange(0.002,0.062,0.002)
corrdr7flcls=tpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,randfile='/usr3/vstr/yrohin/Downloads/random-DR7-Full.ascii',estimator='ls',cosmology='lc',weights='eq')
corrdr7olcls=tpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,randfile='/usr3/vstr/yrohin/Downloads/random-DR7-Full.ascii',estimator='ls',cosmology='lc',weights='eq',geometry='open')
corrdr7flclsw=tpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,randfile='/usr3/vstr/yrohin/Downloads/random-DR7-Full.ascii',estimator='ls',cosmology='lc',weights=True)
corrdr7flolsw=tpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,randfile='/usr3/vstr/yrohin/Downloads/random-DR7-Full.ascii',estimator='ls',cosmology='lc',weights=True,geometry='open')
corrdr7flcdmls=tpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,randfile='/usr3/vstr/yrohin/Downloads/random-DR7-Full.ascii',estimator='ls',cosmology='lcdm',weights='eq')
corrdr7flcls=tpcf('/usr3/vstr/yrohin/Downloads/DR7-Full.ascii',bins,randfile='/usr3/vstr/yrohin/Downloads/random-DR7-Full.ascii',estimator='ls',cosmology='lcdm',weights=True)
