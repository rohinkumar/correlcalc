__author__ = 'Rohin Kumar Y'
#Calculate anisotropic 2pCF
import tpcf
from scipy.spatial import distance as dist
#antpcf(dat,datR,bins,parmetric,permetric) returns numpy 2d array DD, RR, DR correl
#poserr(xi,DD) returns (1.0+xi)/np.sqrt(DD)

def antpcf(dat,datR,bins,parmetric,permetric,**kwargs):
    rng=np.array([[min(bins), max(bins)], [min(bins), max(bins)]])
    print "Calculating anisotropic DD..."
    dd2d=np.zeros((20,20))
    ddbt=BallTree(dat,metric='pyfunc',func=permetric)
    for i in tqdm(xrange(len(dat))):
        ind=rrbt.query_radius(dat[i].reshape(1,-1),max(bins))
        for j in ind:
            dist0=dist.cdist([dat[i],],dat[j],parmetric)[0]
            dist1=dist.cdist([dat[i],],dat[j],permetric)[0]
            dd2d+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins))[0]
    print dd2d

    print "Calculating anisotropic RR..."
    rr2d=np.zeros((20,20))
    rrbt=BallTree(datR,metric='pyfunc',func=permetric)
    for i in tqdm(xrange(len(datR))):
        ind=rrbt.query_radius(datR[i].reshape(1,-1),max(bins))
        for j in ind:
            dist0=dist.cdist([datR[i],],datR[j],parmetric)[0]
            dist1=dist.cdist([datR[i],],datR[j],permetric)[0]
            rr2d+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins))[0]
    print rr2d

    print "Calculating anisotropic DR..."
    dr2d=np.zeros((20,20))
    for i in tqdm(xrange(len(dat))):
    ind=rrbt.query_radius(dat[i].reshape(1,-1),max(bins))
    for j in ind:
        dist0=d.cdist([dat[i],],datR[j],parmetric)[0]
        dist1=d.cdist([dat[i],],datR[j],permetric)[0]
        dr2d+=np.histogram2d(dist0, dist1,range=rng,bins=(bins,bins))[0]
    print dr2d

    print "Calculating anisotropic 2pCF with Poisson error"
    Nrd=len(datR)
    N=len(dat)
    f=(1.0*Nrd)/N
    if method='ls':
        correl2d=1.0+f**2*dd2d/rr2d-2.0*f*dr2d/rr2d
    else if method='simple':
        correl2d=f**2*dd2d/rr2d-1.0

    correl2derr = poserr(correl2d,dd2d)

    return correl2d, correl2derr


#def poserr2d(xi,dd2d):
#    return (1.0+xi)/np.sqrt(dd2d)
