from libc.math cimport sin, cos, sqrt, acos, sinh, cosh


def flatdistsq (double[:] x, double[:] y):
    cdef double s1 = x[0]
    cdef double s2 = y[0]
    return s1**2+s2**2-2.0*s1*s2*costh(x,y)


def opendistsq (double[:] x, double[:] y):
    cdef double res = 0.0
    cdef double K = -1.0
    cdef double s1 = sinh(x[0]*0.007)
    cdef double s2 = sinh(y[0]*0.007)
    cdef double c1 = cosh(x[0]*0.007)
    cdef double c2 = cosh(y[0]*0.007)
    cdef double costheta = costh(x,y)
    res = s1*s1+s2*s2-2.0*s1*s2*c1*c2*costheta-K*s1*s1*s2*s2*sqrt(1.0+costheta*costheta)
    return res*20408.1633


def closedistsq (double[:] x, double[:] y):
    cdef double res = 0.0
    cdef double K = 1.0
    cdef double s1 = sin(x[0]*0.007)
    cdef double s2 = sin(y[0]*0.007)
    cdef double c1 = cos(x[0]*0.007)
    cdef double c2 = cos(y[0]*0.007)
    cdef double costheta = costh(x,y)
    res = s1*s1+s2*s2-2.0*s1*s2*c1*c2*costheta-K*s1*s1*s2*s2*sqrt(1.0+costheta*costheta)
    return res*20408.1633


def APzdth(double[:] x, double[:] y):
    return (x[0]+y[0])*0.5*acos(costh(x,y))


def APdz(double[:] x, double[:] y):
    return abs(x[0]-y[0])


def costh(double[:] x, double[:] y):
    return sin(x[2])*sin(y[2])+cos(x[2])*cos(y[2])*cos(x[1]-y[1])


def Ezsq(double zv, double Om, double Ol):
    return Om*(1.0+zv)**3+(1.0-Om-Ol)*(1.0+zv)**2+Ol


def sparsqlcdm(double[:] dat1, double[:] dat2):
    cdef double zv = 0
    cdef double Om = 0.3
    cdef double Ol = 0.7
    zv = (dat1[3] + dat2[3])/2.0
    return (dat1[3]-dat2[3])**2/Ezsq(zv, Om, Ol)


def sparsqlc(double[:] dat1, double[:] dat2):
    cdef double zv = 0
    zv = (dat1[3] + dat2[3])/2.0
    return (dat1[3]-dat2[3])**2/(1.0+zv)**2


def musqlcdmf(double[:] dat1, double[:] dat2):
    return sparsqlcdm(dat1, dat2)/flatdistsq(dat1, dat2)


def musqlcf(double[:] dat1, double[:] dat2):
    return sparsqlc(dat1, dat2)/flatdistsq(dat1, dat2)


def musqlcdmo(double[:] dat1, double[:] dat2):
    return sparsqlcdm(dat1, dat2)/opendistsq(dat1, dat2)


def musqlco(double[:] dat1, double[:] dat2):
    return sparsqlc(dat1, dat2)/opendistsq(dat1, dat2)


def musqlcdmc(double[:] dat1, double[:] dat2):
    return sparsqlcdm(dat1, dat2)/closedistsq(dat1, dat2)


def musqlcc(double[:] dat1, double[:] dat2):
    return sparsqlc(dat1, dat2)/closedistsq(dat1, dat2)


def spersqlcdmf(double[:] x, double[:] y):
    return flatdistsq(x, y) - sparsqlcdm(x, y)


def spersqlcf(double[:] x, double[:] y):
    return flatdistsq(x, y) - sparsqlc(x, y)


def spersqlcdmo(double[:] x, double[:] y):
    return opendistsq(x, y) - sparsqlcdm(x, y)


def spersqlco(double[:] x, double[:] y):
    return opendistsq(x, y) - sparsqlc(x, y)


def spersqlcdmc(double[:] x, double[:] y):
    return closedistsq(x, y) - sparsqlcdm(x, y)


def spersqlcc(double[:] x, double[:] y):
    return closedistsq(x, y) - sparsqlc(x, y)

# def sparcsq(double[:] x, double[:] y):
#     return closedistsq(x, y)*(costh(x,y))**2

# def sparosq(double[:] x, double[:] y):
#     return opendistsq(x, y)*(costh(x,y))**2

# def sparfsq(double[:] x, double[:] y):
#     return flatdistsq(x, y)*(costh(x,y))**2

# def sparsq(double[:] dat1, double[:] dat2, double Om, double Ol):
#    cdef double zv = 0
#    zv = (dat1[0] + dat2[0])/2.0
#    return (dat1[0]-dat2[0])**2/Ezsq(zv, Om, Ol)
