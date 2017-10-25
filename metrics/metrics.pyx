from libc.math cimport sin, cos, sqrt, acos, sinh, cosh


def flatdistsq (double[:] x, double[:] y):
    cdef double s1 = x[0]
    cdef double s2 = y[0]
    return s1**2+s2**2-2.0*s1*s2*mu(x,y)


def opendistsq (double[:] x, double[:] y):
    cdef double res = 0.0
    cdef double K = -1.0
    cdef double s1 = sinh(x[0]*0.007)
    cdef double s2 = sinh(y[0]*0.007)
    cdef double c1 = cosh(x[0]*0.007)
    cdef double c2 = cosh(y[0]*0.007)
    cdef double costheta = mu(x,y)
    res = s1*s1+s2*s2-2.0*s1*s2*c1*c2*costheta-K*s1*s1*s2*s2*sqrt(1.0+costheta*costheta)
    return res*20408.1633


def closedistsq (double[:] x, double[:] y):
    cdef double res = 0.0
    cdef double K = 1.0
    cdef double s1 = sin(x[0]*0.007)
    cdef double s2 = sin(y[0]*0.007)
    cdef double c1 = cos(x[0]*0.007)
    cdef double c2 = cos(y[0]*0.007)
    cdef double costheta = mu(x,y)
    res = s1*s1+s2*s2-2.0*s1*s2*c1*c2*costheta-K*s1*s1*s2*s2*sqrt(1.0+costheta*costheta)
    return res*20408.1633


def APzdth(double[:] x, double[:] y):
    return (x[0]+y[0])*0.5*acos(mu(x,y))


def APdz(double[:] x, double[:] y):
    return abs(x[0]-y[0])


def mu(double[:] x, double[:] y):
    return sin(x[2])*sin(y[2])+cos(x[2])*cos(y[2])*cos(x[1]-y[1])


def sparfsq(double[:] x, double[:] y):
    return flatdistsq(x, y)*(mu(x,y))**2


def sperfsq(double[:] x, double[:] y):
    return flatdistsq(x, y)*(1.0-(mu(x,y))**2)


def sparosq(double[:] x, double[:] y):
    return opendistsq(x, y)*(mu(x,y))**2


def sperosq(double[:] x, double[:] y):
    return opendistsq(x, y)*(1.0-(mu(x,y))**2)


def sparcsq(double[:] x, double[:] y):
    return closedistsq(x, y)*(mu(x,y))**2


def spercsq(double[:] x, double[:] y):
    return closedistsq(x, y)*(1.0-(mu(x,y))**2)
