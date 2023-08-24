#############################################################################
# Create nahmdata function
# Practical imports
import numpy as np
from scipy import integrate
from scipy import special
from scipy import pi

from mpmath import hyp2f1, qfrom, jtheta
from sympy import Float
from numpy import real as RR
from numpy import sqrt

import sys

# create_nahmdata is a function that gets passed to monopole_plotting.py. 
# From the arguments provided, it creates the function nahmdata, which returns
# the three matrices Ti(s) for any given input s in (0,2). 
# As arguments it requires the three parameters that define a V3 3-monopole, k, alpha, and sgn. 
def create_nahmdata(args):
    k, alpha, Delta = args
    k = np.float64(k)
    alpha = np.float64(alpha)
    Delta = np.float64(Delta)

    m = k**2

    if Delta > 0:
        q1 = 1 - m + m**2
        q2 = (m - 2)*(m + 1)
    else:
        q1 = 1 - 16*m + 16*m**2
        q2 = 2*(32*m**2 - 32*m - 1)
    q2 *= 2*m - 1

    bts = np.roots([(4 - 2*alpha)/216, 0.0, -2*q1, -4*q2])
    bts = [bti for bti in bts if np.sign(bti) == -np.sign(alpha) and np.abs(bti) < 12*np.sqrt(q1)]
    bt = bts[0]
    K = special.ellipk(m)
    b = RR(bt * K**2 / 3)
    
    if not -alpha*b > 0:
        print("ct not real")
        sys.exit(2)

    ct = RR(sqrt(-alpha*b**3/27))
    a2 = 12 * (K**2/3)**2 * q1 - b**2/12

    if not a2 > 0:
        print("a not real")
        sys.exit(2)

    a = RR(sqrt(np.abs(a2)))
    g2 = a**2 + b**2/12
    g3 = b*(b**2 - 36*a**2)/216 + ct**2/4
    es = np.roots([4.0, 0.0, -g2, -g3])

    if Delta > 0:
        es = sorted(list(es), reverse=True)
    else:
        es = sorted(list(es), reverse=True, key=lambda z: z.imag)

    e1, e2, e3 = es
    Kp = special.ellipk(1 - m)
    I = 1j
    tau = I*Kp/K

    if Delta > 0:
        q = qfrom(tau=tau)
        added_term = e1
    else:
        q = I*qfrom(tau=tau/2)
        added_term = e2

    def wp(s):
        wp = added_term + (pi/2 * 
                jtheta(3, 0, q)*jtheta(4, 0, q)*jtheta(2, pi*s/2, q)/jtheta(1, pi*s/2, q))**2
        return np.complex128(Float(wp.real) + I*Float(wp.imag))

    def dwp(s, eps=1e-10):
        Ds = [sign*sqrt(4*wp(s)**3 - g2*wp(s) - g3) for sign in [-1, 1]]
        fd = (wp(s+eps) - wp(s))/eps
        return min(Ds, key=lambda d: np.abs(d-fd))
    
    def Gt(s):
        return wp(s).real

    def G(s):
        wp = added_term + (pi/2 * 
                jtheta(3, 0, q)*jtheta(4, 0, q)*jtheta(2, pi*s/2, q)/jtheta(1, pi*s/2, q))**2
        return np.float128(wp.real) + (b + 6*a)/12

    a12 = a
    a31 = -(b + 2*a)/4
    c1 = -(b + 6*a)/12
    c2 = c1 + a12
    c3 = c1 - a31
    cs = [c1, c2, c3]
    f = lambda x: 1/sqrt(np.abs(4*x**3 - g2*x - g3))
    omega1 = 1.
    omega2 = tau
    vs = []

    if Delta > 0:
        for x0 in cs:
            if x0 > e1:
                Kappa = integrate.quad(f, x0, np.inf)
                v = Kappa[0]
            elif x0 > e2:
                Kappa = integrate.quad(f, x0, e1)
                v = Kappa[0]*I + omega1
            elif x0 > e3:
                Kappa = integrate.quad(f, e3, x0)
                v = Kappa[0] + omega2
            else:
                Kappa = integrate.quad(f, -np.inf, x0)
                v = Kappa[0]*I
            v *= np.sign((dwp(v)/ct).imag)
            vs.append(v)
    else:
        for x0 in cs:
            if x0 > e2:
                Kappa = integrate.quad(f, x0, np.inf)
                v = Kappa[0]
            else:
                Kappa = integrate.quad(f, x0, e2.real)
                v = Kappa[0]*I + omega1
            v *= np.sign((dwp(v)/ct).imag)
            vs.append(v)

    eta1 = - pi**2/12 * jtheta(1, 0, q, derivative=3)/jtheta(1, 0, q, derivative=1)
    eta1 = np.float128(Float(eta1.real))

    def sigma(s):
        sig = 2*np.exp(eta1*s**2/2)/pi * jtheta(1, pi*s/2, q)/jtheta(1, 0, q, derivative=1)
        return np.complex128(sig)

    def zeta(s):
        z = pi*s/2
        zeta = eta1*s + pi/2 * jtheta(1, z, q, derivative=1)/jtheta(1, z, q)
        return np.complex128(Float(zeta.real) + I*Float(zeta.imag))

    th0 = [pi/2, pi/2, -pi/2]
    zs = [zeta(vi) for vi in vs]
    csts = list(zip(th0, vs, zs))

    def ets(s):
        etis = []
        for ti0, vi, zi in csts:
            ti = ti0 + 1j * ((s - 1)*(zi - eta1*vi) + 
                    0.5*np.log(np.complex128(jtheta(1, pi*(s - vi)/2, q) *  
                                         jtheta(1, pi*(1 + vi)/2, q) / 
                                         jtheta(1, pi*(s + vi)/2, q) / 
                                         jtheta(1, pi*(1 - vi)/2, q))))
            etis.append(np.complex128(np.exp(1j*ti)))
        return etis

    one = 1+0j
    zero = 0j

    if Delta > 0:
        def nahmdata(t):
            sn, cn, _, _ = special.ellipj(K*t, m)
            F2 = e1 + K**2 * (cn/sn)**2 + (b + 6*a)/12
            f1, f2, f3 = sqrt(F2), sqrt(F2-a12), sqrt(F2+a31)
            et1, et2, et3 = ets(t)
            T1 = f1 * np.array([[zero, zero, zero], [zero, zero, -1/et1], [zero, et1, zero]])
            T2 = f2 * np.array([[zero, zero, et2], [zero, zero, zero], [-1/et2, zero, zero]])
            T3 = f3 * np.array([[zero, -1/et3, zero], [et3, zero, zero], [zero, zero, zero]])
            return T1, T2, T3
    else:
        def nahmdata(t):
            _, cn, _, _ = special.ellipj(2*K*t, m)
            F2 = e2 + K**2 * ((1 + cn)/(1 - cn)) + (b + 6*a)/12
            f1, f2, f3 = sqrt(F2), sqrt(F2-a12), sqrt(F2+a31)
            et1, et2, et3 = ets(t)
            T1 = f1 * np.array([[zero, zero, zero], [zero, zero, -1/et1], [zero, et1, zero]])
            T2 = f2 * np.array([[zero, zero, et2], [zero, zero, zero], [-1/et2, zero, zero]])
            T3 = f3 * np.array([[zero, -1/et3, zero], [et3, zero, zero], [zero, zero, zero]])
            return T1, T2, T3


    return nahmdata
#############################################################################

