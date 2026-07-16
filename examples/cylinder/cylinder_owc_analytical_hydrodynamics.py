import numpy as np
from scipy.optimize import brentq
from scipy.special import jv, j1, iv, i1e, k1e, hankel1; hv = hankel1
from numpy.linalg import solve

def hydrocoeffs(
        a, # cylinder_depth # m
        b, # cylinder_radius # m
        h, # water_depth # m
        frequency, # Hz
        N,
        Nk,
        𝜌_w=1000, # water density # kg/m^2
        g = 9.81, # gravitational acceleration # m/s^2
        _p = 1,  # unit pressure in radiation problem # Pa
        _a = 1, # unit wave amplitude in diffraction problem # m
    ):
    ω = 2*np.pi*frequency
    k = ω**2/g
    # kn, Eq 31-32 (14-15)
    def __kn(n, k, h, eps=1e-8):
        if n==None: # Eq 32 (15)
            ai, bi = 0, 2*max(np.sqrt(k/h), k)
            f = lambda x: k-x*np.tanh(x*h)
            _k = brentq(f, ai, bi)
        elif n==0:
            _k = -1j * __kn(None, k, h)
        else: # Eq 31 (14)
            ai = (n-0.5+eps) * np.pi/h
            bi = (n+0.5-eps) * np.pi/h
            f = lambda  x: k + x*np.tan(x*h)
            _k = brentq(f, ai, bi)
        return _k
    _k = __kn(None, k, h)
    _k0 = __kn(0, k, h)
    _kn = np.zeros(Nk)
    for nk in np.arange(Nk):
        _kn[nk] = __kn(nk+1, k, h)
    # Nn, Eq 30 & 33 (13)
    N0 = 0.5*(1 + np.sin(2*h*_k0)/(2*h*_k0))
    kn2h = 2*h*_kn
    Nn = 0.5*(1 + np.sin(kn2h)/kn2h)
    # Kmn, Eq 93 (~72)
    Kmn = np.zeros((N+1, N+1))
    arg1 = _kn * (h-a)
    arg2 = _kn * b
    fac = (1/(Nn*_kn*h*_kn*b))
    for m in range(N+1):
        for n in range(N+1):
            Kmn[m, n] = np.sum(fac * jv(2*m, arg1) * jv(2*n, arg1) / (i1e(arg2) * k1e(arg2)))
    # L1, L2, Eq 94-95 (73-74)
    L1 = np.zeros((N+1, 1))+0j
    L2 = np.zeros((N+1, 1))+0j
    for m in range(N+1):
        L1[m] = 1 if (m==0) else 0 # Eq 95 (74)
        L2[m] = (-1)**m * 1/np.sqrt(N0) * iv(2*m, _k*(h-a)) # Eq 94 (73)
    # a1, a2, Eq 88 (66)
    a1 = solve(Kmn, L1)
    a2 = solve(Kmn, L2)
    # S mat, Eq 89 (67)
    S11 = np.sum(a1*L1)
    S22 = np.sum(a2*L2)
    S12 = np.sum(a1*L2)
    S21 = np.sum(a2*L1)
    # qr, Eq 61 (~49)
    kb = _k*b
    gamma = np.pi * _k * b * _k * h * j1(kb)
    _Del = S11*S22-S12*S21
    _qr_num = -2*np.pi * k**(-1) * b * (gamma * hv(1, kb) * S11 + 2j*_Del)
    _qr_denom = gamma * hv(1,kb) + 2j*S22
    qr = _qr_num / _qr_denom
    qr *= 1j*ω*_p/(𝜌_w*g)
    qr = np.conj(qr) # Conversion: e^(-iωt) ==> e^(iωt)
    # qs , Eq 76 (~58)
    _qs_num = _a*4*np.pi*1j*_k*b*h*j1(kb)*S21
    _qs_num *= -1j*g/ω * np.sqrt(N0)/np.cosh(_k*h) # Correction from comparing incident wave potential in paper with airy wave potential (e.g. in capytaine)
    _qs_denom = _qr_denom
    qs = _qs_num / _qs_denom
    qs = np.conj(qs) # Conversion: e^(-iωt) ==> e^(iωt)
    return qr/_p, qs/_a
