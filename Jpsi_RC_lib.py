
import numpy as np
import matplotlib.pyplot as plt
import math
import os
from scipy.special import spence
from math import log as ln

# Constants

PI = math.pi
alpha = 1. / 137.036
me  = 0.000511  # electron mass in GeV
mmu = 0.10565   # muon mass
M  = 0.93827    # proton   mass in GeV
Mj = 3.096916   # M of J/psi

### F. Ehlotzky, Nuovo Cimento Vol LV A, N 1, (1968)
### M. Vanderghagen Phys. Rev. D97, 076012 (2018) 

def Li2(x):
    return spence(1 - x)

def dNdk0(k0, Mee):  # Formula (23) from CERN paper
    m = me      # electron mass
    E = Mee / 2
    x = k0 / E
    return (2 * alpha / PI) * (2 * np.log(2 * E / m * np.sqrt(1 - x)) - 1) * (1 / k0) * (1 - x + x * x / 2)

def dNdk0_muon(k0, Mee):  # Formula (23) from CERN paper
    m = mmu      # muon mass
    E = Mee / 2
    x = k0 / E
    return (2 * alpha / PI) * (2 * np.log(2 * E / m * np.sqrt(1 - x)) - 1) * (1 / k0) * (1 - x + x * x / 2)

def dNdEg(x):
    #    return -alpha/PI*(1+np.log(m*m/(3.1*3.1)))*2/x*(delta_RC65exp(x,3.1)+1)
    return -alpha / PI * (1 + np.log(me * me / (Mj * Mj))) * 2 / x

def dNdEg1(x, n, eps):
    return np.power(eps, n) / np.power(x, n) * dNdEg(eps)


def dNdk01(x, n, eps):
    return np.power(eps, n) / np.power(x, n) * dNdk0(eps, Mj)


def R(rmax, E):    # Eq.25
    # rmax = kmax/E where kmax is the maximum of the photon energy
    # E       = mV / 2, where mV is the vector meson mass
    m  = me # electron mass
    return (rmax - rmax * rmax / 4 - 3 / 4) * (np.log(2 * E / m * np.sqrt(1 - rmax)) - 0.5)


def delta_H(eps, kmax, mV):  # Eq.24, Hard radiation correction 
    # eps   = minimum photon energy for hard production.
    # kmax = maximum photon energy for hard production.
    m = me  # electron mass
    E = mV / 2
    rmax = kmax / E
    R = (rmax - rmax * rmax / 4 - 3 / 4) * (np.log(2 * E / m * np.sqrt(1. - rmax)) - 0.5)
    return -4 * alpha / PI * ((np.log(2 * E / m) - 0.5) * (np.log(eps / kmax) + 3 / 4) + 0.5 * Li2(rmax) + R)


def delta_R(eps, mV):   # Eq.22 Soft photon RC
    # eps   = minimum photon energy for hard production.
    # mV    = vector meson mass
    m = me  # electron mass
    E = mV / 2
    return -4 * alpha / PI * ((np.log(2 * E / m) - 0.5) * (np.log(E / eps) - 13 / 12) + 17 / 72 - PI * PI / 12)


def delta_HR(kmax, mV):  # = Delta_R + delta_H does not depend on eps
    m = me  # electron mass
    E = mV / 2
    eps = m * 0.01
    return delta_H(eps, kmax, mV) + delta_R(eps, mV)

def delta_H_muons(eps, kmax, mV):  # Eq.24, Hard radiation correction 
    # eps   = minimum photon energy for hard production.
    # kmax = maximum photon energy for hard production.

    m = mmu   # muon mass

    E = mV / 2
    rmax = kmax / E
    R = (rmax - rmax * rmax /4. - 3./4.) * (np.log(2. * E / m * np.sqrt(1. - rmax)) - 0.5)
    return -4 * alpha / PI * ((np.log(2 * E / m) - 0.5) * (np.log(eps / kmax) + 3 / 4) + 0.5 * Li2(rmax) + R)


def delta_R_muons(eps, mV):   # Eq.22 Soft photon RC
    # eps   = minimum photon energy for hard production.
    # mV    = vector meson mass
    E = mV / 2
    return -4 * alpha / PI * ((np.log(2 * E / mmu) - 0.5) * (np.log(E / eps) - 13 / 12) + 17 / 72 - PI * PI / 12 +
                              + 5./18.-1./3.*np.log(2*E/me))

def delta_VM(kmax, mV):  # Eq.26 F. Ehlotzky, Nuovo Cimento Vol LV A, N 1, (1968)
    #     delta_VM(0.80805,Mj) = -2.3878017445380593e-06 = zero. We can fit upto 0.80805 and get delta=0.
    #     It means that fit upto this value will get the number of events without RC
    m = me  # electron mass 
    E = mV / 2
    rmax = kmax / E
    R = (rmax - rmax * rmax / 4 - 3 / 4) * (np.log(2 * E / m * np.sqrt(1. - rmax)) - 0.5)
    d = -4 * alpha / PI * (
                (1 / 2 - np.log(2 * E / m)) * (np.log(rmax) + 1 / 3) + 1 / 2 * (Li2(rmax) - PI * PI / 6) + 17 / 72 + R)
    return d


def delta_VM_muons(kmax, mV):  # F. Ehlotzky, Nuovo Cimento Vol LV A, N 1, (1968)
    me  = 0.000511  # electron mass in GeV    kmax = Delta
    mmu = 0.10565   # muon mass
    E = mV / 2
    rmax = kmax / E
    R = (rmax - rmax * rmax /4. - 3./4.) * (np.log(2. * E / mmu * np.sqrt(1. - rmax)) - 0.5)
    d = -4.*alpha/PI * (
                (1./2. - np.log(2. * E / mmu)) * (np.log(rmax) + 1./3.) + 1./2. * (Li2(rmax) - PI * PI / 6.) + 17./72. + R
                + 5./18.-1./3.*np.log(2*E/me))
    return d


def delta_Eq24(E, eps, kmax):
    rmax = kmax / E
    return -4 * alpha / PI * ((np.log(2 * E / me) - 0.5) * (np.log(eps / kmax) + 3 / 4) + 0.5 * Li2(rmax) + R(rmax, E))


def delta_Eq22(E, eps):
    return -4 * alpha / PI * ((np.log(2 * E / me) - 0.5) * (np.log(E / eps) - 13 / 12) + 17 / 72 - PI * PI / 12)



###   Phy. Rev. D 97, 076012 (2018)

def delta_RC65(Delta, Mee):  # Formula 65
    m = me  # electron mass
    return -(alpha / PI) * (np.log(4 * Delta * Delta / (Mee * Mee)) * (1 + np.log(m * m / (Mee * Mee))) - PI * PI / 3)

def delta_RC65exp(Delta, Mee):  # Formula 65
    return np.exp(delta_RC65(Delta, Mee)) - 1

def delta_RC66(Delta, Mee):  # Formula 66
    m = me  # electron mass
    beta = np.sqrt(1. - 4 * m * m / (Mee * Mee))
    beta2 = beta * beta
    F = 1. - alpha * alpha / 3. * np.power(1. + ((1 + beta2) / (2 * beta)) * np.log((1 - beta) / (1 + beta)), 2)
    exp = -alpha / PI * (np.log(4 * Delta * Delta / (m * m)) + np.log((1 - beta) / (1 + beta))) * (
                1 + ((1 + beta2) / (2 * beta)) * np.log((1 - beta) / (1 + beta)))
    K = 1 - alpha / PI * ((1 - beta) / beta * np.log((1 - beta) / (1 + beta)) + (1 + beta2) / (2 * beta) * (
                4 * Li2(2 * beta / (1 + beta)) - PI * PI))
    return F * np.exp(exp) * K - 1.


def delta_RC64(Delta, Mee):  # Formula 64
    m = me  # electron mass 
    beta = np.sqrt(1. - 4 * m * m / (Mee * Mee))
    beta2 = beta * beta
    return -(alpha / PI) * ((np.log(4 * Delta * Delta / (m * m)) + np.log((1 - beta) / (1 + beta))) * (
                1 + (1 + beta2) / (2 * beta) * np.log((1 - beta) / (1 + beta))) +
                            (1 - beta) / beta * np.log((1 - beta) / (1 + beta)) + (1 + beta2) / (2 * beta) * (
                                        4 * Li2(2 * beta / (1 + beta)) - PI * PI))

def delta_RC64exp(Delta, Mee):  # Formula 65
    return np.exp(delta_RC64(Delta, Mee)) - 1

if __name__ == "__main__":
    print("Jpsi_RC_lib:")
