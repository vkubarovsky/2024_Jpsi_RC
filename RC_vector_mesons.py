
import numpy as np
from scipy.special import spence

# Constants
PI = np.pi
ALPHA = 1. / 137.036
ELECTRON_MASS = 0.000511  # in GeV
MUON_MASS = 0.10565  # in GeV
PROTON_MASS = 0.93827  # in GeV
JPSI_MASS = 3.096916  # M of J/psi


def li2(x):
    """Dilogarithm function."""
    return spence(1 - x)

### F. Ehlotzky, Nuovo Cimento Vol LV A, N 1, (1968)

def photon_spectrum(k, mV, particle_mass):
    """Calculate the differential photon spectrum using Formula (23) from a CERN paper.
       k  - photon energy
       mV - Vector meson mass (rho,omega,phi, J/psi
       particle_mass - electron mass or muon mass for two decay modes: J/psi-->e+e- or mu+mu-
    """
    E = mV / 2
    x = k / E
    return (2 * ALPHA / PI) * (2 * np.log(2 * E / particle_mass * np.sqrt(1 - x)) - 1) * (1 / k) * (1 - x + x ** 2 / 2)


def soft_photon_correction(eps, mV, particle_mass):
    """Soft photon radiative correction Eq.22.
       eps - maximum soft photon energy. I am using eps = 0.001 GeV
    """
    E = mV / 2
    d  = -4 * ALPHA / PI * ((np.log(2 * E / particle_mass) - 0.5) * (np.log(E / eps) - 13 / 12) + 17 / 72 - PI * PI / 12)
    if particle_mass > 0.1:
        d = d - 4 * ALPHA / PI * (5/18 - 1/3 * np.log(2*E/ELECTRON_MASS))
    return d
def hard_photon_correction(eps, kmax, mV, particle_mass):
    """Hard photon correction Eq.24.
       eps  - minimum hard photon energy. I am using eps = 0.001 GeV
       kmax - maximum hard photon energy to calculate RC. In ourcase it is around 0.1 GeV
       mV - Vector meson mass (rho,omega,phi, J/psi)
       particle_mass - electron mass or muon mass for two decay modes: J/psi-->e+e- or mu+mu-
    """
    E = mV / 2.
    rmax = kmax / E
    R = (rmax - rmax **2/4. - 3./4.) * (np.log(2. * E / particle_mass * np.sqrt(1.-rmax)) - 0.5)
    return -4*ALPHA/PI * ((np.log(2 * E / particle_mass) - 0.5) * (np.log(eps / kmax) + 3. / 4.) + 0.5 * li2(rmax) + R)


def delta_vm(kmax, mV, particle_mass):
    """Soft+Hard radiation correction using the method from F. Ehlotzky, Nuovo Cimento.
       kmax - maximum hrd photon energy to calculate RC. In ourcase it is around 0.1 GeV
       mV - Vector meson mass (rho,omega,phi, J/psi
       particle_mass - electron mass or muon mass for two decay modes: J/psi-->e+e- or mu+mu-
    """
    E = mV / 2
    rmax = kmax / E
    R = (rmax - rmax ** 2 / 4 - 3 / 4) * (np.log(2 * E / particle_mass * np.sqrt(1 - rmax)) - 0.5)
    d = -4 * ALPHA / PI * (
            (1 / 2 - np.log(2 * E / particle_mass)) * (np.log(rmax) + 1 / 3) + 1 / 2 * ( li2(rmax) - PI ** 2 / 6) + 17 / 72 + R)
    if particle_mass > 0.1:
        d = d - 4 * ALPHA / PI * (5/18 - 1/3 * np.log(2*E/ELECTRON_MASS))
    return d


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    import os
    from scipy.special import spence
    from math import log as ln
    # Test of the functions
    mV = JPSI_MASS
    eps = 0.001
    kmax = 1.000
    particle_mass = ELECTRON_MASS
    print()
    print('#####   e+e- Decay Mode #####')
    print('particle Mass = ',particle_mass) 
    print('mV            = ',mV) 
    print('eps           = ',eps) 
    print('kmax          = ',kmax) 
    print('mV            = ',mV) 
    print('                                                       ', 'Value from the code      ',       'Correct Value')
    print('photon_spectrum(0.1, mV, particle_mass)             =  ', photon_spectrum(0.1, mV, particle_mass),       '   =  0.7121948026783391')
    print('soft_photon_correction(eps, mV, particle_mass)      = ', soft_photon_correction(eps, mV, particle_mass),   '   = -0.4721754689879045')
    print('hard_photon_correction(eps,kmax, mV, particle_mass) =  ', hard_photon_correction(eps,kmax, mV, particle_mass), '   =  0.4808786384628666 ')
    print('delta_vm(kmax, mV, particle_mass)                   =  ', delta_vm(kmax, mV, particle_mass),               ' =  0.008703169474962141')


    particle_mass = MUON_MASS
    print()
    print('#####   mu+mu- Decay Mode #####')
    print('particle Mass = ',particle_mass) 
    print('mV            = ',mV) 
    print('eps           = ',eps) 
    print('kmax          = ',kmax) 
    print('mV            = ',mV) 
    print('                                                       ', 'Value from the code      ',       'Correct Value')
    print('photon_spectrum(0.1, mV, particle_mass)             =  ', photon_spectrum(0.1, mV, particle_mass),       '   =  0.24778664264217518')
    print('soft_photon_correction(eps, mV, particle_mass)      = ', soft_photon_correction(eps, mV, particle_mass),   '    = -0.1375996297863609')
    print('hard_photon_correction(eps,kmax, mV, particle_mass) =  ', hard_photon_correction(eps,kmax, mV, particle_mass), '   =  0.16551782544258395')
    print('delta_vm(kmax, mV, particle_mass)                   =  ', delta_vm(kmax, mV, particle_mass),               '  =  0.02791819565622305')

    eg = np.linspace(0.01,1.0,1000)
    y1  = photon_spectrum(eg, mV, ELECTRON_MASS)
    y2  = photon_spectrum(eg, mV, MUON_MASS)
    plt.title(r'$\frac{1}{\Gamma_0}\frac{d\Gamma}{dE_\gamma}$ for electrons and muons', fontsize=22)
    plt.plot(eg,y1,color='red', label=r'$e^+e^-$')
    plt.plot(eg,y2,color='blue',label=r'$\mu^+\mu^-$')
    plt.yscale('log')
    plt.xlabel(r'$E_\gamma$ , GeV')
    plt.ylabel(r'$\frac{1}{\Gamma_0}\frac{d\Gamma}{dE_\gamma}$')
    plt.legend(fontsize=20)
    plt.grid()
    plt.show(block=False)
