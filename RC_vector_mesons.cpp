
#include <cmath>
#include <iostream>
#include <vector>
#include <gsl/gsl_sf_dilog.h>

// Constants
const double PI = 3.141592653589793;
const double ALPHA = 1.0 / 137.036;
const double ELECTRON_MASS = 0.000511;  // in GeV
const double MUON_MASS = 0.10565;       // in GeV
const double PROTON_MASS = 0.93827;     // in GeV
const double JPSI_MASS = 3.096916;      // M of J/psi

double li2(double x) {
    // Dilogarithm function
    //    return gsl_sf_dilog(1 - x);   in c++ lib : gsl_sf_dilo(x) = Li2(1-x). So the python and c++ codes are diferent
    return gsl_sf_dilog(x);
}

// F. Ehlotzky, Nuovo Cimento Vol LV A, N 1, (1968)

double photon_spectrum(double k, double mV, double particle_mass) {
    // Calculate the differential photon spectrum using Formula (23) from a CERN paper.
    // k  - photon energy
    // mV - Vector meson mass (rho, omega, phi, J/psi)
    // particle_mass - electron mass or muon mass for two decay modes: J/psi --> e+e- or mu+mu-
    double E = mV / 2;
    double x = k / E;
    return (2 * ALPHA / PI) * (2 * std::log(2 * E / particle_mass * std::sqrt(1 - x)) - 1) * (1 / k) * (1 - x + x * x / 2);
}

double soft_photon_correction(double eps, double mV, double particle_mass) {
    // Soft photon radiative correction Eq.22.
    // eps - maximum soft photon energy. I am using eps = 0.001 GeV
    double E = mV / 2;
    double d = -4 * ALPHA / PI * ((std::log(2 * E / particle_mass) - 0.5) * (std::log(E / eps) - 13. / 12) + 17. / 72 - PI * PI / 12);
    if (particle_mass > 0.1) {
        d -= 4 * ALPHA / PI * (5.0 / 18 - 1.0 / 3 * std::log(2 * E / ELECTRON_MASS));
    }
    return d;
}

double hard_photon_correction(double eps, double kmax, double mV, double particle_mass) {
    // Hard photon correction Eq.24.
    // eps  - minimum hard photon energy. I am using eps = 0.001 GeV
    // kmax - maximum hard photon energy to calculate RC. In our case it is around 0.1 GeV
    // mV - Vector meson mass (rho, omega, phi, J/psi)
    // particle_mass - electron mass or muon mass for two decay modes: J/psi --> e+e- or mu+mu-
    double E = mV / 2;
    double rmax = kmax / E;
    double R = (rmax - rmax * rmax / 4 - 3.0 / 4) * (std::log(2 * E / particle_mass * std::sqrt(1 - rmax)) - 0.5);
    return -4 * ALPHA / PI * ((std::log(2 * E / particle_mass) - 0.5) * (std::log(eps / kmax) + 3.0 / 4) + 0.5 * li2(rmax) + R);
}

double delta_vm(double kmax, double mV, double particle_mass) {
    // Soft+Hard radiation correction using the method from F. Ehlotzky, Nuovo Cimento.
    // kmax - maximum hard photon energy to calculate RC. In our case it is around 0.1 GeV
    // mV - Vector meson mass (rho, omega, phi, J/psi)
    // particle_mass - electron mass or muon mass for two decay modes: J/psi --> e+e- or mu+mu-
    double E = mV / 2;
    double rmax = kmax / E;
    double R = (rmax - rmax * rmax / 4 - 3.0 / 4) * (std::log(2 * E / particle_mass * std::sqrt(1 - rmax)) - 0.5);
    double d = -4 * ALPHA / PI * ((1.0 / 2 - std::log(2 * E / particle_mass)) * (std::log(rmax) + 1.0 / 3) + 1.0 / 2 * (li2(rmax) - PI * PI / 6) + 17.0 / 72 + R);
    if (particle_mass > 0.1) {
        d -= 4 * ALPHA / PI * (5.0 / 18 - 1.0 / 3 * std::log(2 * E / ELECTRON_MASS));
    }
    return d;
}

/*
Uncomment lines below and compile (MAC) with the command

$ g++ -o photon_spectrum RC_vector_mesons.cpp -std=c++11 -lgsl -lgslcblas -lm

Then run the program

$ ./photon_spectrum 

Look to the result to make sure the code is correctly working

*/

/*
int main() {
    double mV = JPSI_MASS;
    double eps = 0.001;
    double kmax = 1.0;
    double particle_mass = ELECTRON_MASS;

    std::cout << "\n#####   e+e- Decay Mode #####\n";
    std::cout << "particle Mass = " << particle_mass << "\n";
    std::cout << "mV            = " << mV << "\n";
    std::cout << "eps           = " << eps << "\n";
    std::cout << "kmax          = " << kmax << "\n";
    std::cout << "mV            = " << mV << "\n";
    std::cout << "                                                     Code          Correct Value   Ratio\n";
    std::cout << "photon_spectrum(0.1, mV, particle_mass)              = " << photon_spectrum(0.1, mV, particle_mass) << "    =  0.712194     "<< photon_spectrum(0.1, mV, particle_mass)/0.712194 <<"\n";
    std::cout << "soft_photon_correction(eps, mV, particle_mass)       = " << soft_photon_correction(eps, mV, particle_mass) << "   = -0.472175     "<< soft_photon_correction(eps, mV, particle_mass)/-0.472175<<"\n";
    std::cout << "hard_photon_correction(eps, kmax, mV, particle_mass) = " << hard_photon_correction(eps, kmax, mV, particle_mass) << "    =  0.480878     "<<hard_photon_correction(eps, kmax, mV, particle_mass)/0.480878<<"\n";
    std::cout << "delta_vm(kmax, mV, particle_mass)                    = " << delta_vm(kmax, mV, particle_mass) << "  =  0.00870316   " << delta_vm(kmax, mV, particle_mass)/0.00870316 << "\n";

    particle_mass = MUON_MASS;
    std::cout << "\n#####   mu+mu- Decay Mode #####\n";
    std::cout << "particle Mass = " << particle_mass << "\n";
    std::cout << "mV            = " << mV << "\n";
    std::cout << "eps           = " << eps << "\n";
    std::cout << "kmax          = " << kmax << "\n";
    std::cout << "mV            = " << mV << "\n";
    std::cout << "                                                     Code          Correct Value   Ratio\n";
    std::cout << "photon_spectrum(0.1, mV, particle_mass)              = " << photon_spectrum(0.1, mV, particle_mass) << "    = 0.2477866     "<< photon_spectrum(0.1, mV, particle_mass)/0.2477866 <<"\n";
    std::cout << "soft_photon_correction(eps, mV, particle_mass)       = " << soft_photon_correction(eps, mV, particle_mass) << "     = -0.137599     "<< soft_photon_correction(eps, mV, particle_mass)/-0.137599<<"\n";
    std::cout << "hard_photon_correction(eps, kmax, mV, particle_mass) = " << hard_photon_correction(eps, kmax, mV, particle_mass) << "    = 0.1655178     "<<hard_photon_correction(eps, kmax, mV, particle_mass)/0.1655178<<"\n";
    std::cout << "delta_vm(kmax, mV, particle_mass)                    = " << delta_vm(kmax, mV, particle_mass) << "   = 0.0279182     " << delta_vm(kmax, mV, particle_mass)/ 0.027918 << "\n";

    return 0;
}
*/
