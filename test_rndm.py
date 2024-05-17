import numpy as np
import matplotlib.pyplot as plt
from RC_vector_mesons import photon_spectrum as f
from scipy.integrate import quad
from math import log as ln
import os

PI = np.pi
ALPHA = 1. / 137.036
ELECTRON_MASS = 0.000511  # in GeV
MUON_MASS = 0.10565  # in GeV
PROTON_MASS = 0.93827  # in GeV
JPSI_MASS = 3.096916  # M of J/psi

mJ = JPSI_MASS
me = ELECTRON_MASS

# Define min and max values
min_Eg, max_Eg = 0.001, 1.0

# Define the function g(x) normalized to 1 from min_Eg to max_Eg
def g(x, min_Eg, max_Eg):
    return 1 / (x * np.log(max_Eg / min_Eg))

# Create an array of x values
x = np.linspace(min_Eg, max_Eg, 10000)

# Calculate the integral of f(x) numerically for normalization
# integral_f is the probabiity to radiate photon in the range (min_Eg,max_Eg)
# Set min_Eg = 0.001 and max-Eg = 1.0 for our case

integral_f, _ = quad(f, min_Eg, max_Eg, args=(mJ, me))

# Normalize f(x)
def f_normalized(x, mJ, me):
    return f(x, mJ, me) / integral_f

# Calculate the ratio f(x) / g(x)
ratio = f_normalized(x, mJ, me) / g(x, min_Eg, max_Eg)

# Find the maximum of the ratio
max_ratio  = np.max(ratio)

# Define the new function G(x) that is > f(x)_normalized in (min_Eg, max_Eg)
def G(x, min_Eg, max_Eg):
    return g(x, min_Eg, max_Eg) * max_ratio

# Save pdf and png
def save(filename):
    file_f = 'Figures/'+filename
    os.system('rm -f '+file_f+'.pdf')
    os.system('rm -f '+file_f+'.png')
    plt.savefig(file_f+'.pdf')
    plt.savefig(file_f+'.png')
    
# Plot the functions for the first range
xx=x[0:1000]
plt.figure(figsize=(12, 6),num='Comparison of f(x), g(x), and G(x)')
plt.plot(xx, f_normalized(xx, mJ, me), label='f(x) = dNdk0 (normalized)', color='blue')
plt.plot(xx, g(xx, min_Eg, max_Eg), label='g(x) = 1 / (x ln(max/min))', color='red')
plt.plot(xx, G(xx, min_Eg, max_Eg), label='G(x) = g(x) * max(f(x)/g(x))', color='green')
plt.xlabel('x')
plt.ylabel('Density')
plt.title('Comparison of f(x), g(x), and G(x) (Range: 0.001 to 1)')
plt.yscale('log')
plt.legend()
save('f_g_G')
plt.show(block=False)

# Generate random samples according to normalized g(x)
def generate_g_samples(n_samples, min_Eg, max_Eg):
    R = np.random.uniform(0, 1, n_samples)  # Generate R uniformly in [0, 1]
    samples = min_Eg * (max_Eg / min_Eg) ** R
    return samples

# Adjust the samples using rejection sampling for f(x)
def rejection_sampling(samples, min_Eg, max_Eg):
    accepted_samples = []
    for x in samples:
        y = np.random.uniform(0, G(x, min_Eg, max_Eg))
        if y < f_normalized(x, mJ, me):  # Accept the sample if y < f(x)
            accepted_samples.append(x)
    return np.array(accepted_samples)

def rejection_sampling2(samples, min_Eg, max_Eg):
    accepted_samples = []
    for x in samples:
        y = np.random.uniform(0, G2(x, min_Eg, max_Eg))
        if y < f_normalized(x, mJ, me):  # Accept the sample if y < f(x)
            accepted_samples.append(x)
    return np.array(accepted_samples)

# Parameters
n_samples = 100000

# Step 1: Generate samples according to normalized g(x)
g_samples = generate_g_samples(n_samples, min_Eg, max_Eg)

# Step 2: Adjust samples to follow f(x) using rejection sampling
f_samples = rejection_sampling(g_samples, min_Eg, max_Eg)

# Plot the histogram of the f_samples
plt.figure(figsize=(12, 6),num='F-normalized in range 0.001,1.0')
#plt.hist(f_samples, bins=500, density=True, alpha=0.6, color='blue', label='f_samples (0.001 to 1)',range=(0.001,1.0))
plt.hist(f_samples, bins=1000,                alpha=0.6, color='blue', label='f_samples (0.001 to 1)',range=(0.001,1.0))
step = (1-0.001)/1000
norm = len(f_samples)*step
plt.plot(x, norm*f_normalized(x, mJ, me), 'magenta', label='f(x)')
plt.xlabel('x')
plt.yscale('log')
plt.ylabel('Density')
plt.title('Rejection Sampling of f(x) (Range: 0.001 to 1)')
plt.legend()
save('dNdk0_rndm')
plt.show(block=False)


