#!/usr/bin/en/python
# Filename: cosmology_functions.py

from __future__ import division
import numpy as np

# Constants

G = 43.020 # Mpc/(1e10 M_solar) km^2 / s^2

# Evolution of the Hubble constant

def E(Omega_m0, Omega_l0, z):
    return np.sqrt( Omega_l0 + Omega_m0 * (1+z)**3 )

# Evolution of Omega_m

def Omz(Omega_m0, Omega_l0, z):
    return Omega_m0 * (1+z)**3 / (E(Omega_m0, Omega_l0, z)**2)

# Evolution of Omega_l

def Olz(Omega_m0, Omega_l0, z):
    return Omega_l0 / (E(Omega_m0, Omega_l0, z)**2)

# Evolution of Rho_critical

def Rhocrit_z(Omega_m0, Omega_l0, z):
    Rhocrit_0 = 3.0/(8.0 * np.pi * G) * 1e4 
    return Rhocrit_0 * (E(Omega_m0, Omega_l0, z)**2)

# Real-space Tophat in Fourier space 

def TopHat(k, R):

    # k: wavenumber
    # R: filter size in Mpc/h

    return 3.0/(k*R)**2 * (np.sin(k*R)/(k*R) - np.cos(k*R))

# Integrand that goes in calculation of Sigma(M)

def SigmaIntegrand(k, Pk, R):
    return k**2 * Pk * TopHat(k, R)**2

def linear_growth_factor(Omega_m0, Omega_l0, z):
    if len(np.atleast_1d(z)) == 2:
        z1    = z[0]
        z2    = z[1] # z2 > z1                                                                                            

    if (len(np.atleast_1d(z)) == 1) or (len(np.atleast_1d(z)) > 2):
        z1    = 0.
        z2    = z    # z2 > z1                                                                                           
             
    Omega_lz1 = Omega_l0 / (Omega_l0 + Omega_m0 * (1.+z1)**3)
    Omega_mz1 = 1. - Omega_lz1
    gz1       = (5./2.) * Omega_mz1 / (Omega_mz1**(4./7.) - Omega_lz1 + (1. + Omega_mz1/2.) * (1. + Omega_lz1/70.))
    Omega_lz2 = Omega_l0 / (Omega_l0 + Omega_m0 * (1.+z2)**3)
    Omega_mz2 = 1. - Omega_lz2
    gz2       = (5./2.) * Omega_mz2 / (Omega_mz2**(4./7.) - Omega_lz2 + (1. + Omega_mz2/2.) * (1. + Omega_lz2/70.))
    return (gz2 / (1.+z2)) / (gz1 / (1+z1))
