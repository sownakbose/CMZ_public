#!/usr/bin/env python
# Filename: sigma_cmz.py

from __future__ import division
from scipy.special import cbrt, gammainc, erf
from scipy.interpolate import interp1d, UnivariateSpline, splrep, splev
import cosmology_functions as cf
import numpy as np
import matplotlib.pyplot as plt
import warnings
import params

######################################
# LOAD THE PARAMETERS FROM PARAMS.PY #
######################################

Database      = params.Database
name          = params.name
redshift      = params.redshift
MassMin       = params.MassMin
MassMax       = params.MassMax
numGridPoints = params.numGridPoints
N_int         = params.N_int

warnings.filterwarnings("ignore", category = RuntimeWarning, append = 1)

###############################
#   LOAD THE POWER SPECTRA    #
###############################

Database = "./PowerSpectra/"

if (name == "WMAP-1"):
    PS_name = Database + "WMAP1_camb_matterpower_z0_extrapolated.dat"
elif (name == "WMAP-3"):
    PS_name = Database + "WMAP3_camb_matterpower_z0_extrapolated.dat"
elif (name == "WMAP-5"):
    PS_name = Database + "WMAP5_camb_matterpower_z0_extrapolated.dat"
elif (name == "WMAP-7"):
    PS_name = Database + "WMAP7_camb_matterpower_z0_extrapolated.dat"
elif (name == "WMAP-9"):
    PS_name = Database + "WMAP9_camb_matterpower_z0_extrapolated.dat"
elif (name == "Planck"):
    PS_name = Database + "Planck_camb_matterpower_z0_extrapolated.dat"
elif (name == "COCO"):
    PS_name = Database + "COCO_camb_matterpower_z0_extrapolated.dat"
elif (name == "Millennium"):
    PS_name = Database + "Millennium_camb_matterpower_z0_extrapolated.dat"
elif (name == "VVV"):
    PS_name = Database + "Planck_VVV.dat"
elif (name == "PlanckNeutralino"):
    PS_name = Database + "PlanckNeutralino_camb_matterpower_z0_extrapolated.dat"

###############################
#   READ THE POWER SPECTRA    #
###############################

# Cosmology info in the first 5 lines

info_file    = open(PS_name, "r")
cosmo        = info_file.readlines()

OmegaBar     = float(cosmo[0])
OmegaMatter  = float(cosmo[1])
hubble       = float(cosmo[2])
n_spec       = float(cosmo[3])
sigma_8      = float(cosmo[4])

info_file.close()

# Calculate the critical density at this redshift

Rhocrit_z  = cf.Rhocrit_z(OmegaMatter, 1.-OmegaMatter, 0.) # units: 1e10 M_solar/Mpc^3/h^2
Rhocrit_z *= 1e10 # units: M_solar/Mpc^3/h^2

# Print cosmology

print "You have selected the ", name, " power spectrum"
print "OmegaB = %f, OmegaM = %f, OmegaL = %f"%(OmegaBar, OmegaMatter, 1.-OmegaMatter)
print "h [100 km/s/Mpc] = %f"%(hubble)
print "Spectral Index = %f"%(n_spec)
print "Sigma (8 Mpc/h) = %f"%(sigma_8)
print "z = %f"%(redshift)

# Read in P(k) and k

PowerSpectrum    = np.genfromtxt(PS_name, skip_header = 5)
Pk_file, k_file  = PowerSpectrum[:,0], PowerSpectrum[:,1]

###############################
#     CALCULATE SIGMA (M)     #
###############################

# Interpolate the power spectrum 

Pk_interp   = interp1d(k_file, Pk_file)

# Omega_m at this redshift

#Omz         = cf.Omz(OmegaMatter, 1.-OmegaMatter, redshift)
Omz          = OmegaMatter


# Mean matter density at this redshift

Rhomean_z   = Rhocrit_z * Omz # M_solar/Mpc/h^2

dlogm       = (np.log10(MassMax) - np.log10(MassMin)) / (numGridPoints-1)
logM0       = np.log10(MassMin) + np.arange(numGridPoints)*dlogm + 0.5*dlogm

filter_Mass  = 10**logM0
R            = cbrt(filter_Mass / (4/3 * np.pi * Rhomean_z)) # Mpc/h

k_min    = k_file.min() * 1.10
k_max    = k_file.max() * 0.90

print k_min, k_max

# Tophat

def TopHat(k, r):
    return 3.0/(k*r)**2 * (np.sin(k*r)/(k*r) - np.cos(k*r))

def SigmaIntegrand(k, r):
    return k**2 * Pk_interp(k) * TopHat(k,r)**2

# Integration function

def integratePk(kmin, kmax, r):

    log_k_min   = np.log10(kmin) 
    log_k_max   = np.log10(kmax) 

    # Size of intervals
    
    dlogk       = (log_k_max - log_k_min)/(N_int - 1)
    tot_sum     = 0. 

    for ii in range(N_int):
        logk     = log_k_min + dlogk*ii
        sum_rect = SigmaIntegrand(10**logk, r) * 10**logk * np.log(10)
        tot_sum  = tot_sum + sum_rect

    # Add contributions from the ends of the integration interval

    log_k_min    = log_k_min - dlogk
    sum_rect_min = SigmaIntegrand(10**log_k_min, r) * 10**log_k_min * np.log(10)
    log_k_max    = log_k_max + dlogk
    sum_rect_max = SigmaIntegrand(10**log_k_max, r) * 10**log_k_max * np.log(10)
    
    sigma_sq     = (tot_sum + 0.5*sum_rect_min + 0.5*sum_rect_max) * dlogk
    sigma_sq    /= (2*np.pi**2)
    return sigma_sq

Sigma_Sq     = np.zeros(len(filter_Mass))

# Sigma^2(M)
Sigma_Sq     = integratePk(k_min, k_max, R)

# Sigma (M)
Sigma        = np.sqrt(Sigma_Sq)

sig_interp   = interp1d(filter_Mass, Sigma)
MassIn8Mpc   = 4/3 * np.pi * 8**3 * Rhomean_z

sig_8        = sig_interp(MassIn8Mpc)

print("Sigma_8 (now): %f"%(sig_8))

normalise    = sig_8 / sigma_8
Sigma       /= normalise

print("Normalisation factor: %f"%(normalise))

# dlog Sigma^2(M) / dlog M

Sigma_Sq     = Sigma**2 
logSigma     = np.log10(Sigma)
logSigma_Sq  = np.log10(Sigma_Sq)
logMass      = np.log10(filter_Mass)

derivSigma   = np.diff(logSigma_Sq) / np.diff(logMass)

###############################
#  CALCULATE c(M) RELATION    #
###############################

A         = 650. / 200
f         = 0.02
delta_sc  = 1.686
alpha_rho = 0.18

c_array    = 10**(np.linspace(-1,2,1000))
delta_sc_0 = delta_sc / cf.linear_growth_factor(OmegaMatter, 1.-OmegaMatter, redshift)
c_ein      = np.zeros(len(filter_Mass))
c_nfw      = np.zeros(len(filter_Mass))

OmegaL      = 1.-OmegaMatter
sig2_interp = splrep(logM0-10., Sigma_Sq, k=1)

for jj in range(numGridPoints):

    M2          = gammainc(3.0 / alpha_rho, 2.0 / alpha_rho) / gammainc(3.0 / alpha_rho, 2.0 * c_array**alpha_rho / alpha_rho)
    rho_2       = 200. * c_array**3 * M2
    rhoc        = rho_2 / (200. * A)
    z2          = (1. / OmegaMatter *(rhoc* (OmegaMatter*(1+redshift)**3 + OmegaL) - OmegaL))**0.3333 - 1.
    delta_sc_z2 = delta_sc / cf.linear_growth_factor(OmegaMatter, OmegaL, z2)

    sig2fM      = splev(logM0[jj] -10. + np.log10(f), sig2_interp)
    sig2M       = Sigma_Sq[jj]

    sig2Min     = splev(np.log10(M2), sig2_interp)

    arg         = A*rhoc/c_array**3 - (1.-erf( (delta_sc_z2-delta_sc_0) / np.sqrt(2.*(sig2fM-sig2M)) ))
    mask        = np.isinf(arg) | np.isnan(arg)
    arg         = arg[mask == False]
    c_array     = c_array[mask==False]
    c_ein[jj]   = np.interp(0.0, arg, c_array)

    # NFW

    M2          = (np.log(2.)-0.5) / (np.log(1.+c_array)-c_array/(1.+c_array))

    rho_2       = 200. * c_array**3 * M2
    rhoc        = rho_2 / (200. * A)
    z2          = (1. / OmegaMatter * (rhoc * (OmegaMatter*(1+redshift)**3 + OmegaL) - OmegaL))**0.33333 - 1.
    delta_sc_z2 = delta_sc / cf.linear_growth_factor(OmegaMatter, OmegaL, z2)
    arg         = A*rhoc/(c_array**3) - (1.-erf((delta_sc_z2-delta_sc_0) / (np.sqrt(2.*(sig2fM-sig2M)))))
    mask        = np.isnan(arg) | np.isinf(arg)
    arg         = arg[mask==False]
    c_array     = c_array[mask==False]
    c_nfw[jj]   = np.interp(0.0, arg, c_array)

###########################
# PLANCK FITTING FUNCTION #
###########################

D_z        = cf.linear_growth_factor(OmegaMatter, OmegaL, redshift)
xi         = 1./(filter_Mass / 1e10)
sig_z_fit  = 22.26 * xi**(0.292) / (1. + 1.53 * xi**0.275 + 3.36 * xi**0.198)
delta_sc_z = delta_sc
nu_z       = delta_sc_z / sig_z_fit
inva       = (1. + redshift)
nu0        = 4.135 - 0.564 * inva - 0.210 * inva**2 + 0.0557 * inva**3 - 0.00348 * inva**4
c0         = 3.395 * inva**(-0.215)
beta       = 0.307 * inva**(0.540)
gamma1     = 0.628 * inva**(-0.047)
gamma2     = 0.317 * inva**(-0.893)
c_nu       = c0 * (nu_z/nu0)**(-gamma1) * (1. + (nu_z/nu0)**(1./beta))**(-beta*(gamma2-gamma1))

####################
# PLOT THE RESULTS #
####################

params = {"text.usetex": True,
          "text.latex.unicode": True,
          "xtick.minor.size": 3.5,
          "xtick.minor.width": 1,
          "ytick.minor.size": 3.5,
          "ytick.minor.width": 1,
          "xtick.major.size": 7.5,
          "xtick.major.width": 2,
          "ytick.major.size": 7.5,
          "ytick.major.width": 2,
          "xtick.labelsize": 25,
          "ytick.labelsize": 25,
          "figure.subplot.left": 0.16,
          "figure.subplot.right": 0.94}

plt.rc("font", family = "serif")
plt.rcParams.update(params)

fig = plt.figure(figsize = (8,8))
ax  = fig.add_subplot(111)
ax.plot(logM0, logSigma, marker = "D", ms = 6, color= "#3366FF", alpha = 0.6, lw = 2.)
ax.set_xlim(-8,17)
ax.set_ylim(-0.5,3.5)
ax.set_xlabel(r"$\log \, M~\left[M_\odot/h\right]$", fontsize = 22)
ax.set_ylabel(r"$\log \, \sigma\left(M\right)$", fontsize = 22)
ax.text(0.2, 0.2, name, fontsize = 25, color = "k", transform = ax.transAxes)
if name == "Planck":
    ax.plot(logM0, np.log10(sig_z_fit), color = "grey", label = "Planck fitting function", lw = 2.)
    ax.legend(loc = "upper right", prop = {"size": 20})
plt.minorticks_on()

fig = plt.figure(figsize=(8,8))
ax  = fig.add_subplot(111)
ax.plot(np.log10(filter_Mass), np.log10(c_ein), lw = 2., marker = "^", alpha = 0.6, ms = 6, color = "#FF6600", label = "Einasto")
ax.plot(np.log10(filter_Mass), np.log10(c_nfw), lw = 2., marker = "o", alpha = 0.6, ms = 6,  color = "#3366FF", label = "NFW")
ax.set_xlabel(r"$\log \, M~\left[M_\odot/h\right]$", fontsize = 22)
ax.set_ylabel(r"$\log \, c_{200}$", fontsize = 22)
if name == "Planck":
    ax.plot(np.log10(filter_Mass), np.log10(c_nu), color = "grey", lw = 2., label = "Planck fitting function")
ax.legend(loc = "upper right", prop = {"size": 20})
ax.set_xlim(-8,17)
ax.set_ylim(0.2,2)
ax.text(0.2, 0.2, name, fontsize = 25, color = "k", transform = ax.transAxes)
ax.text(0.2, 0.1, r"$z=%3.2f$"%(redshift), fontsize = 25, color = "k", transform = ax.transAxes)
plt.minorticks_on()

out = np.column_stack((filter_Mass, c_nfw))
np.savetxt("ludlow_Planck13_z0p00.dat", out, fmt = "%e %f")

plt.show()


