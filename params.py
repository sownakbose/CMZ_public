# Power spectra

# Path to power spectra

Database = "./PowerSpec/"

# Name of power spectrum

#name     = "WMAP-1"
#name     = "WMAP-3"
#name     = "WMAP-5"
#name     = "WMAP-7"
#name     = "WMAP-9"
#name     = "Planck"
#name     = "COCO"
#name     = "Millennium"
#name     = "VVV"
name      = "PlanckNeutralino"

# Choose redshift at which to calculate c(M) relation

redshift  = 0.

# Range of masses to calculate c(M) relation over

MassMin   = 1e-8
MassMax   = 1e16

# Number of mass grid points to generate for c(M) relation

numGridPoints = 5000

# Number of intervals to integrate P(k) (recommended >= 500)

N_int        = 2000
