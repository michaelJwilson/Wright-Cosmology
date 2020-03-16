#!/usr/bin/env python
  
import sys
import argparse
import numpy as np

from   math      import *


parser = argparse.ArgumentParser(description='Ned Wright Cosmology Calcalculator.')

parser.add_argument('--z',  metavar='z',  type=np.float, nargs='?', default=3.0,                              help='redshift to calculate')
parser.add_argument('--H0', metavar='H0', type=np.float, nargs='?', default=69.6, help='Hubble constant today.')
parser.add_argument('--Om', metavar='Om', type=np.float, nargs='?', default=0.286, help='Om.')

##  WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528.
parser.add_argument('--Ov', metavar='-Ov', type=np.float, nargs='?', default=1.0 - 0.286 - 0.4165/(69.6*69.6), help='Ov.')

parser.add_argument('--verbose', action='store_true', help='Be verbose.', default=1)

args = parser.parse_args()

# Unpack
z       = args.z
H0      = args.H0
WM      = args.Om
WV      = args.Ov
verbose = args.verbose

# initialize constants
WR  = 0.         # Omega(radiation)
WK  = 0.         # Omega curvaturve = 1-Omega(total)
c   = 299792.458 # velocity of light in km/sec
Tyr = 977.8      # coefficent for converting 1/H into Gyr
DTT = 0.5        # time from z to now in units of 1/H0
DTT_Gyr = 0.0    # value of DTT in Gyr
age = 0.5        # age of Universe in units of 1/H0
age_Gyr = 0.0    # value of age in Gyr
zage = 0.1       # age of Universe at redshift z in units of 1/H0
zage_Gyr = 0.0   # value of zage in Gyr
DCMR = 0.0       # comoving radial distance in units of c/H0
DCMR_Mpc = 0.0 
DCMR_Gyr = 0.0
DA = 0.0         # angular size distance
DA_Mpc = 0.0
DA_Gyr = 0.0
kpc_DA = 0.0
DL = 0.0         # luminosity distance
DL_Mpc = 0.0
DL_Gyr = 0.0     # DL in units of billions of light years
V_Gpc = 0.0
a = 1.0          # 1/(1+z), the scale factor of the Universe
az = 0.5         # 1/(1+z(object))

h = H0/100.
WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
WK = 1-WM-WR-WV
az = 1.0/(1+1.0*z)
age = 0.
n=1000         # number of points in integrals

for i in range(n):
    a = az*(i+0.5)/n
    adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
    age = age + 1./adot

zage = az*age/n
zage_Gyr = (Tyr/H0)*zage
DTT = 0.0
DCMR = 0.0

# do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
for i in range(n):
    a = az+(1-az)*(i+0.5)/n
    adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
    DTT = DTT + 1./adot
    DCMR = DCMR + 1./(a*adot)

DTT = (1.-az)*DTT/n
DCMR = (1.-az)*DCMR/n
age = DTT+zage
age_Gyr = age*(Tyr/H0)
DTT_Gyr = (Tyr/H0)*DTT
DCMR_Gyr = (Tyr/H0)*DCMR
DCMR_Mpc = (c/H0)*DCMR

# tangential comoving distance
ratio = 1.00
x = sqrt(abs(WK))*DCMR
if x > 0.1:
  if WK > 0:
    ratio =  0.5*(exp(x)-exp(-x))/x 
  else:
    ratio = sin(x)/x
else:
  y = x*x
  if WK < 0:
      y = -y
  ratio = 1. + y/6. + y*y/120.
DCMT = ratio*DCMR
DA = az*DCMT
DA_Mpc = (c/H0)*DA
kpc_DA = DA_Mpc/206.264806
DA_Gyr = (Tyr/H0)*DA
DL = DA/(az*az)
DL_Mpc = (c/H0)*DL
DL_Gyr = (Tyr/H0)*DL

# comoving volume computation
ratio = 1.00
x = sqrt(abs(WK))*DCMR
if x > 0.1:
  if WK > 0:
    ratio = (0.125*(exp(2.*x)-exp(-2.*x))-x/2.)/(x*x*x/3.)
  else:
    ratio = (x/2. - sin(2.*x)/4.)/(x*x*x/3.)
else:
  y = x*x
  if WK < 0:
    y = -y
  ratio = 1. + y/5. + (2./105.)*y*y
VCM = ratio*DCMR*DCMR*DCMR/3.
V_Gpc = 4.*pi*((0.001*c/H0)**3)*VCM

print('\n\n--------------  Ned Wright Cosmology calculator  -------------\n\n')
print('For H0 = ' + '%1.3f' % H0 + ', Om = ' + '%1.3f' % WM + ', Ov = {:1.3f}, z= {:1.3f}.\n'.format(WV, z))
print('It is now {:.3f} Gyr since the Big Bang.\n'.format(age_Gyr))
print('The age at redshift z was {:.3f} Gyr.\n'.format(zage_Gyr))
print('The comoving radial distance is {:.3f} Mpc/h.\n'.format(h * DCMR_Mpc))
print('The comoving volume within redshift z is {:.3f} (Gpc/h)^3.\n'.format(h*h*h*V_Gpc))
print('The angular diameter distance is {:.3f} Mpc/h.\n'.format(h * DA_Mpc))
print('This gives a scale of {:.3f} kpc/h/".\n'.format(h * kpc_DA))
print('The luminosity distance is {:.3f} Mpc/h.\n'.format(h * DL_Mpc))
print('The distance modulus (m-M) is {:.3f}.'.format(5.*np.log10(DL_Mpc*1.e6)-5.))
print('\n\nDone.\n\n')
