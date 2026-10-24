#!/usr/bin/env python
import os
import sys
#sys.path.insert(0,'/data/AllPulsars/Denis/astroutils/') # for met, mjd conversions.
from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
import argparse
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS
import astropy.units as u
sys.path.insert(0,'/data/AllPulsars/Denis/astroutils/')
import conv_dates_lib as conv

# March 2020 -- David doesn't figure out how to use the Professional packages to do this
from coords import ec2eq

'''
Calculate Fermi LAT photon weights using the simple method of Bruel (A&A submitted, Fall 2018) equation 11.
Also described by equation 2 in Smith et al. (ApJ submitted, Fall 2018).
This script originally written by D.A. Smith, L. Guillemot, and P. Ray , Fall 2016.
Questions, comments to smith_AT_cenbg.in2p3.fr 

# das, 12 June 2017 -- mask phases == -1 which arise if your FT1 is newer than your FT2.

# das, 12 March 2021 -- ported to centos7 / python3. 
    sys.path.insert(0,'/data/presto/python/')
    from presto.parfile import psr_par
Had to diddle, see /data/presto/python/README_das
'''

# Define constant
LogESig = 0.5

# Define gaussian model function, with constant sigma of LogESig
def gauss_loge(x, *p):
    A, mu = p
    return A*np.exp(-(x-mu)**2/(2.*LogESig**2))

# This should really be imported from uw.pulsar.stats
def hmw(phases,weights,m=20,c=4):
    """ Calculate the H statistic (de Jager et al. 1989) and weight each
        sine/cosine with the weights in the argument.  The distribution
        is corrected such that the CLT still applies, i.e., it maintains
        the same calibration as the unweighted version."""

    phases = np.asarray(phases)*(2*np.pi) #phase in radians

    s = (np.asarray([(weights*np.cos(k*phases)).sum() for k in range(1,m+1)]))**2 +\
        (np.asarray([(weights*np.sin(k*phases)).sum() for k in range(1,m+1)]))**2

    return ( (2./(weights**2).sum()) * np.cumsum(s) - c*np.arange(0,m) ).max()
    
def calc_lat_weights(energies, angseps, logeref=4.1, logesig=0.5):
    ''' Compute LAT photon weights
        Based on an approximation of the LAT energy-dependent PSF.
        Bruel (A&A submitted, Fall 2018)
	Smith et al. (ApJ submitted, Fall 2018)
        energies are in MeV
        angseps are in degrees
        logeref is the mean of the log gaussian
        logesig is the sigma of the log gaussian'''
    psfpar0  =  5.445
    psfpar1  =  0.848
    psfpar2  =  0.084
    norm     = 1.
    gam      = 2.
    scalepsf = 3.

    logE     = np.log10(energies)
    sigma    = np.sqrt(psfpar0*psfpar0*np.power(100./energies, 2.*psfpar1) + psfpar2*psfpar2)/scalepsf
    fgeom    = norm*np.power(1+angseps*angseps/2./gam/sigma/sigma, -gam)
    return fgeom * np.exp(-np.power((logE-logeref)/np.sqrt(2.)/logesig,2.))
    
desc = """Compute photon weights for a LAT FT1 file, based on an empirical approximation
to the LAT PSF. Optionally, search over LogERef to maximize the weights, as Philippe Bruel's
SearchPulsation code does.  Writes an output FT1 file with a weight column, if desired.

Coordinates can be read from a .par file, or specified with --ra and --dec options.
--ra and --dec can be in the form --ra=HH:MM:SS.s --dec=DD:MM:SS.d or --ra=ddd.d --dec=ddd.d.

NOTE: The FT1 file must have a pulse phase column if you want to fit for the best LogERef.
This is usually done with the Tempo2 fermi plugin."""

parser=argparse.ArgumentParser(description=desc)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--par", default=None, help="Name of par file to read pulsar position from")
group.add_argument("-r","--ra",default=None,help="RA of the pulsar")
parser.add_argument("-d","--dec",default=None,help="DEC of the pulsar")
parser.add_argument("-e","--logeref",type=float,default=None,help="LogERef to use, default is to find optimal")
parser.add_argument("-n","--neref",type=int,default=26,help="Number of LogERef trials")
parser.add_argument("-emin","--erefmin",type=int,default=2.0,help="Starting LogERef value")
parser.add_argument("-emax","--erefmax",type=int,default=4.5,help="Last LogERef value")
parser.add_argument("-w","--weightcol",default=None,help="Name of weight column in FT1 file")
parser.add_argument("-o","--outfile",default=None,help="Output file name, if desired")
parser.add_argument("--phasecol",default="PULSE_PHASE",help="Name of pulse phase column in FT1 file (default PULSE_PHASE)")
parser.add_argument("--plot",default=True,action="store_true")
parser.add_argument("--tmin",type=float,default=54680.,help="Earliest MJD date used.")
parser.add_argument("--tmax",type=float,default=61700.,help="Latest MJD date used.")
parser.add_argument("FT1",help="Input FT1 file name")

args = parser.parse_args()

tmin = conv.mjd2met(args.tmin)
tmax = conv.mjd2met(args.tmax) 

# Get pulsar position
if args.par is not None:
    print("Reading coords from  par file",args.par)
    sys.path.insert(0,'/data/presto/python/')
    from presto.parfile import psr_par
    par = psr_par(args.par)
    if hasattr(par, 'RAJ'):  psr_coord = SkyCoord(par.RAJ, par.DECJ, unit=(u.hourangle,u.deg),frame=ICRS)
    elif hasattr(par, 'ELONG'): 
       rajd, decjd = ec2eq(par.ELONG, par.ELAT)
       psr_coord = SkyCoord(rajd, decjd, unit="deg", frame=ICRS)
       print("Smith converted ecliptic", par.ELONG, par.ELAT, "to equatorial ", rajd, decjd)
else:
    if ":" in args.ra:
        psr_coord = SkyCoord(args.ra,args.dec,unit=(u.hourangle,u.deg),frame=ICRS)
    else:
        psr_coord = SkyCoord(args.ra,args.dec,unit=u.deg,frame=ICRS)        
    print("Pulsar ",psr_coord)
#

hdulist   = fits.open(args.FT1)
tbdata    = hdulist[1].data
hdulist.close()

ras       = tbdata['RA']
decs      = tbdata['DEC']
times     = tbdata['TIME']
eventcoords = SkyCoord(ras,decs,unit=u.deg,frame=ICRS)
angseps   = psr_coord.separation(eventcoords).degree
energies  = tbdata['ENERGY']


if args.logeref is None:
    phases = tbdata[args.phasecol]
    idx = (phases > 0.) & (times > tmin) & (times < tmax)  #  & (energies > 100.) & (angseps < 2.)
    phases = phases[idx]
    energies = energies[idx]
    angseps = angseps[idx]

    logerefs = np.linspace(args.erefmin, args.erefmax, args.neref)
    h = []
    for le in logerefs:
        weights   = calc_lat_weights(energies, angseps, le, LogESig)
        h.append(hmw(phases,weights))
        print("LogEREF = {0:.2f},  H = {1}".format(le,h[-1]))

    h = np.array(h)
    i = np.argmax(h)
    hmax = h[i]
    logerefmax = logerefs[i]

    # Parameter array is amplitude, mean, sigma of the gaussian
    p0 = [hmax, logerefmax]
    # Fit gaussian to find best LogERef
    coeff, var_matrix = curve_fit(gauss_loge, logerefs, h, p0=p0)
    logerefs_fine = np.linspace(1.5,5.0,100)
    hfit = gauss_loge(logerefs_fine, *coeff)
    logeref_fitmax = coeff[1]
    print("Fitted best LogERef = {0}".format(logeref_fitmax))
    
    weights   = calc_lat_weights(energies, angseps, logeref_fitmax, LogESig)
    h_fitmax = hmw(phases,weights)
    print("H at best fit LogERef = ",h_fitmax)

    if h_fitmax < hmax:
        bestLogERef = logerefmax
        print("Fit did not improve H-test ({1} vs {2}), using best LogERef ({0}) from trials.".format(bestLogERef,h_fitmax,hmax))
    else:
        bestLogERef = logeref_fitmax

    # Now calc weights with chosen logERef
    weights   = calc_lat_weights(energies, angseps, bestLogERef, LogESig)

    if args.plot:
        import matplotlib.pyplot as plt
        f, ax = plt.subplots()
        ax.plot(logerefs,h,'o')
        ax.plot(logerefs_fine,hfit)
        ax.grid(True)
        ax.set_xlabel('LogERef (MeV)')
        ax.set_ylabel('Weighted H-test')
        plt.show()

else:
    weights   = calc_lat_weights(energies, angseps, args.logeref, LogESig)


# Write output file, if requested
if args.outfile is not None:
    hdulist   = fits.open(args.FT1)
    names     = hdulist[1].data.names

    if args.weightcol is None:
        raise RuntimeError('Weight column must be defined, use "-w" option!')

    if args.weightcol in names:
        hdulist[1].data.field(args.weightcol)[:] = weights
    else:
        col   = fits.Column(name=args.weightcol,format='E',array=weights)
        tbhdu = fits.BinTableHDU.from_columns(hdulist[1].data.columns+col,header=hdulist[1].header)
        hdulist[1] = tbhdu

    hdulist.writeto(args.outfile,clobber=True)
    hdulist.close()

print('Done!')
