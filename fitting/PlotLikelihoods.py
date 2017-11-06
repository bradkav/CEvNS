"""
LimitCalculator.py - 13/03/2017

Summary: 
Tool for calculating rough limits on New Physics models 
from measurements ofCoherent Elastic Neutrino Nucleus
Scattering (CEvNS).

Requires numpy, scipy and CEvNS.py.

Author: Bradley J Kavanagh
Please report any problems to: bradkav@gmail.com
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.integrate import quad, cumtrapz
import profile
import CEvNS
import sys
import emcee

from tqdm import * #Progress Bar

#----Plotting parameters----

import matplotlib as mpl
font = { 'size'   : 16, 'family': 'serif', 'serif':'Times New Roman'}
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 2
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 2
mpl.rcParams['ytick.minor.width'] = 1
mpl.rc('font', **font)

import matplotlib.pyplot as pl

#----Notes-----




#Detector parameters

exposure = 1000 #kg-days
A_Ge = 73
Z_Ge = 32
E_min = 0.1 #keV
E_max = 1  #keV


#----Functions-----

#Calculate number of signal events
def CalcNevents(gsq=0.0, m_med=1000.0, tab=False):
    
    integ = lambda x: CEvNS.differentialRate_CEvNS(x,\
         A_Ge, Z_Ge, gsq, m_med, tab)
    return exposure*quad(integ, E_min, E_max, epsrel=1e-4)[0]
    
def CalcLikelihood(gsq, m_med, event_list):

    #Poisson likelihood
    No = len(event_list)
    Ne = CalcNevents(gsq, m_med, tab=True)
    PL = -Ne + No*np.log(Ne) 
    
    for i in range(No):
        PL += np.log(exposure*CEvNS.differentialRate_CEvNS(event_list[i],\
         A_Ge, Z_Ge, gsq, m_med,tab=True)/Ne)
    
    return -PL
    
#Log prior
def lnprior(lQ, lQmin=-15.0, lQmax=10.0):
    if lQmin < lQ < lQmax:
        return 0.0
    return -np.inf
    
#Wrapper for emcee likelihood
def lnprobBG(lQ, m_med, obs_events):
    lp = lnprior(lQ)
    if not np.isfinite(lp):
        return -np.inf
    return lp + CalcLikelihood(10**lQ/(3.0*A_Ge), m_med, obs_events)
    
def lnprob(lQ, m_med, obs_events, L0, lQmin=-15.0, lQmax=10.0):
    #95%
    lp = lnprior(lQ, lQmin, lQmax)
    if not np.isfinite(lp):
        return -np.inf
    X2crit = 2.71
    return lp - (2.0*(CalcLikelihood(10**lQ/(3.0*A_Ge), m_med, obs_events)-L0) - X2crit)**2.0
    
    
def GenerateEvents(N_exp):
    No = np.random.poisson(lam=N_exp)
    print " N_obs = ", No
    events = np.zeros(No)
    
    Rmax = CEvNS.differentialRate_CEvNS(E_min,A_Ge, Z_Ge)
    
    #Get the inverse cumulative distribution
    Evals = np.logspace(np.log10(E_min), np.log10(E_max), 100)
    dRdE = np.vectorize(CEvNS.differentialRate_CEvNS)
    Rvals = dRdE(Evals,A_Ge, Z_Ge)
    Rcum = cumtrapz(Rvals, Evals, initial=0.0)
    Rcum = Rcum/Rcum[-1]
    
    fx = interp1d(Rcum, Evals)
    
    xvals = np.random.rand(No)
    return fx(xvals)

    #Slow rejection method
    """
    for i in range(No):
        x = E_min
        y = Rmax*100
        while (y > CEvNS.differentialRate_CEvNS(x,A_Ge, Z_Ge)):
            x = E_min + np.random.rand(1)*(E_max - E_min)
            y = np.random.rand(1)*Rmax
        events[i] = x
       return events
    """
    
def GetLimit_brute(m_med, event_list):
    Nsamps = 50
    Qmid = np.log10(1e-11 + 1e-12*m_med**2)
    
    Q_list = np.logspace(Qmid-3, Qmid+3, Nsamps)

    #glim = GetLimit(m_V, obs_events)
    #print glim

    Lvals = Q_list*0.0
    for i in range(Nsamps):
        Lvals[i] = CalcLikelihood(Q_list[i]/(3.0*A_Ge), m_med, obs_events)
    
    lnL = -2*(Lvals-np.max(Lvals))
    
    #pl.figure()
    #pl.loglog(Q_list, lnL)
    #pl.axhline(3.84, linestyle=":")
    #pl.show()
    
    return np.min(Q_list[lnL > 3.84])
    

    
#----Main Procedure----
print "****************************"
print "*    PlotLikelihoods.py    *"
print "****************************"
print " "
print " Plotting some example likelihoods..."

# See 1701.07443 - perturbativity limit!
Qmax = 4.0*np.pi*3*A_Ge

Nobs1 =  CalcNevents()
print " Number of events (SM-only): ", Nobs1
mV = 1e6
CEvNS.differentialRate_tabulate(E_min, E_max, A_Ge, Z_Ge, m_med=mV)

#Calculate Qstar
Nplus = CalcNevents(gsq=1.0, tab=True)
Nminus = CalcNevents(gsq=-1.0, tab=True)
Qstar = -3*A_Ge*0.25*(Nplus - Nminus)/(0.5*(Nplus + Nminus) - CalcNevents(gsq=0.0, tab=True))
#print Qstar

obs_events = GenerateEvents(Nobs1+50)
#Qvals = np.logspace(-6.0, np.log10(Qmax), 200)
Qvals = np.logspace(-1.0, np.log10(Qmax), 200)
Ctab = 0.0*Qvals
Nlist = 0.0*Qvals
#Ctrue = 0.0*Evals
for i, Q in enumerate(Qvals):
    Ctab[i] = CalcLikelihood(Q/(3.0*A_Ge), mV, obs_events)
    Nlist[i] = CalcNevents(Q/(3.0*A_Ge), m_med=mV, tab=True)
    
fig, ax1 = pl.subplots()

Ctab = 2.0*(Ctab - np.min(Ctab))

ax1.semilogx(Qvals, Ctab, 'b')
ax1.set_ylim(-1, 1000)
ax1.set_ylabel("-2 lnL", color='b')
ax1.set_xlabel(r"$Q_{Z'}$")
ax1.tick_params('y', colors='b')
ax1.axhline(0, linestyle="--")
ax1.axvline(4.0*np.pi*3*A_Ge, linestyle='-', color='k', linewidth=2.0)
ax1.axvline(Qstar, linestyle='--', color='k', linewidth=1.0)

ax1.fill_between(Qvals, Ctab, where=Ctab > 2.71, alpha=0.5)

ax2 = ax1.twinx()
ax2.semilogx(Qvals, Nlist, 'r')
ax2.set_ylim(0, 1000)
ax2.set_ylabel("# of events", color='r')
ax2.tick_params('y', colors='r')

ax1.set_title("Mediator mass = 1e7 MeV; Typical", fontsize=12)

pl.tight_layout()

#pl.savefig("Typical.pdf")

pl.show()