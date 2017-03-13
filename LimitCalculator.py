"""
LimitCalculator.py - 03/03/2017

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

import CEvNS

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
def CalcNevents(gsq=0.0, m_med=1000.0):
    
    integ = lambda x: CEvNS.differentialRate_CEvNS(x,\
         A_Ge, Z_Ge, gsq, m_med)
    return exposure*quad(integ, E_min, E_max, epsrel=1e-4)[0]
    
def CalcLikelihood(gsq, m_med, event_list):

    #Poisson likelihood
    No = len(event_list)
    Ne = CalcNevents(gsq, m_med)
    PL = -Ne + No*np.log(Ne) 
    
    for i in range(No):
        PL += np.log(exposure*CEvNS.differentialRate_CEvNS(event_list[i],\
         A_Ge, Z_Ge, gsq, m_med)/Ne)
    
    return PL
    
def GenerateEvents(N_exp):
    No = np.random.poisson(lam=N_exp)
    events = np.zeros(No)
    
    Rmax = CEvNS.differentialRate_CEvNS(E_min,A_Ge, Z_Ge)
    
    #Get the inverse cumulative distribution
    Evals = np.linspace(E_min, E_max, 100)
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
 
def GetLimit(m_med, event_list):
    L0 = CalcLikelihood(0, m_med, event_list)
    
    #Critical value for one sided test is:
    
    #95%
    X2crit = 2.71
    
    fobj = lambda lg: (2*(L0-CalcLikelihood((10**lg)**2, m_med, event_list))- X2crit)**2
    res = minimize(fobj,-6, method="Nelder-Mead") 
    return 10**res.x
    
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
print "*    LimitCalculator.py    *"
print "****************************"
print " "
print " Calculating Z-prime limits..."

Nobs1 =  CalcNevents()
print " Number of events (SM-only): ", Nobs1
print " Calculate (and plots) likelihood as a function of g..."

m_V = 1000 #MeV

obs_events = GenerateEvents(Nobs1)

Elist = np.linspace(0.1, 1.0, 100)

print CalcNevents(gsq=1.456363040986020291e-05/(3.0*A_Ge), m_med=1000.0)

#Nbins = 9
#pl.figure()
#pl.hist(obs_events, bins=np.linspace(0.1, 1.0, Nbins+1))
#dRdE = np.vectorize(CEvNS.differentialRate_CEvNS)
#Rvals = dRdE(Elist,A_Ge, Z_Ge)
#pl.plot(Elist,1000*Rvals*1.0/Nbins)
#pl.show()

Nvals = 11

mVlist = np.logspace(-3, 7, Nvals)
limlist = mVlist*0.0
for i in tqdm(range(Nvals)):
    limlist[i] = GetLimit_brute(mVlist[i], obs_events)

np.savetxt("ApproxLimits.txt", zip(mVlist, limlist), header="Ricochet Z-prime limit, 0.1keV Threshold, 1000 kg day Exposure")

#Nsamps = 50
#g_list = np.logspace(-9, -1, Nsamps)

#glim = GetLimit(m_V, obs_events)
#print glim

#L0 = CalcLikelihood(0, m_V, obs_events)

#Lvals = g_list*0.0
#for i in tqdm(range(Nsamps)):
#    Lvals[i] = CalcLikelihood(g_list[i]**2.0, m_V, obs_events)

#print Lvals

#print " lnL(g = 0) = ", L0
#print " lnL_max = ", np.max(Lvals)

#pl.figure()
#pl.loglog(g_list, -2*(Lvals-np.max(Lvals)))
#pl.axhline(3.84, linestyle=":")
#pl.show()
