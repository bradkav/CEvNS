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

#Sampling parameters
ndim, nwalkers = 1, 2
nsteps = 100

Nvals = 1

#Number of trials per mediator mass
Nsamps = 50

Tchain = 1.0


#Other stuff...

# See 1701.07443 - perturbativity limit!
Qmax = 4.0*np.pi*3*A_Ge

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
    
    return PL*1.0/Tchain
    
#Log prior
#NB: perturbativity limit for couplings!
def lnprior(lQ, lQmin=-10.0, lQmax=np.log10(Qmax)):
    if lQmin < lQ < lQmax:
        return 0.0
    return -np.inf
    
#Wrapper for emcee likelihood
def lnprobBG(lQ, m_med, obs_events, lQmin=-15.0, lQmax=np.log10(Qmax)):
    lp = lnprior(lQ, lQmin,lQmax)
    if not np.isfinite(lp):
        return -np.inf
    return lp + CalcLikelihood(10**lQ/(3.0*A_Ge), m_med, obs_events)

def lnprobBG_negative(lQ, m_med, obs_events, lQmin=-15.0, lQmax=np.log10(Qmax)):
    lp = lnprior(lQ, lQmin,lQmax)
    if not np.isfinite(lp):
        return -np.inf
    return lp + CalcLikelihood(-1.0*(10**lQ)/(3.0*A_Ge), m_med, obs_events)
    
def lnprob(lQ, m_med, obs_events, L0, lQmin=-10.0, lQmax=np.log10(Qmax)):
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
print " "
print "****************************"
print "*    LimitCalculator.py    *"
print "****************************"
print " "
print " Calculating Z-prime limits..."
print " NB: perturbativity limit is Q_Z'=", 12.0*np.pi*A_Ge

"""
Qv = (A_Ge-Z_Ge)
G_FERMI = 1.1664e-5
SQRT2 = np.sqrt(2.0)
qsq = 2*0.9315*A_Ge
def G_V(QvNP, m_med):
    return 1e6*(SQRT2/G_FERMI)*(QvNP/Qv)*1.0/(qsq + m_med*m_med)


m_medlist = np.logspace(0.0, 7.0)

Qlist = np.logspace(-10, np.log10(12.0*np.pi*A_Ge))

#Glist = G_V(12.0*np.pi*A_Ge, m_medlist)
#print m_medlist, Glist

Rlist = 0.0*Qlist
Rlist_tab = 0.0*Qlist
mV = 1e3
CEvNS.differentialRate_tabulate(E_min, E_max, A_Ge, Z_Ge, m_med=mV)
for i, Qi in enumerate(Qlist):
    Rlist[i] = CEvNS.differentialRate_CEvNS(0.5, A_Ge, Z_Ge, gsq=Qi/(3.0*A_Ge), m_med=mV, tab=False)
    Rlist_tab[i] = CEvNS.differentialRate_CEvNS(0.5, A_Ge, Z_Ge, gsq=Qi/(3.0*A_Ge), m_med=mV, tab=True)

pl.figure()
pl.loglog(Qlist, Rlist, '-')
pl.loglog(Qlist, Rlist_tab,'bs')
#pl.axhline(1.0, linestyle="--")
pl.show()
"""



Nobs1 =  CalcNevents()
print " Number of events (SM-only): ", Nobs1


#Initialise emcee


mVlist = np.logspace(-3, 7, Nvals)
limlist = mVlist*0.0
samps_pos = np.zeros(Nsamps)
samps_neg = np.zeros(Nsamps)

for i in range(Nvals):
    mV = mVlist[i]
    mV = 1.0e5   
    print " Calculating for m_med =", mV, "MeV..."
    CEvNS.differentialRate_tabulate(E_min, E_max, A_Ge, Z_Ge, m_med=mV)
    
    #Calculate Qstar
    Nplus = CalcNevents(gsq=1.0, tab=True)
    Nminus = CalcNevents(gsq=-1.0, tab=True)
    Qstar = -3*A_Ge*0.25*(Nplus - Nminus)/(0.5*(Nplus + Nminus) - CalcNevents(gsq=0.0, tab=True))
    
    for j in range(Nsamps):
        print " Sample", j+1, "of", Nsamps
        obs_events = GenerateEvents(Nobs1)
    
    
        pos0 = [np.log10(Qstar) + 1*np.random.randn(ndim) for i in range(nwalkers)]

        #Negative Qzp
        sampler_neg = emcee.EnsembleSampler(nwalkers, ndim, lnprobBG_negative, args=(mV, obs_events))
        for i, result in enumerate(tqdm(sampler_neg.sample(pos0, iterations=nsteps), total=nsteps)):
            pass

        res_neg = np.column_stack((sampler_neg.flatchain[:],sampler_neg.flatlnprobability))
        res_neg = res_neg[res_neg[:,0].argsort()]
        
        res_neg[:,1] *= Tchain
        res_neg[:,1] = -2*(res_neg[:,1] - np.max(res_neg[:,1]))
        #print res_neg
        #pl.figure()
        #pl.plot(res_neg[:,0], res_neg[:,1])
        #pl.plot(res_neg[:,0], res_neg[:,1],'s')
        #pl.show()
        
        like_neg = interp1d(res_neg[:,0], res_neg[:,1])
        x1 = np.linspace(np.min(res_neg[:,0]), np.max(res_neg[:,0]), 1000)
        #pos = [np.log10(3.0*A_Ge) + 0.5*np.random.randn(ndim) for i in range(nwalkers)]
        #print pos
        
        limvals_neg = x1[np.where(like_neg(x1) < 2.71)]
        if (len(limvals_neg) == 0):
            samps_neg[j] = 1e-15
        else:
            samps_neg[j] = 10**np.max(limvals_neg)
    
        pos = [np.log10(Qstar)-2 + 1*np.random.randn(ndim) for i in range(nwalkers)]
       
    
        sampler_L0 = emcee.EnsembleSampler(nwalkers, ndim, lnprobBG, args=(mV, obs_events, -15, np.log10(Qstar)))
        for i, result in enumerate(tqdm(sampler_L0.sample(pos, iterations=nsteps), total=nsteps)):
            pass
    

        res = np.column_stack((sampler_L0.flatchain[:],sampler_L0.flatlnprobability))
        #res = np.sort(res,axis=0)
        res = res[res[:,0].argsort()]
        
        #Cool the chain!
        res[:,1] *= Tchain

        
        #L0 = np.nanmax(res[:,1])
        #res[:,1] = -2*(res[:,1] - L0)
        
        pos = [np.log10(Qstar)+2 + 1*np.random.randn(ndim) for i in range(nwalkers)]
        
        sampler_L1 = emcee.EnsembleSampler(nwalkers, ndim, lnprobBG, args=(mV, obs_events, np.log10(Qstar), np.log10(Qmax)))
        for i, result in enumerate(tqdm(sampler_L1.sample(pos, iterations=nsteps), total=nsteps)):
            pass
    

        res1 = np.column_stack((sampler_L1.flatchain[:],sampler_L1.flatlnprobability))
        res1 = res1[res1[:,0].argsort()]
        
        #Cool the chain!
        res1[:,1] *= Tchain

        L0 = np.maximum(np.nanmax(res[:,1]), np.nanmax(res1[:,1]))
        
        #L0 = np.nanmax(res[:,1])
        res[:,1] = -2*(res[:,1] - L0)
        res1[:,1] = -2*(res1[:,1] - L0)
        
        like = interp1d(res[:,0], res[:,1])
        x = np.linspace(np.min(res[:,0]), np.max(res[:,0]), 1000)
        
        limvals = x[np.where(like(x) < 2.71)]
        if (len(limvals) == 0):
            samps_pos[j] = 1e-15
        else:
            samps_pos[j] = 10**np.max(limvals)
        #samps[j] = 10**lim
        
        """
        pl.figure()
        pl.semilogx(10**res[:,0], res[:,1], 'bs')
        pl.semilogx(10**res[:,0], res[:,1], 'b-')
        
        pl.semilogx(10**res1[:,0], res1[:,1], 'rs')
        pl.semilogx(10**res1[:,0], res1[:,1], 'r-')
        
        pl.axvline(Qstar, color='k', linestyle='--')
        
        pl.ylim(-1,20)
        pl.show()
        """
        
        
        """

        #L0 = np.min(res[:,1])
        x0 = res[np.argmin(res[:,1]),0]
        #print L0
        #res[:,1] -= L0
    
        L_inv = interp1d(res[:,1], res[:,0])
        x = L_inv(2.71)
        print " Min val:", x0
        print " Estimated val:",x
        
        if (x < x0):
            x = x0
 
        
        pos1 = [x + 0.5*np.random.randn(ndim) for i in range(nwalkers)]
        #pos1 = [np.array([np.linspace(-10, 8, nwalkers)[i]]) for i in range(nwalkers)]
        #pos1 = np.append(pos1, np.array([res[np.argmin(res[:,1]),0]]))
        #pos1 = [res[np.argmin(res[:,1]),0] + 0.1*np.random.randn(ndim) for i in range(nwalkers)]
        #print pos1
    
        sampler_L1 = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(mV, obs_events, L0, -15.0, x0))
        for i, result in enumerate(tqdm(sampler_L1.sample(pos1, iterations=nsteps), total=nsteps)):
            pass
    
        res2 = np.column_stack((sampler_L1.flatchain[:],sampler_L1.flatlnprobability))
        #res2 = np.sort(res2, axis=0)
        res2 = res2[res2[:,0].argsort()]
        #print np.max(sampler_L1.flatlnprobability)
        #samp_lower = 10**sampler_L1.flatchain[np.argmax(sampler_L1.flatlnprobability)]
        
        sampler_L2 = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(mV, obs_events, L0, x0, 10.0))
        for i, result in enumerate(tqdm(sampler_L2.sample(pos1, iterations=nsteps), total=nsteps)):
            pass
    
        res2 = np.column_stack((sampler_L2.flatchain[:],sampler_L2.flatlnprobability))
        #res2 = np.sort(res2, axis=0)
        res2 = res2[res2[:,0].argsort()]
        #print np.max(sampler_L1.flatlnprobability)
        #samp_upper = 10**sampler_L2.flatchain[np.argmax(sampler_L2.flatlnprobability)]
        
        if (np.max(sampler_L1.flatlnprobability) > np.max(sampler_L2.flatlnprobability)):
            samps[j] = 10**sampler_L1.flatchain[np.argmax(sampler_L1.flatlnprobability)]
        else:
            samps[j] = 10**sampler_L2.flatchain[np.argmax(sampler_L2.flatlnprobability)]
        
        print "95% CL on Q_Z': ", samps[j]
        
        #50
        if (samps[j] < 100000):
        
            pl.figure()
            pl.plot(res[:,0], 2*res[:,1], 'b+')
            pl.plot(res[:,0], 2*res[:,1], 'b-', label="L0")
            #pl.scatter(sampler_L0.flatchain[:],sampler_L0.flatlnprobability, color='blue')
            #pl.scatter(sampler_L1.flatchain[:],sampler_L1.flatlnprobability, color='red')
            pl.plot(res2[:,0], -res2[:,1], 'r+')
            pl.plot(res2[:,0], -res2[:,1], 'r-', label="L1")
            pl.ylim(0, 100)
            pl.axvline(np.log10(samps[j]), linestyle=':')
            pl.axhline(2.71, linestyle=':')
            pl.legend(loc='best')
            pl.show()
        """
        
    print samps_pos
    print "Median limit (+ve):", np.median(samps_pos)
    print samps_neg
    print "Median limit (-ve):", np.median(samps_neg)
    
    #pl.figure()
    #pl.hist(samps)
    #pl.show()

#np.savetxt("ApproxLimits.txt", zip(mVlist, limlist), header="Ricochet Z-prime limit, 0.1keV Threshold, 1000 kg day Exposure")

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
