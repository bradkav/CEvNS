"""
LimitCalculator_MC.py - 23/03/2017

Summary: 
Tool for calculating limits on New Physics models 
from measurements of Coherent Elastic Neutrino Nucleus
Scattering (CEvNS). Uses MCMC to sample the likelihood.

Requires numpy, scipy and CEvNS.py.
Also requires emcee - http://dan.iel.fm/emcee/current/.

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


#----Detector parameters----

exposure = 1000 #kg-days
A_Ge = 73
Z_Ge = 32
E_min = 0.01 #keV
E_max = 1  #keV

#----Sampling parameters----
ndim, nwalkers = 1, 2
nsteps = 50

#Number of trials per mediator mass
Nsamps = 10

#Chain temperature
Tchain = 4.0

#----Misc. constants----

#Perturbativity limit for Q_Zp
#See 1701.07443
Qmax = 4.0*np.pi*3*A_Ge


#----Functions-----

#Calculate number of signal events
def CalcNevents(gsq=0.0, m_med=1000.0, tab=False):
    
    integ = lambda x: CEvNS.differentialRate_CEvNS(x,\
         A_Ge, Z_Ge, gsq, m_med, tab)
    return exposure*quad(integ, E_min, E_max, epsrel=1e-4)[0]
    
    
#Event generator
def GenerateEvents(N_exp):
    No = np.random.poisson(lam=N_exp)
    #print " N_obs = ", No
    events = np.zeros(No)
        
    #Get the inverse cumulative distribution
    Evals = np.logspace(np.log10(E_min), np.log10(E_max), 100)
    dRdE = np.vectorize(CEvNS.differentialRate_CEvNS)
    Rvals = dRdE(Evals,A_Ge, Z_Ge)
    Rcum = cumtrapz(Rvals, Evals, initial=0.0)
    Rcum = Rcum/Rcum[-1]
    
    fx = interp1d(Rcum, Evals)
    
    xvals = np.random.rand(No)
    return fx(xvals)
      

#Calculate Likelihood (given events in event_list)
def CalcLikelihood(gsq, m_med, event_list):

    #Poisson likelihood
    No = len(event_list)
    Ne = CalcNevents(gsq, m_med, tab=True)
    PL = -Ne + No*np.log(Ne) 
    
    #Event-by-event likelihood
    for i in range(No):
        PL += np.log(exposure*CEvNS.differentialRate_CEvNS(event_list[i],\
         A_Ge, Z_Ge, gsq, m_med,tab=True)/Ne)
    
    return PL*1.0/Tchain
    
#Log prior
def lnprior(lQ, lQmin=-10.0, lQmax=np.log10(Qmax)):
    if lQmin < lQ < lQmax:
        return 0.0
    return -np.inf
    
#Wrappers for emcee likelihood calculations
#Q_Zp positive
def lnprob_positive(lQ, m_med, event_list, lQmin=-15.0, lQmax=np.log10(Qmax)):
    lp = lnprior(lQ, lQmin,lQmax)
    if not np.isfinite(lp):
        return -np.inf
    return lp + CalcLikelihood(10**lQ/(3.0*A_Ge), m_med, event_list)

#Q_Zp negative
def lnprob_negative(lQ, m_med, event_list, lQmin=-15.0, lQmax=np.log10(Qmax)):
    lp = lnprior(lQ, lQmin,lQmax)
    if not np.isfinite(lp):
        return -np.inf
    return lp + CalcLikelihood(-1.0*(10**lQ)/(3.0*A_Ge), m_med, event_list)
    
#Get a (sorted) chain with two columns - values and log-likelihoods
def ProcessChain(sampler):
    chain = np.column_stack((sampler.flatchain[:],sampler.flatlnprobability))
    chain = chain[chain[:,0].argsort()]
    return chain
    
#Take a chain, interpolate it and calculate the upper limit
#Final column should be -2*delta-log-like
def CalcLimit(chain):
    #Interpolate
    like_interp = interp1d(chain[:,0], chain[:,1])
    
    #Resample
    x = np.linspace(np.min(chain[:,0]), np.max(chain[:,0]), 1000)
    
    #Get allowed values (95% CL)
    limvals = x[np.where(like_interp(x) < 2.71)]
    
    #Return upper limit
    if (len(limvals) == 0):
        return 1e-15
    else:
        return 10**np.max(limvals)
    
    
    
#----Main Procedure----
print " "
print "*********************************"
print "*    LimitCalculator_MCMC.py    *"
print "*********************************"
print " "
print " Calculating Z-prime limits (with MCMC)..."


Nobs1 =  CalcNevents()
print " Number of events (SM-only): ",'%.2f' % Nobs1

Nvals = 11 #Number of values of m_Z' to sample
mVlist = np.logspace(-3, 7, Nvals)

limits_pos = np.zeros(Nvals)
limits_neg = np.zeros(Nvals)

for i in range(Nvals):
    mV = mVlist[i]
    
    #Tabulate for a given mediator mass
    CEvNS.differentialRate_tabulate(E_min, E_max, A_Ge, Z_Ge, m_med=mV)
    
    #Calculate Qstar - where the number of events is minimal
    Nplus = CalcNevents(gsq=1.0, tab=True)
    Nminus = CalcNevents(gsq=-1.0, tab=True)
    Qstar = -3*A_Ge*0.25*(Nplus - Nminus)/(0.5*(Nplus + Nminus) - CalcNevents(gsq=0.0, tab=True))
    
    samples_pos = np.zeros(Nsamps)
    samples_neg = np.zeros(Nsamps)
    
    for j in range(Nsamps):
        print " Sample", j+1, "of", Nsamps
        
        #Generate a sample of events
        obs_events = GenerateEvents(Nobs1)
   
        #Set initial points for the random walkers
        #Start high up in the likelihood, near Q*
        pos = [np.log10(Qstar)-1 + 0.5*np.random.randn(ndim) for i in range(nwalkers)]
       
        #Sample for positive values of the couplings
        sampler_positive = emcee.EnsembleSampler(nwalkers, ndim, lnprob_positive, args=(mV, obs_events, -15, np.log10(Qstar)))
        sampler_positive.run_mcmc(pos, nsteps)
        chain = ProcessChain(sampler_positive)  #Get the chain from the sampler 
        chain[:,-1] *= Tchain   #Cool the chain!
        
        #Calculate -2*delta-log-like and then the limit
        chain[:,-1] = -2.0*(chain[:,-1] - np.max(chain[:,-1]))
        samples_pos[j] = CalcLimit(chain)
    
        #Sample for negative values of the couplings
        sampler_negative = emcee.EnsembleSampler(nwalkers, ndim, lnprob_negative, args=(mV, obs_events, -15, np.log10(Qstar)))
        sampler_negative.run_mcmc(pos, nsteps)
        chain = ProcessChain(sampler_negative) #Get the chain from the sampler
        chain[:,-1] *= Tchain   #Cool the chain!
        
        #Calculate -2*delta-log-like and then the limit
        chain[:,-1] = -2.0*(chain[:,-1] - np.max(chain[:,-1]))
        samples_neg[j] = CalcLimit(chain)
        
    limits_pos[i] = np.median(samples_pos)
    limits_neg[i] = np.median(samples_neg)
    print "Median limit (+ve):", np.median(samples_pos)
    print "Median limit (-ve):", np.median(samples_neg)
    print " "

#Save to file
#print mVlist.shape, limits_neg.shape, limits_pos.shape
output = np.column_stack([mVlist, limits_neg, limits_pos])

headtxt = '\n'.join(["Ricochet Z-prime limit on Q_Z', " + str(E_min) + "keV Threshold, 1000 kg day Exposure",\
                     "Columns: m_Z' (MeV), 95% CL Q_Z' < 0, 95% CL Q_Z' > 0"])

np.savetxt("RicochetLimits_Zprime_NoBG_Eth=" + str(E_min) + "keV.txt", output, header=headtxt)