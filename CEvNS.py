"""
CEvNS.py - Version 1.6 - 06/11/2017

Summary: 
Code for calculating differential cross section
for Coherent Elastic Neutrino Nucleus Scattering (CEvNS).

Version 1.6: Updated and checked code for Scalar cross section and
             differential rate. Can now calculate scalar xsec
             in terms of product of scalar coupling to quarks
             and neutrinos (assuming universal quark couplings)
Version 1.5: Added separate fluxes for different neutrino species
             You can now specify the flavor of neutrino you're 
             interested in, using the nu_flavor variable in the
             differential rate functions ('e', 'eb', 'mu', 'mub', 'tau', 'taub', or 'all')
Version 1.4: Added neutrino fluxes for COHERENT@SNS

Cross sections mainly taken from arXiv:1604.01025.
See also arXiv:1701.07443.

Magnetic moment - arXiv:1508.07981

Requires numpy and scipy.

Author: Bradley J Kavanagh
Please report any problems to: bradkav@gmail.com
"""

import numpy as np
from scipy.interpolate import interp1d, UnivariateSpline,InterpolatedUnivariateSpline
from scipy.integrate import quad

#----Notes-----
"""
We're assuming that the Z-mass is much larger
than the typical momentum transfer (which is fine
for solar and reactor neutrinos).
"""


#----Constants----
G_FERMI = 1.1664e-5     #Fermi Constant in GeV^-2
SIN2THETAW = 0.2387     #Sine-squared of the weak angle
ALPHA_EM = 0.007297353  #EM Fine structure constant
m_e = 0.5109989461e-3   #Electron mass in GeV  
SQRT2 = np.sqrt(2.0)  


flav = {'e': 0, 'eb': 1, 'mu': 2, 'mub': 3, 'tau': 4, 'taub': 5}

#----Module-scope variables----
neutrino_flux_tot = None
neutrino_flux_list = None

diffrate_A = None
diffrate_B = None

nu_source = None

Enu_min = 0
Enu_max = 0

#----Some functions----

#Helm Form Factor
def HelmFormFactor(E, A):
    #Define conversion factor from amu-->keV
    amu = 931.5*1e3

    #Convert recoil energy to momentum transfer q in keV
    q1 = np.sqrt(2*A*amu*E)

    #Convert q into fm^-1
    q2 = q1*(1e-12/1.97e-7)
    
    #Calculate nuclear parameters
    s = 0.9
    a = 0.52
    c = 1.23*(A**(1.0/3.0)) - 0.60
    R1 = np.sqrt(c*c + 7*np.pi*np.pi*a*a/3.0 - 5*s*s)
 
    #Calculate form factor
    x = q2*R1
    J1 = np.sin(x)/(x*x) - np.cos(x)/x
    F = 3*J1/x
    return (F*F)*(np.exp(-(q2*q2*s*s)))


#Maximum nuclear recoil energy (in keV)
def ERmax(E_nu, A):

    #Nuclear mass in MeV
    m_A_MeV = A*0.9315e3
    
    #return 1e3*2.0*E_nu**2/m_N_MeV
    
    return 1e3*(2.0*E_nu*E_nu)/(m_A_MeV + 2*E_nu)


#----Main cross section calculation----

def xsec_CEvNS(E_R, E_nu, A, Z, gsq=0.0, m_med=1000.0):
    """
    Calculates the differential cross section for
    Coherent Elastic Neutrino-Nucleus Scattering
    including a new Z' mediator.
    
    Note: currently assuming couplings of Z'
    to u and d quarks are equal.
    
    Parameters
    ----------
    E_R : float
        Recoil energy (in keV)
    E_nu : float
        Neutrino energy (in MeV)
    A   : int
        Mass number of target nucleus
    Z   : int
        Atomic number of target nucleus

    gsq   : float, optional
        Product of couplings of the new mediator to quark
        vector current and LH-neutrino VECTOR current. 
        Set to zero by default.
    m_med   : float, optional
        Mass of new vector mediator (in MeV). Set to 1000 MeV
        by default.
    
    
    Returns
    -------
    float
        Differential scattering cross section 
        (in cm^2/keV)
    """
    
    m_A = A*0.9315 #Mass of target nucleus (in GeV)
    q = np.sqrt(2.0*E_R*m_A) #Recoil momentum (in MeV)
    #Note: m_A in GeV, E_R in keV, E_nu in MeV
    
    #Calculate SM contribution
    Qv = (A-Z) - (1.0-4.0*SIN2THETAW)*Z #Coherence factor
    
    xsec_SM = (G_FERMI*G_FERMI/(4.0*np.pi))*Qv*Qv*m_A*   \
        (1.0-(q*q)/(4.0*E_nu*E_nu))
    
    #Calculate New-Physics correction from Z' coupling
    #Assume universal coupling to quarks (u and d)
    QvNP = 3.0*A*gsq

    #Factor of 1e6 from (GeV/MeV)^2
    G_V = 1 - 1e6*(SQRT2/G_FERMI)*(QvNP/Qv)*1.0/(q*q + m_med*m_med)
    
    #Convert from (GeV^-3) to (cm^2/keV)
    #and multiply by form factor and New Physics correction
    return G_V*G_V*xsec_SM*1e-6*(1.98e-14)*(1.98e-14)*HelmFormFactor(E_R, A)
    
    
def xsec_magneticNS(E_R, E_nu, A, Z, mu_nu=0.0):
    """
    Calculates the differential cross section for
    Coherent Elastic Neutrino-Nucleus Scattering from
    neutrino magnetic moment - nuclear charge coupling.
    
    Note that in principle there could also
    be a dipole-dipole coupling between the neutrino
    and the magnetic dipole of odd-A nuclei.
    
    
    Parameters
    ----------
    E_R : float
        Recoil energy (in keV)
    E_nu : float
        Neutrino energy (in MeV)
    A   : int
        Mass number of target nucleus
    Z   : int
        Atomic number of target nucleus

    mu_nu   : float, optional
        Neutrino magnetic moment (in units of the bohr
        magneton mu_B = e/2m_e).
  
    
    Returns
    -------
    float
        Differential scattering cross section 
        (in cm^2/keV)
    """
    
    #Note: m_A in GeV, E_R in keV, E_nu in MeV
    
    #in GeV^-3
    xsec_mag = ((np.pi*ALPHA_EM**2*mu_nu**2)/m_e**2)* \
            (1.0/(1e-6*E_R) - 1.0/(E_nu*1e-3) + (1e-6*E_R)/(1e-6*4.0*E_nu**2))
    
    #Convert from (GeV^-3) to (cm^2/keV)
    #and multiply by form factor and coherent charge enhancement
    return 1e-6*(1.98e-14)**2*xsec_mag*(Z**2*HelmFormFactor(E_R, A))
    
def xsec_scalar(E_R, E_nu, A, Z, gsq=0.0, m_S=1000.0):
    """
    Calculates contribution to the differential cross section for
    Coherent Elastic Neutrino-Nucleus Scattering from
    a new scalar mediator, S. 
    
    NB: we assume equal coupling to all quarks.
    
    See arXiv:1701.07443 - Eq. (3.6).
    See arXiv:1604.01025 - Tab. IV.
    
    Parameters
    ----------
    E_R : float
        Recoil energy (in keV)
    E_nu : float
        Neutrino energy (in MeV)
    A   : int
        Mass number of target nucleus
    Z   : int
        Atomic number of target nucleus

    gsq   : float, optional
        Product of couplings of the new mediator to quark
        scalar current and neutrino scalar current. 
        For now, we assume that the coupling to all quarks
        is equal.
        Set to zero by default.
    m_S   : float, optional
        Mass of new scalar mediator (in MeV). Set to 1000 MeV
        by default.
    
    Returns
    -------
    float
        Differential scattering cross section 
        (in cm^2/keV)
    """
    
    #Note: m_A in GeV, E_R in keV, E_nu in MeV
    
    m_A = A*0.9315 #Mass of target nucleus (in GeV)
    q = np.sqrt(2.0*E_R*m_A) #Recoil momentum (in MeV)
    #Note: m_A in GeV, E_R in keV, E_nu in MeV
    
    #Calculate total scalar charge of the nucleus
    Q_S = gsq*(14.0*A + 1.1*Z)
    
    xsec_scalar = (E_R/(4.0*np.pi))*m_A**2/(E_nu**2*(q**2 + m_S**2)**2)
    
    # keV/MeV = 1e-3
        
    # (MeV)^-1 = (1e-3 GeV)^-1 = 1e3 GeV^-1
    
    # MeV^-2 keV GeV^-2 = GeV^-3
    
    #Convert from (keV GeV^2 MeV^-6) to (cm^2/keV)
    #and multiply by form factor and New Physics correction
    return Q_S*Q_S*xsec_scalar*1e6*(1.98e-14)*(1.98e-14)*HelmFormFactor(E_R, A)
    
    
def xsec_NSI(E_R, E_nu, A, Z, Eps_u_e, Eps_d_e, Eps_u_mu=0, Eps_d_mu=0, Eps_u_tau=0, Eps_d_tau=0):
    """
    Calculates the differential cross section for
    Coherent Elastic Neutrino-Nucleus Scattering
    including Non-standard Neutrino interactions (NSI)
    
    
    Parameters
    ----------
    E_R : float
        Recoil energy (in keV)
    E_nu : float
        Neutrino energy (in MeV)
    A   : int
        Mass number of target nucleus
    Z   : int
        Atomic number of target nucleus

    Eps_X_Y   : float
        NSI coupling epsilson to X type quarks (X = u, d).
        Y = e, mu, tau specifies the final neutrino flavour. 
        For flavour-conserving (i.e. nu_e -> nu_e) Y = e.
        For flavour-changing (e.g. nu_e -> nu_tau) Y = mu, tau.
    
    
    Returns
    -------
    float
        Differential scattering cross section 
        (in cm^2/keV)
    """
    
    m_A = A*0.9315 #Mass of target nucleus (in GeV)
    q = np.sqrt(2.0*E_R*m_A) #Recoil momentum (in MeV)
    #Note: m_A in GeV, E_R in keV, E_nu in MeV
    
    #Calculate NSI charge
    Q_NSI_sq = 4*((A-Z)*(-0.5 + Eps_u_e + 2*Eps_d_e) + Z*(0.5 - 2.0*SIN2THETAW + 2*Eps_u_e + Eps_d_e))**2
    Q_NSI_sq += 4*((A-Z)*(Eps_u_mu + 2*Eps_d_mu) + Z*(2*Eps_u_mu + Eps_d_mu))**2
    Q_NSI_sq += 4*((A-Z)*(Eps_u_tau + 2*Eps_d_tau) + Z*(2*Eps_u_tau + Eps_d_tau))**2
    
    xsec_NSI = (G_FERMI*G_FERMI/(4.0*np.pi))*Q_NSI_sq*m_A*   \
        (1.0-(q*q)/(4.0*E_nu*E_nu))
    
    #Convert from (GeV^-3) to (cm^2/keV)
    #and multiply by form factor and New Physics correction
    return xsec_NSI*1e-6*(1.98e-14)*(1.98e-14)*HelmFormFactor(E_R, A)
    
    
    
#----Nuclear-scattering rate calculation----

#Load neutrino fluxes and save as interpolation function
#(neutrino_flux) in units of neutrinos/cm^2/s/MeV 
#(with input E, in MeV)
def loadNeutrinoFlux(source="CHOOZ"):
    global neutrino_flux_tot
    global neutrino_flux_list
    global Enu_min
    global Enu_max
    global nu_source
    
    nu_source = source
    
    #Initialise the list of neutrino fluxes to all return zero
    #NB: the list corresponds to (e, e-bar, mu, mu-bar, tau, tau-bar)
    neutrino_flux_list = [lambda x: 1e-30 for i in range(6)]
    
    if (source == "CHOOZ"):
        data = np.loadtxt("DataFiles/neutrino_spectrum.txt")
    
        #Factor of 0.9 in normalisation comes from Reactor up-time
        normalisation = np.loadtxt("DataFiles/ScaleConstants.txt")[0]
        #Factors of 1e-3 to convert from keV to MeV
        #Electron anti-neutrinos
        neutrino_flux_list[flav['eb']] = InterpolatedUnivariateSpline(1e-3*data[:,0], 1e3*data[:,1]*normalisation, k = 1)
                #bounds_error=False, fill_value=1e-20)
        Enu_min = 1e-3*np.min(data[:,0])
        Enu_max = 1e-3*np.max(data[:,0])
        #print Enu_min, Enu_max
        
    elif (source == "SNS"):
        #0.08 decay-at-rest (DAR) neutrinos per proton
        #5e20 protons per day (divided by 24*60*60 to get per second)
        #COHERENT detector at a distance of 19.3m from neutrino production point
        #so 4pi surface area is 4pi*19.3m^2
        Nflux = 0.08*5e20/(24.0*60.0*60.0*4*np.pi*1930.0**2)
        
        E_mub, data_mub = np.loadtxt("DataFiles/SNSflux_numub.txt", unpack=True)
        rawflux_mub = InterpolatedUnivariateSpline(E_mub, data_mub, k = 1)
        
        #Normalise the raw fluxes (DAR neutrinos only go up to the end point
        # of the Michel spectrum, ~ 53.5 MeV)
        DAR_norm = quad(rawflux_mub, np.min(E_mub), 54.0)[0]
        flux_mub = InterpolatedUnivariateSpline(E_mub, data_mub*Nflux/DAR_norm, k = 1)    
        neutrino_flux_list[flav['mub']] = flux_mub
        
        E_e, data_e = np.loadtxt("DataFiles/SNSflux_nue.txt", unpack=True)
        flux_e = InterpolatedUnivariateSpline(E_e, data_e*Nflux/DAR_norm, k = 1)
        neutrino_flux_list[flav['e']] = flux_e
        
        E_mu, data_mu = np.loadtxt("DataFiles/SNSflux_numu.txt", unpack=True)
        flux_mu = InterpolatedUnivariateSpline(E_mu, data_mu*Nflux/DAR_norm, k = 1)
        neutrino_flux_list[flav['mu']] = flux_mu
    
        Enu_min = 1.0
        Enu_max = 300.0
    
    #Now tabulate the total neutrino flux
    Evals = np.logspace(np.log10(Enu_min), np.log10(Enu_max), 1000)
    flux_tab = 0.0*Evals
    for j in range(6):
        flux_tab += neutrino_flux_list[j](Evals)
    neutrino_flux_tot = InterpolatedUnivariateSpline(Evals,flux_tab, k = 1)
    
    
#Calculate recoil rate (in events/kg/keV/day)
def differentialRate_full(E_R, A, Z, gsq=0.0, m_med=1000.0, mu_nu=0.0, nu_flavor ="all"):
    """
    Calculates the differential recoil rate for
    Coherent Elastic Neutrino-Nucleus Scattering
    from Chooz reactor neutrinos.
    
    Includes contribution from both neutral current vector
    exchange and neutrino magnetic moment interaction.
    
    Checks to see whether the neutrino flux table
    has been loaded (and loads it if not...)
    
    Note: currently assuming couplings of Z'
    to u and d quarks are equal.
    
    Parameters
    ----------
    E_R : float
        Recoil energy (in keV)
    A   : int
        Mass number of target nucleus
    Z   : int
        Atomic number of target nucleus

    gsq   : float, optional
        Product of couplings of the new mediator to quark
        vector current and LH-neutrino VECTOR current. 
        Set to zero by default.
    m_med   : float, optional
        Mass of new mediator (in MeV). Set to 1000 MeV
        by default.
    mu_nu   : float, optional
        Neutrino magnetic moment (in units of the bohr
        magneton mu_B = e/2m_e).
        
    
    Returns
    -------
    float
        Differential recoil rate
        (in /kg/keV/day)
    """
    
    
    #This is not the most efficient thing, as some of the calculation
    #is repeated is repeated is repeated in the two 'differentialRate'
    #functions...but it's the easiest for now.
    return differentialRate_CEvNS(E_R, A, Z, gsq, m_med, nu_flavor = nu_flavor) + \
            differentialRate_magnetic(E_R, A, Z, mu_nu, nu_flavor = nu_flavor)
    
    
    
#Calculate recoil rate (in events/kg/keV/day)
def differentialRate_CEvNS(E_R, A, Z, gsq=0.0, m_med=1000.0, tab=False, nu_flavor="all"):
    """
    Calculates the differential recoil rate for
    Coherent Elastic Neutrino-Nucleus Scattering
    from Chooz reactor neutrinos.
    
    Checks to see whether the neutrino flux table
    has been loaded (and loads it if not...)
    
    Note: currently assuming couplings of Z'
    to u and d quarks are equal.
    
    Parameters
    ----------
    E_R : float
        Recoil energy (in keV)
    A   : int
        Mass number of target nucleus
    Z   : int
        Atomic number of target nucleus

    gsq   : float, optional
        Product of couplings of the new mediator to quark
        vector current and LH-neutrino VECTOR current. 
        Set to zero by default.
    m_med   : float, optional
        Mass of new mediator (in MeV). Set to 1000 MeV
        by default.
    
    
    Returns
    -------
    float
        Differential recoil rate
        (in /kg/keV/day)
    """
    
    #Use tabulated values if requested. 
    #Note: no bounds checking, or checking if tabulation is alread done!
    if (tab):
        A1 = diffrate_A(E_R)
        B1 = diffrate_B(E_R)
        return (A1 + B1*gsq)**2.0
    
    
    #First, check that the neutrino flux has been loaded
    if (neutrino_flux_tot == None):
        print " CEvNS.py: Loading neutrino flux for the first time..."
        loadNeutrinoFlux()
    
    if (nu_flavor == "all"):
        integrand = lambda E_nu: xsec_CEvNS(E_R, E_nu, A, Z, gsq, m_med)\
                        *neutrino_flux_tot(E_nu)
    else:
        integrand = lambda E_nu: xsec_CEvNS(E_R, E_nu, A, Z, gsq, m_med)\
                        *neutrino_flux_list[flav[nu_flavor]](E_nu)
    
    #Minimum neutrino energy required (in MeV)
    E_min = np.sqrt(A*0.9315*E_R/2)
    
    E_min = np.maximum(E_min, Enu_min)
    
    #For reactor neutrinos, set E_max:
    E_max = Enu_max
    
    if (E_min > E_max):
        return 0
    
    m_N = A*1.66054e-27 #Nucleus mass in kg
    rate = quad(integrand, E_min, E_max, epsrel=1e-4)[0]/m_N
    
    if (nu_source == "SNS" and (nu_flavor == "all" or nu_flavor == "mu")):
        if (E_min < 29.65): #Add in delta function neutrino flux from muon decay
            #Check the loadNeutrinoFlux() for a descripion of this normalisation
            Nflux = 0.08*5e20/(24.0*60.0*60.0*4*np.pi*1930.0**2)
            rate += Nflux*xsec_CEvNS(E_R, 29.65, A, Z, gsq, m_med)/m_N
    
    return 86400.0*rate #Convert from (per second) to (per day)
    
#Calculate recoil rate (in events/kg/keV/day)
def differentialRate_NSI(E_R, A, Z, Eps_u_e, Eps_d_e, Eps_u_mu=0, Eps_d_mu=0, Eps_u_tau=0, Eps_d_tau=0, nu_flavor="all"):
    """
    Calculates the differential recoil rate for
    Coherent Elastic Neutrino-Nucleus Scattering
    from Chooz reactor neutrinos (including 
    Non-standard Neutrino interactions (NSI)).
    
    Checks to see whether the neutrino flux table
    has been loaded (and loads it if not...)
    
    Parameters
    ----------
    E_R : float
        Recoil energy (in keV)
    A   : int
        Mass number of target nucleus
    Z   : int
        Atomic number of target nucleus

    Eps_X_Y   : float
        NSI coupling epsilson to X type quarks (X = u, d).
        Y = e, mu, tau specifies the final neutrino flavour. 
        For flavour-conserving (i.e. nu_e -> nu_e) Y = e.
        For flavour-changing (e.g. nu_e -> nu_tau) Y = mu, tau.
    
    
    Returns
    -------
    float
        Differential recoil rate
        (in /kg/keV/day)
    """    
    
    #First, check that the neutrino flux has been loaded
    if (neutrino_flux_tot == None):
        print " CEvNS.py: Loading neutrino flux for the first time..."
        loadNeutrinoFlux()
    
    if (nu_flavor == "all"):
        integrand = lambda E_nu: xsec_NSI(E_R, E_nu, A, Z, Eps_u_e, Eps_d_e,\
                                 Eps_u_mu, Eps_d_mu, Eps_u_tau, Eps_d_tau)\
                                 *neutrino_flux_tot(E_nu)
    else:
        integrand = lambda E_nu: xsec_NSI(E_R, E_nu, A, Z, Eps_u_e, Eps_d_e,\
                                 Eps_u_mu, Eps_d_mu, Eps_u_tau, Eps_d_tau)\
                                 *neutrino_flux_list[flav[nu_flavor]](E_nu)
    
    #Minimum neutrino energy required (in MeV)
    E_min = np.sqrt(A*0.9315*E_R/2)
    
    E_min = np.maximum(E_min, Enu_min)
    
    #For reactor neutrinos, set E_max:
    E_max = Enu_max
    
    if (E_min > E_max):
        return 0
    
    m_N = A*1.66054e-27 #Nucleus mass in kg
    rate = quad(integrand, E_min, E_max, epsrel=1e-4)[0]/m_N
    
    if (nu_source == "SNS" and (nu_flavor == "all" or nu_flavor == "mu")):
        if (E_min < 29.65): #Add in delta function neutrino flux from muon decay
            #Check the loadNeutrinoFlux() for a descripion of this normalisation
            Nflux = 0.08*5e20/(24.0*60.0*60.0*4*np.pi*1930.0**2)
            rate += Nflux*xsec_NSI(E_R, 29.65, A, Z, Eps_u_e, Eps_d_e,\
                                 Eps_u_mu, Eps_d_mu, Eps_u_tau, Eps_d_tau)/m_N
            
    return 86400.0*rate #Convert from (per second) to (per day)
    
#Calculate recoil rate (in events/kg/keV/day)
def differentialRate_scalar(E_R, A, Z, gsq=0.0, m_S=1000.0, tab=False, nu_flavor="all"):
    """
    Calculates the differential recoil rate for
    Coherent Elastic Neutrino-Nucleus Scattering 
    from a new scalar mediator, S. 
    
    NB: we assume equal coupling to all quarks.
    
    See arXiv:1701.07443 - Eq. (3.6).
    See arXiv:1604.01025 - Tab. IV.
    
    Checks to see whether the neutrino flux table
    has been loaded (and loads it if not...)
    
    Parameters
    ----------
    E_R : float
        Recoil energy (in keV)
    A   : int
        Mass number of target nucleus
    Z   : int
        Atomic number of target nucleus

    gsq   : float, optional
        Product of couplings of the new mediator to quark
        scalar current and neutrino scalar current. 
        For now, we assume that the coupling to all quarks
        is equal.
        Set to zero by default.
    m_S   : float, optional
        Mass of new scalar mediator (in MeV). Set to 1000 MeV
        by default.
    
    Returns
    -------
    float
        Differential scattering cross section 
        (in cm^2/keV)
    """
    
    #First, check that the neutrino flux has been loaded
    if (neutrino_flux_tot == None):
        print " CEvNS.py: Loading neutrino flux for the first time..."
        loadNeutrinoFlux()
    
    if (nu_flavor == "all"):
        integrand = lambda E_nu: xsec_scalar(E_R, E_nu, A, Z, gsq, m_S)\
                        *neutrino_flux_tot(E_nu)
    else:
        integrand = lambda E_nu: xsec_scalar(E_R, E_nu, A, Z, gsq, m_S)\
                        *neutrino_flux_list[flav[nu_flavor]](E_nu)
    
    #Minimum neutrino energy required (in MeV)
    E_min = np.sqrt(A*0.9315*E_R/2)
    
    E_min = np.maximum(E_min, Enu_min)
    
    #For reactor neutrinos, set E_max:
    E_max = Enu_max
    
    if (E_min > E_max):
        return 0
    
    m_N = A*1.66054e-27 #Nucleus mass in kg
    rate = quad(integrand, E_min, E_max, epsrel=1e-4)[0]/m_N
    
    if (nu_source == "SNS" and (nu_flavor == "all" or nu_flavor == "mu")):
        if (E_min < 29.65): #Add in delta function neutrino flux from muon decay
            #Check the loadNeutrinoFlux() for a descripion of this normalisation
            Nflux = 0.08*5e20/(24.0*60.0*60.0*4*np.pi*1930.0**2)
            rate += Nflux*xsec_scalar(E_R, 29.65, A, Z, gsq, m_S)/m_N
            
    return 86400.0*rate #Convert from (per second) to (per day)
    
#Calculate recoil rate (in events/kg/keV/day)
def differentialRate_magnetic(E_R, A, Z, mu_nu=0.0, nu_flavor="all"):
    """
    Calculates the differential recoil rate for
    Coherent Elastic Neutrino-Nucleus Scattering
    from Chooz reactor neutrinos.
    
    Checks to see whether the neutrino flux table
    has been loaded (and loads it if not...)
    
    Note: currently assuming couplings of Z'
    to u and d quarks are equal.
    
    Parameters
    ----------
    E_R : float
        Recoil energy (in keV)
    A   : int
        Mass number of target nucleus
    Z   : int
        Atomic number of target nucleus

    mu_nu   : float, optional
        Neutrino magnetic moment (in units of the bohr
        magneton mu_B = e/2m_e).
    
    
    Returns
    -------
    float
        Differential recoil rate
        (in /kg/keV/day)
    """
    
    
    #First, check that the neutrino flux has been loaded
    if (neutrino_flux_tot == None):
        loadNeutrinoFlux()
    
    if (nu_flavor == "all"):
        integrand = lambda E_nu: xsec_magneticNS(E_R, E_nu, A, Z, mu_nu)\
                        *neutrino_flux_tot(E_nu)
    else:
        integrand = lambda E_nu: xsec_magneticNS(E_R, E_nu, A, Z, mu_nu)\
                        *neutrino_flux_list[flav[nu_flavor]](E_nu)  
    
    #Minimum neutrino energy required (in MeV)
    E_min = np.sqrt(A*0.9315*E_R/2)
    
    E_min = np.maximum(E_min, Enu_min)
    
    #For reactor neutrinos, set E_max:
    E_max = Enu_max
    
    if (E_min > E_max):
        return 0
    
    m_N = A*1.66054e-27 #Nucleus mass in kg
    rate = quad(integrand, E_min, E_max)[0]/m_N
    
    if (nu_source == "SNS" and (nu_flavor == "all" or nu_flavor == "mu")):
        if (E_min < 29.65): #Add in delta function neutrino flux from muon decay
            #Check the loadNeutrinoFlux() for a descripion of this normalisation
            Nflux = 0.08*5e20/(24.0*60.0*60.0*4*np.pi*1930.0**2)
            rate += Nflux*xsec_magneticNS(E_R, 29.65, A, Z, mu_nu)/m_N
    
    return 86400*rate #Convert from (per second) to (per day)
    
#-----------------
def differentialRate_tabulate(E_min, E_max, A, Z, m_med=1000.0, nu_flavor='all'):
    global diffrate_A, diffrate_B
    
    print " Tabulating dR/dE for m_med =", '%.2e' % m_med, "MeV..."
    Evals = np.logspace(np.log10(E_min), np.log10(E_max),1000)
    
    Rvals_A = 0.0*Evals
    Rvals_B = 0.0*Evals
    alpha = 1.0
    for i, Ei in enumerate(Evals):
        Rvals_A[i] = np.sqrt(differentialRate_CEvNS(Ei, A, Z, gsq=0.0, m_med=m_med,nu_flavor=nu_flavor))
        Rvals_B[i] = (1.0/(4.0*alpha*Rvals_A[i]))*(\
                        differentialRate_CEvNS(Ei, A, Z, gsq=alpha, m_med=m_med,nu_flavor=nu_flavor) - \
                        differentialRate_CEvNS(Ei, A, Z, gsq=-alpha, m_med=m_med,nu_flavor=nu_flavor))
        #Check signs!
        #Rvals_B[i] = np.sqrt(differentialRate_CEvNS(Ei, A, Z, gsq=Rvals_A[i], m_med=m_med))*1.0/Rvals_A[i] - 1.0      
   
    diffrate_A = InterpolatedUnivariateSpline(Evals, Rvals_A, k = 1)
    diffrate_B = InterpolatedUnivariateSpline(Evals, Rvals_B, k = 1)
    
    
    
    
    