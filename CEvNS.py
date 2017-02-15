"""
CEvNS.py - Version 1.2 - 15/02/2017

Summary: 
Code for calculating differential cross section
for Coherent Elastic Neutrino Nucleus Scattering (CEvNS).

Cross sections mainly taken from arXiv:1604.01025.
See also arXiv:1701.07443.

Magnetic moment - arXiv:1508.07981

Requires numpy and scipy.

Author: Bradley J Kavanagh
Please report any problems to: bradkav@gmail.com
"""

import numpy as np
from scipy.interpolate import interp1d
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


#----Module-scope variables----
neutrino_flux = None


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
    J1 = np.sin(x)/x**2 - np.cos(x)/x
    F = 3*J1/x
    return (F**2)*(np.exp(-(q2*s)**2))


#Maximum nuclear recoil energy (in keV)
def ERmax(E_nu, A):

    #Nuclear mass in MeV
    m_A_MeV = A*0.9315e3
    
    #return 1e3*2.0*E_nu**2/m_N_MeV
    
    return 1e3*(2*E_nu**2)/(m_A_MeV + 2*E_nu)


#----Main cross section calculation----

def xsec_CEvNS(E_R, E_nu, A, Z, gsq=0.0, m_med=1000.0):
    """
    Calculates the differential cross section for
    Coherent Elastic Neutrino-Nucleus Scattering.
    
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
        Mass of new mediator (in MeV). Set to 1000 MeV
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
    
    xsec_SM = (G_FERMI**2/(4.0*np.pi))*Qv**2*m_A*   \
        (1.0-(q**2)/(4.0*E_nu**2))
    
    #Calculate New-Physics correction from Z' coupling
    #Assume universal coupling to quarks (u and d)
    QvNP = 3.0*A*gsq

    #Factor of 1e6 from (GeV/MeV)^2
    G_V = 1 - 1e6*(np.sqrt(2)/G_FERMI)*(QvNP/Qv)*1.0/(q**2 + m_med**2)
    
    #Convert from (GeV^-3) to (cm^2/keV)
    #and multiply by form factor and New Physics correction
    return G_V**2.0*xsec_SM*1e-6*(1.98e-14)**2*HelmFormFactor(E_R, A)
    
    
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
    
    
#----Nuclear-scattering rate calculation----

#Load neutrino fluxes and save as interpolation function
#(neutrino_flux) in units of neutrinos/cm^2/s/MeV 
#(with input E, in MeV)
def loadNeutrinoFlux():
    global neutrino_flux
    
    data = np.loadtxt("DataFiles/neutrino_spectrum.txt")
    normalisation = np.loadtxt("DataFiles/ScaleConstants.txt")[0]
    
    #Factors of 1e-3 to convert from keV to MeV
    neutrino_flux = interp1d(1e-3*data[:,0], 1e3*data[:,1]*normalisation,\
            bounds_error=False, fill_value=0)
    
    
#Calculate recoil rate (in events/kg/keV/day)
def differentialRate_full(E_R, A, Z, gsq=0.0, m_med=1000.0, mu_nu=3.2e-11):
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
    return differentialRate_CEvNS(E_R, A, Z, gsq, m_med) + \
            differentialRate_magnetic(E_R, A, Z, mu_nu)
    
    
    
#Calculate recoil rate (in events/kg/keV/day)
def differentialRate_CEvNS(E_R, A, Z, gsq=0.0, m_med=1000.0):
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
    
    
    #First, check that the neutrino flux has been loaded
    if (neutrino_flux == None):
        loadNeutrinoFlux()
    
    integrand = lambda E_nu: xsec_CEvNS(E_R, E_nu, A, Z, gsq, m_med)\
                        *neutrino_flux(E_nu)
    
    #Minimum neutrino energy required (in MeV)
    E_min = np.sqrt(A*0.9315*E_R/2)
    
    #For reactor neutrinos, set E_max:
    E_max = 11.0
    
    m_N = A*1.66054e-27 #Nucleus mass in kg
    rate = quad(integrand, E_min, E_max)[0]/m_N
    return 86400*rate #Convert from (per second) to (per day)
    
#Calculate recoil rate (in events/kg/keV/day)
def differentialRate_magnetic(E_R, A, Z, mu_nu=3.2e-11):
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
    if (neutrino_flux == None):
        loadNeutrinoFlux()
    
    integrand = lambda E_nu: xsec_magneticNS(E_R, E_nu, A, Z, mu_nu)\
                        *neutrino_flux(E_nu)
    
    #Minimum neutrino energy required (in MeV)
    E_min = np.sqrt(A*0.9315*E_R/2)
    
    #For reactor neutrinos, set E_max:
    E_max = 11.0
    
    m_N = A*1.66054e-27 #Nucleus mass in kg
    rate = quad(integrand, E_min, E_max)[0]/m_N
    return 86400*rate #Convert from (per second) to (per day)
    