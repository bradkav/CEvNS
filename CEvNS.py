"""
CEvNS.py - Version 1 - 25/01/2017

Summary: 
Code for calculating differential cross section
for Coherent Elastic Neutrino Nucleus Scattering (CEvNS).

Cross sections taken from arXiv:1604.01025

Author: Bradley J Kavanagh
Please report any problems to: bradkav@gmail.com
"""

import numpy as np

#----Notes-----
"""
We're assuming that the Z-mass is much larger
than the typical momentum transfer (which is fine
for solar and reactor neutrinos).
"""


#----Constants----
G_Fermi = 1.1664e-5     #Fermi Constant in GeV^-2
sin2thetaW = 0.2387     #Sine-squared of the weak angle


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

def xsec_CEvNS(E_R, E_nu, A, Z, g_med=0.0, m_med=1000.0):
    """
    Calculates the differential cross section for
    Coherent Elastic Neutrino-Nucleus Scattering.
    
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

    g_med   : float, optional
        Coupling of the new mediator to quarks and neutrinos
        (assumed equal). Set to zero by default.
    m_med   : float, optional
        Mass of new mediator (in MeV). Set to 1000 MeV
        by default.
    
    
    Returns
    -------
    float
        Differential scattering cross section 
        (in cm^2/keV)
    """
    
    #Mass of target nucleus (in GeV)
    m_A = A*0.9315
    #Note: m_A in GeV, E_R in keV, E_nu in MeV
    #so m_A*E_R is in MeV too.
    
    #Calculate SM contribution
    Qv = (A-Z) - (1.0-4.0*sin2thetaW)*Z #Coherence factor
    
    xsec_SM = (G_Fermi**2/(4.0*np.pi))*Qv**2*m_A*   \
        (1.0-(m_A*E_R)/(2.0*E_nu**2))
    
    #Calculate Z' contribution
    Qvp = 3.0*A*g_med**2
    xsec_NP = -1e6*G_Fermi*m_A*Qv*Qvp*(2.0*E_nu**2 - E_R*m_A)  \
        /(2*np.sqrt(2)*np.pi*E_nu**2*(2.0*E_R*m_A+m_med**2))
    #Factor of 1e6 to get it into units of 1/GeV^3
    
    xsec_NP += 1e12*Qvp**2*m_A*(2.0*E_nu**2 - E_R*m_A)   \
        /(4*np.pi*E_nu**2*(2.0*E_R*m_A+m_med**2)**2)
    #Factor of 1e12 to get it into units of 1/GeV^3
    
    #Need to double check whether form factor appears for 
    #the interference terms
    
    #Convert from (GeV^-3) to (cm^2/keV)
    #and multiply by form factor
    return (xsec_SM + xsec_NP)*1e-6*(1.98e-14)**2*HelmFormFactor(E_R, A)
    
    