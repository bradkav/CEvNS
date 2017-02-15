# CEvNS

[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org:/repo/bradkav/cevns)

Code for calculating Coherent Elastic Neutrino-Nucleus Scattering (CEvNS) cross sections. Cross sections taken largely from [arXiv:1604.01025](https://arxiv.org/abs/1604.01025) and [arXiv:1701.07443](https://arxiv.org/abs/1701.07443).

See the **example code** in the iPython notebook (`index.ipynb`) [here](https://nbviewer.jupyter.org/github/bradkav/CEvNS/blob/master/index.ipynb).

Currently includes the Standard Model contribution to the CEvNS cross section, along with the contribution from the following Simplified Model Lagrangian:

<img src="/L1.png" width="400">

Note that the CEvNS cross section depends only on the product of g_v and g_q, so in the code the only coupling is called gsq and is the product of these two.

**Version 1.2 (15/02/2017):** Updated to include contribution to nuclear scattering from 'neutrino magnetic moment'-'nuclear charge' coupling.
**Version 1.1 (31/01/2017):** Updated to include *reactor neutrino flux* and calculation of the differential neutrino-nucleus scattering rate. 
**Version 1.0 (25/01/2017):** Initial code include SM contribution and vector current Simplified Model.  


###Using the code

Import the module `CEvNS.py`. This gives you access to the following functions for calculating the cross section and differential recoil rate. 

#### Cross sections

Calculate differential cross sections in units of cm^2/keV.

Contribution from vector exchange (Z and possible new Z' mediator):  
`xsec_CEvNS(E_R, E_nu, A, Z, g_med=0.0, m_med=1000.0)`

where  
- `E_R` is the nuclear recoil energy (in keV)
- `E_nu` is the neutrino energy (in MeV)
- (`A`, `Z`) are the mass and atomic number of the target nucleus
- `gsq` is the product mediator coupling to LH neutrinos and quark vector currents (dimensionless)
- `m_med` is the mediator mass (in MeV).

Contribution from neutrino magnetic moment:  
`xsec_magneticNS(E_R, E_nu, A, Z, mu_nu=0.0)`

where  
- `E_R` is the nuclear recoil energy (in keV)
- `E_nu` is the neutrino energy (in MeV)
- (`A`, `Z`) are the mass and atomic number of the target nucleus
- `mu_nu` is the neutrino magnetic moment in units of the Bohr Magneton mu_B.

#### Differential rates

Calculate differential rates in units of counts/keV/kg/day (parameter definitions as above).

Contribution from vector exchange (Z and possible new Z' mediator):  
`differentialRate_CEvNS(E_R, A, Z, gsq=0.0, m_med=1000.0)`

Contribution from neutrino magnetic moment:  
`differentialRate_magnetic(E_R, A, Z, mu_nu=0.0)`

Full contribution from vector exchange and magnetic moment:  
`differentialRate_full(E_R, A, Z, gsq=0.0, m_med=1000.0,mu_nu=0.0)`
