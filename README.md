# CEvNS

[![Binder](http://mybinder.org/badge.svg)](https://beta.mybinder.org/v2/gh/bradkav/CEvNS/master)

Code for calculating Coherent Elastic Neutrino-Nucleus Scattering (CEvNS) cross sections. Cross sections taken largely from [arXiv:1604.01025](https://arxiv.org/abs/1604.01025) and [arXiv:1701.07443](https://arxiv.org/abs/1701.07443).

See the **example code** in the iPython notebook (`CEvNS-examples.ipynb`) [here](https://nbviewer.jupyter.org/github/bradkav/CEvNS/blob/master/CEvNS-examples.ipynb) or click the 'Launch Binder' button to view an interactive notebook. For working with the **recent COHERENT data** ([arXiv:1708.01294](https://arxiv.org/abs/1708.01294)), check out the notebook `COHERENT.ipynb` [here](https://nbviewer.jupyter.org/github/bradkav/CEvNS/blob/master/COHERENT.ipynb).

Currently includes (amonth other things) the Standard Model contribution to the CEvNS cross section, along with the contribution from the following Simplified Model Lagrangian:

<img src="/L1.png" width="400">

Note that the CEvNS cross section depends only on the product of g_v and g_q, so in the code the only coupling is called gsq and is the product of these two. More interactions (magnetic moments, scalar mediators, etc.) are being added.

If you are using the NSI rate, you do not need to add the Standard Model rate too - the NSI rate includes (and modifies) the normalisation of Standard Model CEvNS.

### Version History

**Version 1.5 (21/09/2017):** Now stores separate fluxes for different neutrino species (specify in the differential rate functions using nu_flavor=('e', 'eb', 'mu', 'mub', 'tau', 'taub', or 'all')).  
**Version 1.4 (20/09/2017):** Now includes COHERENT@SNS fluxes (simply use `loadNeutrinoFlux("SNS")` to initialise). Results from the COHERENT experiment are reproduced in the notebook `COHERENT.ipynb`.  
**Version 1.3 (13/09/2017):** Updated to calculate rate including NSI interactions.  
**Version 1.2 (15/02/2017):** Updated to include contribution to nuclear scattering from 'neutrino magnetic moment'-'nuclear charge' coupling.  
**Version 1.1 (31/01/2017):** Updated to include *reactor neutrino flux* and calculation of the differential neutrino-nucleus scattering rate.  
**Version 1.0 (25/01/2017):** Initial code include SM contribution and vector current Simplified Model.  


### Using the code

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

Contribution from NSI (essentially a modification of SM rate):  
`xsec_NSI(E_R, E_nu, A, Z, Eps_u_e, Eps_d_e, Eps_u_mu=0, Eps_d_mu=0, Eps_u_tau=0, Eps_d_tau=0)`

where  
- `E_R` is the nuclear recoil energy (in keV)
- `E_nu` is the neutrino energy (in MeV)
- (`A`, `Z`) are the mass and atomic number of the target nucleus
- `Eps_X_Y` is the NSI coupling epsilson to X type quarks (X = u, d). Y = e, mu, tau specifies the final neutrino flavour. For flavour-conserving (i.e. nu_e -> nu_e) Y = e. For flavour-changing (e.g. nu_e -> nu_tau) Y = mu, tau.

#### Differential rates

Calculate differential rates in units of counts/keV/kg/day (parameter definitions as above).

Contribution from vector exchange (Z and possible new Z' mediator):  
`differentialRate_CEvNS(E_R, A, Z, gsq=0.0, m_med=1000.0, nu_flavor='all')`

Contribution from NSI (modification of SM rate):
`differentialRate_NSI(E_R, A, Z, Eps_u_e, Eps_d_e, Eps_u_mu=0, Eps_d_mu=0, Eps_u_tau=0, Eps_d_tau=0, nu_flavor='all')`

Contribution from neutrino magnetic moment:  
`differentialRate_magnetic(E_R, A, Z, mu_nu=0.0, nu_flavor='all')`

Full contribution from vector exchange and magnetic moment:  
`differentialRate_full(E_R, A, Z, gsq=0.0, m_med=1000.0,mu_nu=0.0, nu_flavor='all')`

#### Neutrino fluxes

Before you can calculate the differential rates, you need to load an appropriate set of neutrino fluxes using `loadNeutrinoFlux(source)`. There are two options for `source`: `"CHOOZ"` (for the CHOOZ reactor) and `"SNS"` (for the Spallation Neutron Source where COHERENT is located). `"CHOOZ"` is the default choice, and if you don't call `loadNeutrinoFlux` before you try to calculate a differential rate, the CHOOZ neutrino flux is automatically loaded.

The fluxes for the 6 different neutrino species (e, eb, mu, mub, tau, taub) are stored separately, as well as the total neutrino flux. When you call a differential rate function, you can specify which species you're interested in with nu_flavor=('e', 'eb', 'mu', 'mub', 'tau', 'taub', or 'all'). The default is 'all', which sums all the neutrino fluxes.
