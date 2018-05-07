# CEvNS

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/bradkav/CEvNS/master?filepath=CEvNS-examples.ipynb) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1229055.svg)](https://doi.org/10.5281/zenodo.1229055) [![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php)

*Code for calculating Coherent Elastic Neutrino-Nucleus Scattering (CEvNS) cross sections and recoil spectra.*

Currently includes (among other things) the Standard Model contribution to the CEvNS cross section, along with the contribution from Simplified Models with new vector or scalar mediators. The code also covers neutrino magnetic moments and non-standard contact neutrino interactions (NSI). Cross sections taken largely from [arXiv:1604.01025](https://arxiv.org/abs/1604.01025) and [arXiv:1701.07443](https://arxiv.org/abs/1701.07443).

#### Requirements

The code is compatible with Python2.7 and requires [numpy](http://www.numpy.org) and [scipy](https://www.scipy.org). 

Some of the associated notebooks require [matplotlib](https://matplotlib.org) and [tqdm](https://github.com/tqdm/tqdm).

#### Getting started

To get started, see the **example code** in the jupyter notebook [`CEvNS-examples.ipynb`](https://nbviewer.jupyter.org/github/bradkav/CEvNS/blob/master/CEvNS-examples.ipynb) or click the 'Launch Binder' button to view an interactive notebook. A more detailed description of the code is also shown below.

For working with the **recent COHERENT data** ([arXiv:1708.01294](https://arxiv.org/abs/1708.01294)), check out the notebook [`COHERENT.ipynb`](https://nbviewer.jupyter.org/github/bradkav/CEvNS/blob/master/COHERENT.ipynb).


#### Related projects

This code was used to generate plots and perform some calculations for the following projects:

> J Billard, J Johnston & B J Kavanagh, "*Prospects for exploring New Physics in
 Coherent Elastic Neutrino-Nucleus Scattering*" (2018), [arXiv:1805.01798](https://arxiv.org/abs/1805.01798)

#### Contact

 If you have any questions, comments, feature requests, bug-reports etc., please contact Bradley Kavanagh at bradkav@gmail.com. 

# Using the code

Import the module `CEvNS.py`:

```python
import CEvNS
```

This gives you access to the following functions for calculating the cross section and differential recoil rate. 

#### Cross sections

Calculate differential cross sections in units of cm^2/keV.

*Contribution from vector exchange (Z and possible new Z' mediator):*  
`xsec_CEvNS(E_R, E_nu, A, Z, gsq=0.0, m_med=1000.0)`

where  
- `E_R` is the nuclear recoil energy (in keV)
- `E_nu` is the neutrino energy (in MeV)
- (`A`, `Z`) are the mass and atomic number of the target nucleus
- `gsq` is the product of the mediator coupling to LH neutrinos and quark vector currents (g_nu*g_q, dimensionless)
- `m_med` is the mediator mass (in MeV).

*Contribution from scalar exchange (a new mediator S):*  
`xsec_scalar(E_R, E_nu, A, Z, gsq=0.0, m_S=1000.0)`

where  
- `E_R` is the nuclear recoil energy (in keV)
- `E_nu` is the neutrino energy (in MeV)
- (`A`, `Z`) are the mass and atomic number of the target nucleus
- `gsq` is the product of scalar coupling to neutrinos and quark scalar currents (g_nu*g_q, dimensionless)
- `m_S` is the mediator mass (in MeV).

*Contribution from neutrino magnetic moment:*  
`xsec_magneticNS(E_R, E_nu, A, Z, mu_nu=0.0)`

where  
- `E_R` is the nuclear recoil energy (in keV)
- `E_nu` is the neutrino energy (in MeV)
- (`A`, `Z`) are the mass and atomic number of the target nucleus
- `mu_nu` is the neutrino magnetic moment in units of the Bohr Magneton mu_B.

*Contribution from NSI (essentially a modification of SM rate):*  
`xsec_NSI(E_R, E_nu, A, Z, Eps_u_e, Eps_d_e, Eps_u_mu=0, Eps_d_mu=0, Eps_u_tau=0, Eps_d_tau=0)`

where  
- `E_R` is the nuclear recoil energy (in keV)
- `E_nu` is the neutrino energy (in MeV)
- (`A`, `Z`) are the mass and atomic number of the target nucleus
- `Eps_X_Y` is the NSI coupling epsilson to X type quarks (X = u, d). Y = e, mu, tau specifies the final neutrino flavour. For flavour-conserving (i.e. nu_e -> nu_e) Y = e. For flavour-changing (e.g. nu_e -> nu_tau) Y = mu, tau.

#### Neutrino fluxes

Before you can calculate the differential rates, you need to load an appropriate set of neutrino fluxes using `loadNeutrinoFlux(source)`. There are two options for `source`: `"CHOOZ"` (for the CHOOZ reactor, Near Site) and `"SNS"` (for the Spallation Neutron Source where COHERENT is located). `"CHOOZ"` is the default choice, and if you don't call `loadNeutrinoFlux` before you try to calculate a differential rate, the CHOOZ neutrino flux is automatically loaded.

The fluxes for the 6 different neutrino species (e, eb, mu, mub, tau, taub) are stored separately, as well as the total neutrino flux. When you call a differential rate function, you can specify which species you're interested in with nu_flavor=('e', 'eb', 'mu', 'mub', 'tau', 'taub', or 'all'). The default is 'all', which sums the contributions from all the neutrino flavours.

#### Differential rates

Calculate differential recoil rates (per unit detector mass) in units of counts/keV/kg/day (parameter definitions as above).

*Contribution from vector exchange (Z and possible new Z' mediator):*  
`differentialRate_CEvNS(E_R, A, Z, gsq=0.0, m_med=1000.0, nu_flavor='all')`

*Contribution from scalar exchange (new S mediator):*  
`differentialRate_scalar(E_R, A, Z, gsq=0.0, m_S=1000.0, nu_flavor='all')`

*Contribution from NSI (modification of SM rate):*  
`differentialRate_NSI(E_R, A, Z, Eps_u_e, Eps_d_e, Eps_u_mu=0, Eps_d_mu=0, Eps_u_tau=0, Eps_d_tau=0, nu_flavor='all')`

*Contribution from neutrino magnetic moment:*  
`differentialRate_magnetic(E_R, A, Z, mu_nu=0.0, nu_flavor='all')`

*Full contribution from vector exchange and magnetic moment:*  
`differentialRate_full(E_R, A, Z, gsq=0.0, m_med=1000.0,mu_nu=0.0, nu_flavor='all')`

Note that if you are using the NSI or vector mediator rates, you do not need to add the Standard Model rate too - the NSI and vector rates include (and modify) the normalisation of Standard Model CEvNS. For the scalar mediator and neutrino magnetic moments, you need to add also the contribution of the Standard Model rate.

# Version History

**Release Version 1.0 (05/05/2018):** Release version, including plotting and examples routines, associated with [arXiv:1805.?????](https://arxiv.org/abs/1805.?????).

**Pre-release Version 1.6 (06/11/2017):** Updated and checked code for Scalar cross section and differential rate.  
**Pre-release Version 1.5 (21/09/2017):** Now stores separate fluxes for different neutrino species (specify in the differential rate functions using nu_flavor=('e', 'eb', 'mu', 'mub', 'tau', 'taub', or 'all')).  
**Pre-release Version 1.4 (20/09/2017):** Now includes COHERENT@SNS fluxes (simply use `loadNeutrinoFlux("SNS")` to initialise). Results from the COHERENT experiment are reproduced in the notebook `COHERENT.ipynb`.  
**Pre-release Version 1.3 (13/09/2017):** Updated to calculate rate including NSI interactions.  
**Pre-release Version 1.2 (15/02/2017):** Updated to include contribution to nuclear scattering from 'neutrino magnetic moment'-'nuclear charge' coupling.  
**Pre-release Version 1.1 (31/01/2017):** Updated to include *reactor neutrino flux* and calculation of the differential neutrino-nucleus scattering rate.  
**Pre-release Version 1.0 (25/01/2017):** Initial code include SM contribution and vector current Simplified Model.  


# License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
