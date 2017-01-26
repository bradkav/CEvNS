# CEvNS

Code for calculating Coherent Elastic Neutrino-Nucleus Scattering (CEvNS) cross sections. Cross sections taken largely from [arXiv:1604.01025](https://arxiv.org/abs/1604.01025).

See the **example code** in the iPython notebook (`CEvNS_example.ipynb`) [here](https://nbviewer.jupyter.org/github/bradkav/CEvNS/blob/master/CEvNS_example.ipynb).

Currently includes the Standard Model contribution to the CEvNS cross section, along with the contribution from the following Simplified Model Lagrangian:

<img src="/L1.png" width="400">

Note that the CEvNS cross section depends only on the product of g_v and g_q, so in the code these are set to be equal (and called g_med).

The code will be updated soon to include more general (e.g. axial-vector current) interactions.

**Version 1.0 (25/01/2017):** Initial code include SM contribution and vector current Simplified Model.

###Using the code

Import the module `CEvNS.py`. This gives you access to the function:

`xsec_CEvNS(E_R, E_nu, A, Z, g_med=0.0, m_med=1000.0)`

where 
- `E_R` is the nuclear recoil energy (in keV)
- `E_nu` is the neutrino energy (in MeV)
- (`A`, `Z`) are the mass and atomic number of the target nucleus
- `g_med` is the mediator coupling to neutrinos and quarks (dimensionless)
- `m_med` is the mediator mass (in MeV).
