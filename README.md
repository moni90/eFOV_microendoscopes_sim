# eFOV_microendoscopes_sim
Simulations code for Antonini et al, 2020, bioRxiv

This code generates artificial time series from synthetic data. It randomly places somas in a given volumes and simulates their spiking activity. Then, simulates the imaging of the artificial volume with two different types of endoscopes, one without correction for aberrations and one with aberration correction.

The accompanying paper is:

**Extended field-of-view ultrathin microendoscopes for high-resolution two-photon imaging with minimal invasiveness in awake mice**

Andrea Antonini, Andrea Sattin, Monica Moroni, Serena Bovetti, Claudio Moretti, Francesca Succol, Angelo Forli, Dania Vecchia, Vijayakumar P. Rajamanickam, Andrea Bertoncini, Stefano Panzeri, Carlo Liberale, Tommaso Fellin

2020, eLife, doi: https://doi.org/10.7554/eLife.58882

Before running the code for simulations, some auxiliary files must be downloaded from DOI: 10.17632/wm6c5wzs4c.2 (dataset_3, folders imaging_ellipsoids, intensity_mask, psf_estimate, ruler_FOV) and added to the main folder.

Parameters for the simulations are set in main_simulate_TSeries_only.m 

Details:
- required software: Matlab R2019b
