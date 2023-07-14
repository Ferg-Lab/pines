# Discovery of Many-Body and Solvent-Inclusive Collective Variables Using Permutationally Invariant Networks for Enhanced Sampling (PINES)

Implementation of PINES in PLUMED

## Requirements

* PLUMED2 (Tested using 2.8.1)
* MD Engine (Testing using GROMACS 2021.6)

## Installation

If PLUMED2+MD Engine already exist,
* Download this repository
* Execute `patch_piv_mesa_plumed.sh <path_to_PLUMED>`
* Re-install PLUMED2
  
else,
* Download MD Engine & PLUMED2.
* Download this repository
* Execute `patch_piv_mesa_plumed.sh <path_to_PLUMED>`
* Follow standard MD Engine + PLUMED2 installation process
  

## Acknowledgments

Developed by Nicholas S.M. Herringer and Siva Dasetty in Ferguson Lab.
Copyright (c) 2023.
