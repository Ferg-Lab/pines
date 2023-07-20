# Discovery of Many-Body and Solvent-Inclusive Collective Variables Using Permutationally Invariant Networks for Enhanced Sampling (PINES)

Implementation of PINES in PLUMED2

## Requirements

* PLUMED2
* MD Engine (GROMACS)

## Installation

If PLUMED2+MD Engine already exist,
* Download this repository
* Execute `patch_piv_mesa_plumed.sh <path_to_PLUMED>`
* Re-install PLUMED2 with `--enable-modules=all` or `--enable-modules=pines`

  
else,
* Download MD Engine & PLUMED2.
* Download this repository
* Execute `patch_piv_mesa_plumed.sh <path_to_PLUMED>`
* Follow standard MD Engine + PLUMED2 installation process with `--enable-modules=all` or `--enable-modules=pines`
  
## Examples

- [Argon 13 cluster](examples/SimulationFiles/Ar13)
- [C45 in water](examples/SimulationFiles/C45/)
- [NaCl in water](examples/SimulationFiles/NaCl)

## Parameters

Please refer [parameters.md](parameters.md) for a template PLUMED file along with annotations explaining each parameter. In case if the annotations are not helpful, please open an issue and explain your problem.
