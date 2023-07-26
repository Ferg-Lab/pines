# Discovery of Many-Body and Solvent-Inclusive Collective Variables Using Permutationally Invariant Networks for Enhanced Sampling (PINES)

Implementation of PINES in PLUMED2

## Requirements

* PLUMED2
* MD Engine (GROMACS)

## Installation

If PLUMED2+MD engine exist,
  
-  Follow this for reinstalling PLUMED with PINES:
   * Download this repository
   * Execute `patch_piv_mesa_plumed.sh <path_to_PLUMED>`
   * Re-install PLUMED2 with `--enable-modules=all` or `--enable-modules=pines+annbfunc`
- Alternative usage without reinstalling PLUMED:
  * use `LOAD` command in PLUMED to link in plumed.dat file with the `PINES` codes in this repository. Example usage can be found in [examples/driverTests](examples/driverTests).

else,
* Download this repository
* Download MD Engine & PLUMED2.
* Execute `patch_piv_mesa_plumed.sh <path_to_PLUMED>`
* Follow standard MD Engine + PLUMED2 installation process with `--enable-modules=all` or `--enable-modules=pines+annbfunc`

  
## Examples

- [Argon 13 cluster](examples/SimulationFiles/Ar13)
- [C45 in water](examples/SimulationFiles/C45/)
- [NaCl in water](examples/SimulationFiles/NaCl)

## Parameters

Please refer [parameters.md](parameters.md) for a template PLUMED file along with annotations explaining each parameter. In case if the annotations are not helpful, please open an issue and explain your problem.
