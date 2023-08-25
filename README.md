# Permutationally Invariant Networks for Enhanced Sampling (PINES): Discovery of Multi-Molecular and Solvent-Inclusive Collective Variables


Implementation of PINES in PLUMED2

## Requirements

* PLUMED2
* MD Engine (GROMACS)

## Installation

Case 1: With installed PLUMED2+MD,
  
-  Follow this for reinstalling PLUMED with PINES:
   * Download this repository
   * Execute `patch_piv_mesa_plumed.sh <path_to_PLUMED>`
   * Re-install PLUMED2 with `--enable-modules=all` or `--enable-modules=pines+annbfunc`
- Alternative usage without reinstalling PLUMED:
  * use `LOAD` command in PLUMED to link in plumed.dat file with the `PINES` codes in this repository. Example usage can be found in [examples/driverTests](examples/driverTests).

Case 2: New installation of PLUMED2+MD,
* Download this repository
* Download MD Engine & PLUMED2.
* Execute `patch_piv_mesa_plumed.sh <path_to_PLUMED>`
* Follow standard MD Engine + PLUMED2 installation process with `--enable-modules=all` or `--enable-modules=pines+annbfunc`

  
## Examples

- [Argon 13 cluster](examples/SimulationFiles/Ar13)
- [C45 in water](examples/SimulationFiles/C45/)
- [NaCl in water](examples/SimulationFiles/NaCl)

## Parameters

Please refer [parameters.md](parameters.md) for a template PLUMED file along with annotations explaining each parameter. In case the annotations are not clear, please let us know by opening an issue in this repo.

## Citation

If you use this code in your work, please cite:

Herringer, Nicholas SM; Dasetty, Siva; Gandhi, Diya; Lee, Junhee; Ferguson, Andrew L. "Permutationally Invariant Networks for Enhanced Sampling (PINES): Discovery of Multi-Molecular and Solvent-Inclusive Collective Variables." arXiv preprint arXiv:2308.08680. DOI: [https://doi.org/10.48550/arXiv.2308.08680](https://doi.org/10.48550/arXiv.2308.08680)

```
@misc{herringer2023permutationally,
      title={Permutationally Invariant Networks for Enhanced Sampling (PINES): Discovery of Multi-Molecular and Solvent-Inclusive Collective Variables}, 
      author={Nicholas S. M. Herringer and Siva Dasetty and Diya Gandhi and Junhee Lee and Andrew L. Ferguson},
      year={2023},
      eprint={2308.08680},
      archivePrefix={arXiv},
      primaryClass={q-bio.BM}
}
```
