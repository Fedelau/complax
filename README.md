![COMPLAX](complax-high-resolution-logo.png) 
[![Version](https://img.shields.io/badge/version-1.0.0-green.svg)](https://github.com/federicalauria/complax/releases) ![PyPI](https://img.shields.io/pypi/v/complax.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) 

Complax is a Python tool designed to position solvent molecules (provided as `.xyz` files) around a user-specified atom of another molecule, also supplied in `.xyz` format.  
It was primarily developed for the microsolvation of small organic molecules, but it can also be applied to other molecular systems, for instance as a starting point for further modeling or optimization studies.

After the solvent placement, the program performs a geometry optimization using the semiempirical [xTB method](https://wires.onlinelibrary.wiley.com/doi/10.1002/wcms.1493) <sup>1</sup>.  
Optionally, the user can request an evaluation of the solvation effect in terms of potential energy differences.

## Installation

### Install via PyPI

You can install COMPLAX directly from PyPI using `pip`:

```bash
pip install complax
```

### Installation from GitHub:

Manually cloning the repository:

```bash
git clone https://github.com/Fedelau/complax.git
```
Then adding the location of the Complax directory to the PYTHONPATH environment variable.

## Requirements

- Python 3.8 or higher  
- External dependencies:
  - `numpy`
  - `ase`
  - `colorama`
  - `tqdm`
  - `tabulate`

Tested with Python 3.11.2

All dependencies can be installed with:

```bash
pip install -r requirements.txt
```
### External software 

COMPLAX interfaces with the external program [xTB](https://xtb-docs.readthedocs.io/en/latest/), which must be installed and accessible from your system’s PATH.
See https://github.com/grimme-lab/xtb for installation instructions.

xTB is developed by the Grimme group and distributed under the Creative Commons Attribution-ShareAlike 4.0 International License (CC BY-SA 4.0).

## Quick Start

Example input files are provided in the `examples/` folder.

This example shows how to position 3 tetrahydrofuran molecules (`thf.xyz`) around a lithium atom in methyllithium (`meli.xyz`).

1. First, identify the indices of the atoms to be used:
   - Lithium atom in `meli.xyz`: 5
   - Oxygen atom in `thf.xyz`: 2

2. Then run the command:
```terminal
python3 complax.py meli.xyz thf.xyz -a 5 2 -c 3
```
3. What happens next? 
   - The program first generates a file called `complax_input.xyz` containing the solute with the positioned solvent molecules.
   - Then, a geometry optimization is performed using xTB for each numbers of solvent molecules.
   - After optimization, the final geometry is saved in the `/outplax` directory named `complax_input_{n}solvent.xtbopt.xyz`, where {n} is the number of solvent molecules, choose with `-c`. 

## Commandline Usage

```terminal
python3 complax.py molecule.xyz solvent.xyz [options]
```
__Mandatory arguments:__

` -a (MOLECULE ATOM) (SOLVENT ATOM)` : Atom numbers of molecule and solvent. Format: (MOLECULE ATOM) (SOLVENT ATOM) using 1-based indexing.

` -c INT ` :  Number of file2 copies to be placed around the selected atom of file1. Default is 1.

__Optional arguments:__

` -t INT ` : Target distace from ```MOLECULE ATOM```, in Ångstrom. Default is 2.0 Ångstrom.

` --alpb SOLVENT` : Analytical linearized Poisson-Boltzman (ALPB) model, available solvents on xTB are acetone, acetonitrile, aniline, benzaldehyde, benzene, ch2cl2, chcl3, cs2, dioxane, dmf, dmso, ether, ethylacetate, furane, hexandecane, hexane, methanol, nitromethane, octanol, woctanol, phenol, toluene, thf, water. [[2]](https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00471)

` --gbsa SOLVENT` : The generalized Born solvation model (GBSA) is a simplified version of ALPB. Available solvents are acetone, acetonitrile, benzene (only GFN1-xTB), CH2Cl2, CHCl3, CS2, DMF (only GFN2-xTB), DMSO, ether, H2O, methanol, n-hexane (only GFN2-xTB), THF and toluene. [[2]](https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00471)

`-p INT` : Number of parallel processes. Default=1. During the initial optimization, the program will use the specified number of parallel processes. For subsequent optimizations (one for each solvent configuration), it will automatically launch as many parallel calculations as the number of solvent molecules selected.

`-lev LEVEL` : Level of theory for the optimization. Default is --gnf2. Other options include --gfn0, --gfn1, --gfn2, --gfnff.

`-chrg INT` : Molecular charge. Default is 0.

`-u, --uhf INT` : Number of unpaired electrons. Default is 1.


`--maxtries INT` : Maximum number of attempts to place each solvent molecule without overlaps. If the desired number of solvent molecules cannot be positioned, try increasing this value. However, if placement remains difficult, it is likely due to steric hindrance between the molecules. Default is 1000.

`--solvfx` : If specified, do a evaluation of the effect of the solvation in term of potential energy among the different systems with an increasing number of solvent molecules.

## Tips and Troubleshooting 

- The idea is to use a optimized geometry from a DFT calculations. COMPLAX keeps the internal geometry of both solute and solvent constraints, maintaining its original interatomic distances, and finding the best solvent molecules coordination.

- COMPLAX positions solvent molecules randomly around the selected atom, avoiding overlaps and taking into account only the distance between the chosen atoms. Sometimes the molecules may not be positioned exactly as expected. It is recommended to check the result visually to ensure it is correct.

- If the positioning is not satisfactory, it is worth running the program again to obtain a better configuration.

- **Important**: COMPLAX writes output files in the outplax folder. If you want to keep multiple results, move the files out of outplax before running a new calculation, otherwise the previous files will be overwritten and lost.

## Authors

**Federica Lauria**, University of Turin, Torino, Italy 


- [@Fedelau](https://github.com/Fedelau)
- [ORCID](https://orcid.org/0009-0004-0692-085X) 

**Andrea Maranzana**, University of Turin, Torino, Italy

- [ORCID](https://orcid.org/0000-0002-5524-8068)

## References

[[1]](https://wires.onlinelibrary.wiley.com/doi/10.1002/wcms.1493) C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme WIREs Comput. Mol. Sci., 2020, 11, e01493. DOI: 10.1002/wcms.1493

[[2]](https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00471) S. Ehlert, M. Stahn, S. Spicher, S. Grimme, J. Chem. Theory Comput., 2021, 17, 4250-4261 DOI: 10.1021/acs.jctc.1c00471

## License

Complax is licensed under the [MIT](https://choosealicense.com/licenses/mit/) License. 

See the [LICENSE](LICENSE) file for the full text.