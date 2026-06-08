![COMPLAX](complax-high-resolution-logo.png) 
[![Python package](https://github.com/Fedelau/complax/actions/workflows/python-package.yml/badge.svg)](https://github.com/Fedelau/complax/actions/workflows/python-package.yml) ![PyPI](https://img.shields.io/pypi/v/complax.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) 

COMPLAX is an automated Python workflow designed to position solvent molecules (provided as `.xyz` files) around a user-specified atom of another molecule, also supplied in `.xyz` format.  
It was primarily developed for the microsolvation of small organic molecules, but it can also be applied to other molecular systems, serving as a robust starting point for high-level quantum mechanical modeling (e.g., DFT).

Following an unbiased stochastic solvent placement, the program performs two-stage geometry optimization using the semiempirical [xTB](https://xtb-docs.readthedocs.io/en/latest/) methods. [[1]](https://wires.onlinelibrary.wiley.com/doi/10.1002/wcms.1493)

A constrained pre-relaxation prevents unphysical solvent detachment, which is then automatically followed by a full structural equilibration. 

Optionally, the user can perform a stochastic sampling to explore the conformational space and request an on-the-fly evaluation of the solvation effect in terms of potential energy differences.

## Installation

### Install via PyPI

You can install COMPLAX directly from PyPI using `pip`:

```bash
pip install complax
```

### Installation from GitHub:

Manually clone the repository and install it via `pip`:

```bash
git clone [https://github.com/Fedelau/complax.git](https://github.com/Fedelau/complax.git)
cd complax
pip install .
```

## Requirements

- Python 3.8 or higher  
- External dependencies (automatically installed via `pip`):
  - `numpy`
  - `ase`
  - `colorama`
  - `tqdm`
  - `tabulate`

*Tested with Python 3.11.2*

### External Software 

COMPLAX interfaces with the external program xTB, which must be installed and accessible from your system’s PATH.
See https://github.com/grimme-lab/xtb for installation instructions. Since this version of COMPLAX utilizes the recent g-xTB semi-empirical electronic structure method, to run a calculation with this level of theory, the modified xTB binaries are required, available at this link https://github.com/grimme-lab/g-xtb. [[2]](https://doi.org/10.26434/chemrxiv-2025-bjxvt)

xTB is developed by the Grimme group and distributed under the Creative Commons Attribution-ShareAlike 4.0 International License (CC BY-SA 4.0).

## Quick Start

Example input files are provided in the `examples/` folder. (If you installed `complax` via `pip`, you can access these files by cloning the repository or downloading them directly from the [GitHub `examples` directory](https://github.com/Fedelau/complax/tree/main/examples)).

This example shows how to position 3 tetrahydrofuran molecules (`thf.xyz`) around a lithium atom in methyllithium (`meli.xyz`).

1. First, identify the indices of the atoms to be used:
   - Lithium atom in `meli.xyz`: 5
   - Oxygen atom in `thf.xyz`: 2

2. Then run the command:
```terminal
complax meli.xyz thf.xyz -a 5 2 -c 3
```
3. What happens next? 
   - The program places the molecules avoiding steric clashes.
   - Then, a robust two-stage geometry optimization is performed using xTB for each incremental number of solvent molecules. During the entire process, the solute coordinates and the internal geometries of the solvent molecules (intramolecular distances) are kept completely frozen.
     - Stage 1 (Pre-relaxation): The intermolecular target distance (e.g., Li-O) is harmonically constrained for 50 cycles. This crucial step dissipates initial spatial clashes and prevents the solvent molecules from artificially detaching or flying away.
     - Stage 2 (Full equilibration): The intermolecular distance constraint is released. The rigid solvent molecules are now completely free to translate and reorient around the fixed solute core to find their global thermodynamic minimum.
   - After optimization, the final geometry is saved in the `outplax/` directory named `complax_struct_{n}solvent_{s}.xtbopt.xyz`, where `{n}` is the number of solvent molecules and `{s}` is the structure index.

### Conformational Sampling and Energy Evaluation Example

To get the most out of COMPLAX, you can combine multiple flags to run a fully automated conformational search and thermodynamic analysis in a single line. 

This command generates 40 different starting configurations (`--nstruct 40`), retains only the top 3 most stable geometries for each solvation state (`--keep-best 3`), and computes the stabilization energies (`--solvfx`):

```terminal
complax meli.xyz thf.xyz -a 5 2 -c 3 --nstruct 40 --keep-best 3 --solvfx
```

What happens here?

- COMPLAX explores the coordination space by optimizing 40 random clusters.
- It automatically ranks them by total energy, discarding the highest-energy structures to save disk space and keeping only the top 3 conformers per state.
- Finally, it prints a clean summary report on the terminal (and saves it in `complax_summary.txt`), displaying both the energy ranking and the stepwise solvation effect:

```terminal
[...]
--- 3 solvent molecule(s) ---
complax_struct_3solvent_7.xtbopt.xyz          -103.573414 Eh
complax_struct_3solvent_23.xtbopt.xyz         -103.573345 Eh
complax_struct_3solvent_39.xtbopt.xyz         -103.573296 Eh
[...]

                     ***************************
                     * Effect of the Solvation *
                     ***************************

+--------------------------------+------------+------------+------------+
|                                |          1 |          2 |          3 |
+================================+============+============+============+
| Sum of reag + solv (Hartree)   | -20.609217 | -37.322501 | -54.035785 |
+--------------------------------+------------+------------+------------+
| Tot Energy (Complex) (Hartree) | -20.631873 | -37.364828 | -54.095412 |
+--------------------------------+------------+------------+------------+
| ΔE (kcal/mol)                  |     -14.22 |     -26.56 |     -37.42 |
+--------------------------------+------------+------------+------------+
```



## Commandline Usage

```terminal
complax molecule.xyz solvent.xyz [options]
```

__Mandatory arguments:__

| Argument | Description | Default |
| :--- | :--- | :--- |
| `-a (MOL) (SOLV)` | Atom numbers of molecule and solvent. Format: `(MOLECULE ATOM) (SOLVENT ATOM)` using 1-based indexing. | *Required* |
| `-c INT` | Number of solvent copies (`solvent.xyz`) to be placed around the selected atom of the solute. | `1` |

<br>

__Optional arguments:__

| Argument | Description | Default |
| :--- | :--- | :--- |
| `-t FLOAT` | Target distance from the selected molecule atom, in Ångstrom. | `2.0` |
| `--alpb SOLVENT` | Analytical linearized Poisson-Boltzman (ALPB) model. Available solvents on xTB are acetone, acetonitrile, aniline, benzaldehyde, benzene, ch2cl2, chcl3, cs2, dioxane, dmf, dmso, ether, ethylacetate, furane, hexandecane, hexane, methanol, nitromethane, octanol, woctanol, phenol, toluene, thf, water. [[3]](https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00471) | *None* |
| `--gbsa SOLVENT` | Generalized Born solvation model (GBSA) is a simplified version of ALPB. Available solvents are acetone, acetonitrile, benzene (only GFN1-xTB), CH2Cl2, CHCl3, CS2, DMF (only GFN2-xTB), DMSO, ether, H2O, methanol, n-hexane (only GFN2-xTB), THF and toluene. [[3]](https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00471) | *None* |
| `--gbe SOLVENT` | Generalized Born model with finite epsilon for solvation. Requires `--lev gxtb`. | *None* |
| `-p INT` | Number of parallel processes. During the initial optimization, the program uses the specified number of processes. For subsequent optimizations, it automatically launches as many parallel calculations as the number of solvent molecules selected. | `1` |
| `--lev LEVEL` | Level of theory for the optimization. Options include `gfn0`, `gfn1`, `gfn2`, `gfnff`, `gxtb`. | `gfn2` |
| `--chrg INT` | Molecular charge. | `0` |
| `-u`, `--uhf INT` | Number of unpaired electrons. | `0` |
| `--res INT` | Number of points in the Fibonacci spherical grid. | `10000` |
| `--nstruct INT` | Number of stochastic, non-overlapping starting structures to generate and optimize. Ideal for exploring the conformational space of the solvated complex. | *None* |
| `--keep-best INT` | Filters the output to keep only the `N` lowest-energy optimized structures for each solvent configuration, automatically deleting the less stable ones. | *None* |
| `--solvfx` | Evaluates the effect of solvation in terms of potential energy differences among systems with an increasing number of solvent molecules. | *False* |
| `--cutoff FLOAT` | Steric overlap distance threshold in Ångstrom. | `1.6` |
| `--seed INT` | Seed for the pseudo-random number generator to ensure reproducibility. Set to `-1` for pure stochasticity. | `42` |

## Tips and Troubleshooting 

- **Stochastic Sampling**: For complex solvent environments, use `--nstruct 10 --keep-best 3`. COMPLAX will generate 10 different random orientations of the solvent cluster, optimize all of them, and retain only the 3 most thermodynamically stable configurations, saving you a massive amount of manual sorting.

- **Pre-optimized geometries**: The ideal workflow uses pre-optimized input geometries from prior DFT calculations. COMPLAX keeps the internal geometry of both solute and solvent rigidly constrained during the entire xTB optimization process. This maintains their original intramolecular distances while exclusively allowing the exploration of the best intermolecular solvent coordination.

- **Visual Check**: COMPLAX stochastically positions solvent molecules around the selected atom, avoiding steric overlaps based on the chosen target distance. Due to the unbiased nature of the sampling, especially in highly congested spatial environments, it is always recommended to visually inspect the final optimized geometries to ensure the resulting coordination aligns with your chemical intuition.

- **Important**: COMPLAX writes output files in the `outplax/` folder. If you want to keep multiple results, move the files out of outplax before running a new calculation or run complax in another folder, otherwise the previous files will be overwritten and lost.

## Authors

**Federica Lauria**, Department of Chemistry, University of Turin, Torino, Italy 


- [@Fedelau](https://github.com/Fedelau)
- [ORCID](https://orcid.org/0009-0004-0692-085X) 

**Andrea Maranzana**, Deparment of Chemistry, University of Turin, Torino, Italy

- [ORCID](https://orcid.org/0000-0002-5524-8068)

## References

[[1]](https://wires.onlinelibrary.wiley.com/doi/10.1002/wcms.1493) C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme WIREs Comput. Mol. Sci., 2020, 11, e01493. DOI: 10.1002/wcms.1493

[[2]](https://doi.org/10.26434/chemrxiv-2025-bjxvt) T. Froitzheim, M. Müller, A. Hansen, et al., ChemRxiv. 23 June 2025. DOI: 10.26434/chemrxiv-2025-bjxvt

[[3]](https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00471) S. Ehlert, M. Stahn, S. Spicher, S. Grimme, J. Chem. Theory Comput., 2021, 17, 4250-4261 DOI: 10.1021/acs.jctc.1c00471

## License

Complax is licensed under the [MIT](https://choosealicense.com/licenses/mit/) License. 

See the [LICENSE](LICENSE) file for the full text.

## Citation

If you use COMPLAX in your research, please cite our upcoming paper:
*(Placeholder for the journal reference - currently under review).*
