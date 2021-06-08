# Change Log
------------

## v0.4.0: (2021 June)
(Version numbering system switched from two- to three-digits. All two-digit versions are older.)
- Switch QCC engine from ORCA (closed source, academic license, difficult to install using command line) to pySCF (open source, pure python API, easy to install). This was done to facilitate containerization and ease of use in python.
- Overhaul of code structure:
  - Pared down functionality: MorphCT is expected to only calculate mobility (no coarse graining/fine graining/ or device calculations).
  - The only input required is now an atomistic GSD sanpshot with lengths converted to Angstroms and the indices of the atoms in each chromophore. A helper function for determining these indices with SMARTS/SMILES grammar can be found in `morphct.chromophores.get_chromo_ids_smiles`.
  - Much of the API has been changed, but these changes are documented with docstrings for every function and an example notebook (See `examples/morphct-workflow.ipynb`).

## v3.1: 
- Many quality-of-life improvements, issue resolutions, and bug fixes including but not limited to: the ability to write ORCA jobs to a RAM disk for 30x faster calculations, significantly more unit tests to increase code coverage to 75%, full testing the device simulations, better python formatting, and migration from Coveralls to Codecov.io.
## v3.0: 
- MorphCT codebase brought up to PEP8 standard (with line limit of 120 characters), refactored to work as a package, added extensive unit tests to the pipeline (pytest) and added continuous integration support using Shippable and Coveralls.
## v2.2: 
- Additional funcionality added for blend morphologies, multiple donor/acceptor species, variable-range hopping, Voronoi neighbourhood analysis and other performance tweaks. Results in preparation to be submitted in Q2 2018.
## v2.1: 
- MorphCT updated from python 2.7 to python 3.5
## v2.0: 
- Hardcode removed, utility of MorphCT expanded to any organic molecule with more customizable features (support for small molecules included)
## v1.0: 
- MorphCT released under GPL. Hardcoded mobility results for P3HT, results published in [Molecular Simulation](https://doi.org/10.1080/08927022.2017.1296958) 
