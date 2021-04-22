[<img src=".github/logo.png" width="500">](https://github.com/cmelab/morphct)

[![DOI](https://zenodo.org/badge/100152281.svg)](https://zenodo.org/badge/latestdoi/100152281)
[![build](https://github.com/cmelab/morphct/actions/workflows/build.yml/badge.svg)](https://github.com/cmelab/morphct/actions/workflows/build.yml)
[![pytest](https://github.com/cmelab/morphct/actions/workflows/pytest.yml/badge.svg)](https://github.com/cmelab/morphct/actions/workflows/pytest.yml)
[![codecov](https://codecov.io/gh/cmelab/morphct/branch/master/graph/badge.svg?token=PhAfcr15av)](https://codecov.io/gh/cmelab/morphct)

# Intellectual Property #

MorphCT has been released under a GPL3 license (see LICENSE). Please read and cite our [Molecular Simulations paper](https://doi.org/10.1080/08927022.2017.1296958) and our [zenodo DOI](https://zenodo.org/badge/latestdoi/100152281) if you use our code in your published work.

---

# Code Description #

The intention of this code is to form a modular, computational pipeline that can powerfully and coherently relate organic molecular morphology on the Angstrom lengthscale to electronic device performance over hundreds of nanometers.

MorphCT accomplishes this by:

* Splitting a morphology into electronically-active chromophores that charge carriers (holes in the case of donor materials, electrons for acceptors) are likely to be delocalized along, and can perform quantized charge hops between.
* Performing high-throughput, fast quantum chemical calculations (QCCs) to obtain the energetic landscape caused by conformational disorder, as well as electronic transfer integrals between chromophore pairs.
* Using these calculated electronic properties as inputs into a kinetic Monte Carlo (KMC) algorithm to simulate the motion of charge carriers throughout the device, allowing carrier mobilities to be obtained (a good proxy for device performance)

---

# Getting Started #

## Installation ##

### Using a container
To use MorphCT in a prebuilt container (using [Singularity](https://singularity.lbl.gov/)), run:
```bash
singularity pull docker://cmelab/morphct:latest
singularity exec morphct_latest.sif bash
```

**Or** using [Docker](https://docs.docker.com/), run:
```bash
docker pull cmelab/morphct:latest
docker run -it cmelab/morphct:latest
```

### Custom install

To create a local environment with [conda](https://docs.conda.io/en/latest/miniconda.html), run:
```bash
conda env create -f environment.yml
conda activate morphct
```
And to test your installation, run:
```
pytest
```

---

# Running Jobs #

## Inputs ##

* After setup, MorphCT requires the input molecular system (atomistic), and atomic indices of the chromophore (can be obtained via SMARTS matching).
* An example is provided for the user in `morphct/examples/morphct-workflow.ipynb`.
