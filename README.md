<a href="https://stsci.edu">
  <img src="docs/_static/stsci_logo.png" alt="STScI Logo" width="15%" style="margin-left: auto;"/>
  <img src="docs/_static/stsci_name.png" alt="STScI Logo" width="68%"/>
</a>
<a href="https://science.nasa.gov/mission/webb/">
  <img src="docs/_static/jwst_logo.png" alt="JWST Logo" width="15%" style="margin-right: auto;"/>
</a>

# James Webb Space Telescope Calibration Pipeline

[![DOI](https://zenodo.org/badge/60551519.svg)](https://zenodo.org/badge/latestdoi/60551519)
[![PyPI](https://img.shields.io/pypi/v/jwst.svg)](https://pypi.org/project/jwst)
[![Python Support](https://img.shields.io/pypi/pyversions/jwst)](https://pypi.org/project/jwst/)
[![Powered by STScI](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](https://www.stsci.edu)
[![Powered by Astropy](https://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](https://www.astropy.org/)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![build](https://github.com/spacetelescope/jwst/actions/workflows/build.yml/badge.svg)](https://github.com/spacetelescope/jwst/actions/workflows/build.yml)
[![tests](https://github.com/spacetelescope/jwst/actions/workflows/tests.yml/badge.svg)](https://github.com/spacetelescope/jwst/actions/workflows/tests.yml)
[![readthedocs](https://readthedocs.org/projects/jwst-pipeline/badge/?version=latest)](https://jwst-pipeline.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/spacetelescope/jwst/branch/main/graph/badge.svg?token=Utf5Zs9g7z)](https://codecov.io/gh/spacetelescope/jwst)

This package (`jwst`) processes uncalibrated data from both imagers and spectrographs onboard the [James Webb Space Telescope (JWST)](https://science.nasa.gov/mission/webb/), an orbiting infrared observatory stationed at Earth-Sun L<sub>2</sub>.
The pipeline performs a series of calibration steps that result in standard data products usable for science.

Detailed explanations of specific calibration stages, reference files, and pipeline builds can be found on the [ReadTheDocs pages](https://jwst-pipeline.readthedocs.io) and [JDox](https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline).

> [!NOTE]
> If you have trouble installing this package, have encountered a bug while running the pipeline, or wish to request a new feature,
> please [open an issue on GitHub](https://github.com/spacetelescope/jwst/issues) or [contact the JWST Help Desk](https://jwsthelp.stsci.edu).

<!--toc:start-->
- [Quick Start](#quick-start)
  - [1. Install the Pipeline](#1-install-the-pipeline)
    - [Option: Build Pipeline Directly from Source Code](#option-build-pipeline-directly-from-source-code)
    - [Option: Install Exact Operational Environment](#option-install-exact-operational-environment)
  - [2. Set up the Calibration Reference Data System (CRDS)](#2-set-up-the-calibration-reference-data-system-crds)
  - [3. Run the Pipeline](#3-run-the-pipeline)
- [Code Contributions](#code-contributions)
- [DMS Operational Build Versions](#dms-operational-build-versions)
<!--toc:end-->

## Quick Start

### 1. Install the Pipeline

> [!IMPORTANT]
> The JWST calibration pipeline currently supports Linux and macOS.
> Native Windows builds are **not** currently supported; [use WSL instead](https://stenv.readthedocs.io/en/latest/windows.html).

We recommend using an isolated Python environment to install `jwst`.
Python "environments" are isolated Python installations, confined to a single directory, where you can install packages, dependencies, and tools without cluttering your system Python libraries.
You can manage environments with `mamba` / `conda`, `virtualenv`, `uv`, etc.

These instructions assume you are creating Conda environments with the `mamba` command
(see [Miniforge for installation instructions](https://github.com/conda-forge/miniforge/blob/main/README.md));
to use `conda` instead, simply replace `mamba` with `conda` in the following commands.

First, create an empty environment with Python installed:
```shell
mamba create -n jwst_env python=3.13
```

Then, **activate** that environment (necessary to be able to access this isolated Python installation):
```shell
mamba activate jwst_env
```

Finally, install `jwst` into the environment:
```shell
pip install jwst
```

Without a specified version, `pip` defaults to the latest released version that supports your environment.
To install a specific version of `jwst`, explicitly set that version in your `pip install` command:
```shell
pip install jwst==1.20.2
```

To install a different version of `jwst`, simply create a new environment for that version:
```shell
mamba create -n jwst1.20_env python=3.13
mamba activate jwst1.20_env
pip install jwst==1.20

mamba create -n jwst1.19_env python=3.12
mamba activate jwst1.19_env
pip install jwst==1.19
```

#### Option: Build Pipeline Directly from Source Code

> [!IMPORTANT]
> You need a C compiler to build the JWST calibration pipeline (and dependencies) from source.

To install the latest unreleased (and unstable) development version directly from the source code on GitHub:

```shell
pip install git+https://github.com/spacetelescope/jwst
```

#### Option: Install Exact Operational Environment

There may be occasions where you need to replicate the exact environment used for canonical calibration operations by STScI (e.g. for validation testing or debugging issues).
We package releases for operations [as environment snapshots that specify exact versions for both the pipeline and all dependencies](https://ssb.stsci.edu/stasis/releases/jwst/).

See the [DMS Operational Build Versions](#dms-operational-build-versions) table for the version of the pipeline corresponding to each operational build.
For example, use `jwst==1.17.1` for **DMS build 11.2**.
Also note that Linux and macOS systems require different snapshot files:
```shell
mamba env create --file https://ssb.stsci.edu/stasis/releases/jwst/JWSTDP-1.18.1/delivery/latest-py312-macos-arm64.yml
mamba activate JWSTDP-1.18.1-1-py312-macos-arm64
```

> [!NOTE]
> Starting with `jwst==1.16.1`, the JWST pipeline uses [`stasis`](https://github.com/spacetelescope/stasis) to package environments and deliver releases. If you need a version of `jwst` prior to `1.16.1`, use a slightly different procedure:
> ```shell
> mamba create -n jwstdp-1.16.0 --file https://ssb.stsci.edu/releases/jwstdp/1.16.0/conda_python_macos-stable-deps.txt
> mamba activate jwstdp-1.16.0
> pip install -r https://ssb.stsci.edu/releases/jwstdp/1.16.0/reqs_macos-stable-deps.txt
> ```

### 2. Set up the Calibration Reference Data System (CRDS)

Before running the pipeline, you must first set up your local machine to retrieve files from the [Calibration Reference Data System (CRDS)](https://jwst-crds.stsci.edu/static/users_guide/index.html>)
CRDS provides calibration reference files for several telescopes, including JWST.

Set `CRDS_SERVER_URL` and `CRDS_PATH` to run the pipeline with access to reference files from CRDS:
```shell
export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
export CRDS_PATH=$HOME/data/crds_cache/
```

> [!NOTE]
> The CRDS PUB Server (`https://jwst-crds-pub.stsci.edu`) was decommissioned in March 2023.
> To use historical files from the PUB server, [contact the JWST Help Desk](https://jwsthelp.stsci.edu).

The pipeline will automatically download individual reference files and cache them in the `CRDS_PATH` directory.
**Expect to use 50 gigabytes (or more) of disk space for reference files**, depending on the instrument modes in use.

> [!TIP]
> Users within the STScI network do not need to set `CRDS_PATH` (it defaults to shared network storage).

To use a specific CRDS context other than that [automatically associated with a given pipeline version](https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline/crds-migration-to-quarterly-calibration-updates), explicitly set the `CRDS_CONTEXT` environment variable:
```shell
export CRDS_CONTEXT=jwst_1179.pmap
```

For more information, see [the docs page on JWST CRDS reference files](https://jwst-pipeline.readthedocs.io/en/stable/jwst/user_documentation/reference_files_crds.html#reference-files-crds).

### 3. Run the Pipeline

Once installed, the pipeline allows users to run and configure calibration themselves for custom processing of JWST data,
either [from the command line with `strun`](https://jwst-pipeline.readthedocs.io/en/stable/jwst/user_documentation/running_pipeline_command_line.html)
or from Python with [pipeline and step functions and classes in the `jwst` package](https://jwst-pipeline.readthedocs.io/en/stable/jwst/user_documentation/running_pipeline_python.html)
(see [this curated set of Jupyter notebooks](https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline/jwst-pipeline-notebooks) for example usage).
Additionally, the `jwst` package provides [JWST datamodel classes](https://jwst-pipeline.readthedocs.io/en/stable/jwst/user_documentation/datamodels.html),
the recommended method for reading and writing JWST data files in Python.

## Code Contributions

`jwst` is an open source package written in Python.
The source code is [available on GitHub](https://github.com/spacetelescope/jwst).
New contributions and contributors are very welcome!
Please read [`CONTRIBUTING.md`](CONTRIBUTING.md),
the [public API definition](https://jwst.readthedocs.io/en/latest/jwst/user_documentation/more_information.html#api-public-vs-private),
and the [public API deprecation policy](https://jwst-pipeline.readthedocs.io/en/latest/jwst/user_documentation/more_information.html#api-deprecation-policy)

We strive to provide a welcoming community by abiding with our [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md).

See [`TESTING.md`](./TESTING.md) for instructions on automated testing.

## DMS Operational Build Versions

The table below provides information on each release of the `jwst` package and its relationship to software builds used in STScI JWST DMS operations.
Each `jwst` tag was released on PyPI on the date given in `Released`, and then subsequently installed into operations on the date given in `Ops Install`.

| jwst tag            | DMS build | SDP_VER  | CRDS_CONTEXT | Released   | Ops Install | Notes                                         |
|---------------------|-----------|----------|--------------|------------|-------------|-----------------------------------------------|
| 1.20.2              | B12.1.1   | 2025.4.1 | 1464         | 2025-10-31 |             | Patch release for B12.1.1                     |
| 1.20.1              | B12.1     | 2025.4.0 | 1464         | 2025-10-20 |             | Patch release for B12.1                       |
| 1.20.0              | B12.1     | 2025.4.0 | 1462         | 2025-10-15 |             | First release candidate for B12.1             |
| 1.19.2              | B12.0.2   | 2025.3.0 | 1408         | 2025-09-11 | 2025-10-06  | Patch release for B12.0.2                     |
| 1.19.1              | B12.0.1   | 2025.3.0 | 1408         | 2025-07-21 | 2025-08-26  | Patch release for B12.0.1                     |
| 1.19.0              | B12.0     | 2025.3.0 | 1408         | 2025-06-26 |             | First release candidate for B12.0             |
| 1.18.1              | B11.3.1   | 2025.2.1 | 1364         | 2025-06-10 |             | Patch release for B11.3.1                     |
| 1.18.0              | B11.3     | 2025.2.0 | 1364         | 2025-04-01 | 2025-05-20  | First release for B11.3                       |
| 1.17.1              | B11.2     | 2025.1.0 | 1321         | 2025-01-02 | 2025-03-05  | Final release candidate for B11.2             |
| 1.17.0              | B11.2     | 2025.1.0 | 1321         | 2024-12-20 |             | First release candidate for B11.2             |
| 1.16.1              | B11.1.1   | 2024.3.1 | 1303         | 2024-11-13 | 2024-12-06  | Final release candidate for B11.1             |
| 1.16.0              | B11.1     | 2024.3.0 | 1298         | 2024-09-20 |             | First release candidate for B11.1             |
| 1.15.1              | B11.0     | 2024.2.2 | 1293         | 2024-07-08 | 2024-09-12  | Final release candidate for B11.0             |
| 1.15.0              | B11.0rc1  |          | 1274         | 2024-06-26 |             | First release candidate for B11.0             |
| 1.14.1              |           |          | 1240         | 2024-06-27 |             | PyPI-only release for external users          |
| 1.14.0              | B10.2.1   | 2024.1.1 | 1240         | 2024-03-29 | 2024-06-12  | Final release candidate for B10.2.1           |
| 1.13.4              |           |          | 1210         | 2024-01-25 |             | PyPI-only release for external users          |
| 1.13.3              | B10.1     | 2023.4.0 | 1181         | 2024-01-05 |             | Final release candidate for B10.1             |
| 1.13.2              | B10.1rc3  | 2023.4.0 | 1181         | 2023-12-21 |             | Third release candidate for B10.1             |
| 1.13.1              | B10.1rc2  | 2023.4.0 | 1181         | 2023-12-19 |             | Second release candidate for B10.1            |
| 1.13.0              | B10.1rc1  | 2023.4.0 | 1179         | 2023-12-15 |             | First release candidate for B10.1             |
| 1.12.5              | B10.0.1   | 2023.3.1 | 1166         | 2023-10-19 | 2023-12-05  | Patch release B10.0.1                         |
| 1.12.4              |           |          |              | 2023-10-12 |             | Pinning dependencies for external users       |
| 1.12.3              | B10.0     | 2023.3.0 | 1135         | 2023-10-03 | 2023-12-05  | Final release candidate for B10.0             |
| 1.12.2              | B10.0rc3  |          | 1135         | 2023-10-02 |             | Third release candidate for B10.0             |
| 1.12.1              | B10.0rc2  |          | 1132         | 2023-09-26 |             | Second release candidate for B10.0            |
| 1.12.0              | B10.0rc1  |          | 1130         | 2023-09-18 |             | First release candidate for B10.0             |
| 1.11.4              | B9.3.1    | 2023.2.1 | 1107         | 2023-08-14 | 2023-08-24  | Final release for B9.3.1 patch                |
| 1.11.3              | B9.3      | 2023.2.0 | 1097         | 2023-07-17 |             | Final release candidate for B9.3              |
| 1.11.2              | B9.3rc3   |          | 1097         | 2023-07-12 |             | Third release candidate for B9.3              |
| 1.11.1              | B9.3rc2   |          | 1094         | 2023-06-29 |             | Second release candidate for B9.3             |
| 1.11.0              | B9.3rc1   |          | 1094         | 2023-06-21 |             | First release candidate for B9.3              |
| 1.10.2              |           |          | 1077         | 2023-04-14 |             | Pinning dependencies for external users       |
| 1.10.1              | B9.2.x    | 2023.1.1 | 1077         | 2023-04-13 | 2023-05-23  | Final release candidate for B9.2              |
| 1.10.0              | B9.2rc1   |          | 1075         | 2023-03-31 |             | First release candidate for B9.2              |
| 1.9.6               | B9.1.2    | 2022.5.2 | 1068         | 2023-03-09 | 2023-03-15  | Final release candidate for B9.1.2            |
| 1.9.5               |           |          | 1061         | 2023-03-02 |             | First release candidate for B9.1.2            |
| 1.9.4               | B9.1.1    | 2022.5.1 | 1041         | 2023-01-27 | 2023-02-28  | Final release candidate for B9.1.1            |
| 1.9.3               | B9.1      | 2022.5.0 | 1030         | 2023-01-12 | 2023-02-28  | Final release candidate for B9.1              |
| 1.9.2               | B9.1rc2   |          |              | 2023-01-04 |             | Second release candidate for B9.1 (hotfix)    |
| 1.9.1               | B9.1rc2   |          |              | 2023-01-03 |             | Second release candidate for B9.1             |
| 1.9.0               | B9.1rc1   |          |              | 2022-12-27 |             | First release candidate for B9.1              |
| 1.8.5               | B9.0      |          | 1019         | 2022-12-12 |             | Documentation patch release for B9.0          |
| 1.8.4               | B9.0      |          |              | 2022-11-16 |             | Documentation patch release for B9.0          |
| 1.8.3               | B9.0      |          |              | 2022-11-11 |             | Documentation patch release for B9.0          |
| 1.8.2               | B9.0      | 2022.4.0 | 1017         | 2022-10-19 | 2022-11-17  | Final release candidate for B9.0              |
| 1.8.1               | B9.0rc2   |          |              | 2022-10-17 |             | Second release candidate for B9.0             |
| 1.8.0               | B9.0rc1   |          |              | 2022-10-10 |             | First release candidate for B9.0              |
| 1.7.2               | B8.1.2    | 2022.3.1 | 0984         | 2022-09-12 | 2022-09-21  | Final release candidate for B8.1.2            |
| 1.7.1               | B8.1.2rc2 |          |              | 2022-09-07 |             | Second release candidate for B8.1.2           |
| 1.7.0               | B8.1.2rc1 |          |              | 2022-09-01 |             | First release candidate for B8.1.2            |
| 1.6.2               | B8.1      | 2022.3.0 | 0953         | 2022-07-19 | 2022-08-19  | Final release candidate for B8.1              |
| 1.6.1               | B8.1rc2   |          |              | 2022-07-15 |             | Second release candidate for B8.1             |
| 1.6.0               | B8.1rc1   |          |              | 2022-07-11 |             | First release candidate for B8.1              |
| 1.5.3               | B8.0.1    | 2022.2.1 | 0913         | 2022-06-20 | 2022-06-30  | Patch release B8.0.1                          |
| 1.5.2               | B8.0      | 2022.2.0 | 0874         | 2022-05-20 | 2022-06-16  | Final release candidate for B8.0              |
| 1.5.1               | B8.0rc2   |          |              | 2022-05-17 |             | Second release candidate for B8.0             |
| 1.5.0               | B8.0rc1   |          |              | 2022-05-05 |             | First release candidate for B8.0              |
| 1.4.6               | B7.9.3    | 2022.1.2 | 0800         | 2022-03-25 |             | Final release candidate for B7.9.3            |
| 1.4.5               | B7.9.3rc2 |          |              | 2022-03-23 |             | Second release candidate for B7.9.3           |
| 1.4.4               | B7.9.3rc1 |          |              | 2022-03-16 |             | First release candidate for B7.9.3            |
| 1.4.3               | B7.9.1    | 2022.1.1 | 0800         | 2022-02-03 |             | Final B7.9.1                                  |
| 1.4.2               | B7.9      | 2022.1.0 | 0797         | 2022-01-20 |             | Final release candidate for B7.9              |
| 1.4.1               | B7.9rc2   |          |              | 2022-01-15 |             | Second release candidate for B7.9             |
| 1.4.0               | B7.9rc1   |          |              | 2022-01-10 |             | First release candidate for B7.9              |
| Pre-launch releases |           |          |              |            |             |                                               |
| 1.3.3               | B7.8.2    | 2021.4.0 | 0764         | 2021-10-05 |             | Same as 1.3.2, but with installation bug fix  |
| 1.3.2               | B7.8.2    | 2021.4.0 | 0764         | 2021-09-03 |             | Final release candidate for B7.8.2            |
| 1.3.1               | B7.8.1    | 2021.3.0 | 0742         | 2021-08-09 |             | Final release candidate for B7.8.1            |
| 1.3.0               | B7.8.1rc1 |          | 0741         | 2021-08-02 |             | First release candidate for B7.8.1            |
| 1.2.3               | B7.8      | 2021.2.0 | 0732         | 2021-06-08 |             | Final release candidate for B7.8              |
| 1.2.2               | B7.8rc3   |          |              | 2021-06-08 |             | Third release candidate for B7.8              |
| 1.2.1               | B7.8rc2   |          |              | 2021-06-07 |             | Second release candidate for B7.8             |
| 1.2.0               | B7.8rc1   |          | 0723         | 2021-05-24 |             | First release candidate for B7.8              |
| 1.1.0               | B7.7.1    | 2021.1.0 | 0682         | 2021-02-26 |             | Final release candidate for B7.7.1            |
| 1.0.0               | B7.7.1rc1 |          | 0678         | 2021-02-22 |             | First release candidate for B7.7.1            |
| 0.18.3              | B7.7      | 2020.4.0 | 0670         | 2021-01-25 |             | Final release candidate for B7.7              |
| 0.18.2              | B7.7rc3   |          | 0668         | 2021-01-19 |             | Third release candidate for B7.7              |
| 0.18.1              | B7.7rc2   |          | 0664         | 2021-01-08 |             | Second release candidate for B7.7             |
| 0.18.0              | B7.7rc1   |          | 0645         | 2020-12-21 |             | First release candidate for B7.7              |
| 0.17.1              | B7.6      | 2020.3.0 | 0641         | 2020-09-15 |             | Final release candidate for B7.6              |
| 0.17.0              | B7.6rc1   |          | 0637         | 2020-08-28 |             | First release candidate for B7.6              |
| 0.16.2              | B7.5      | 2020.2.0 | 0619         | 2020-06-10 |             | Same as 0.16.1, but with installation bug fix |
| 0.16.1              | B7.5      | 2020.2.0 | 0619         | 2020-05-19 |             | Final release candidate for B7.5              |
| 0.16.0              | B7.5rc1   |          | 0614         | 2020-05-04 |             | First release candidate for B7.5              |
| 0.15.1              | B7.4.2    | 2020.1.0 | 0586         | 2020-03-10 |             | Final release candidate for B7.4.2            |
| 0.15.0              | B7.4.2rc1 |          | 0585         | 2020-02-28 |             | First release candidate for B7.4.2            |
| 0.14.2              | B7.4      | 2019.3.0 | 0570         | 2019-11-18 |             | Final release candidate for B7.4              |
| 0.14.1              | B7.4rc2   |          | 0568         | 2019-11-11 |             | Second release candidate for B7.4             |
| 0.14.0              | B7.4rc1   |          | 0563         | 2019-10-25 |             | First release candidate for B7.4              |
| 0.13.8              | B7.3.1    | 2019.2.0 | 0541         | 2019-09-05 |             | Patch for Build 7.3 released as Build 7.3.1   |
| 0.13.7              | B7.3      | 2019.1.0 | 0535         | 2019-06-21 |             | Final release candidate for Build 7.3         |
| 0.13.6              | B7.3rc4   |          | 0534         | 2019-06-20 |             | Fourth release candidate for Build 7.3        |
| 0.13.5              | B7.3rc3   |          | 0534         | 2019-06-19 |             | Third release candidate for Build 7.3         |
| 0.13.4              | B7.3rc2   |          | 0534         | 2019-06-18 |             | Second release candidate for Build 7.3        |
| 0.13.3              | B7.3rc1   |          | 0532         | 2019-06-04 |             | First release candidate for Build 7.3         |
| 0.13.2              |           |          | 0500         | 2019-05-14 |             | DMS test, no delivery to I&T                  |
| 0.13.1              |           |          | 0500         | 2019-03-08 |             | DMS test, no delivery to I&T                  |
| 0.13.0              |           |          | 0500         | 2019-02-15 |             | DMS test, no delivery to I&T                  |
| 0.12.3              | B7.2.1    |          | 0500         | 2019-01-15 |             | DMS Build 7.2.1 patch release                 |
| 0.12.2              | B7.2      | 2018_2   | 0495         | 2018-11-07 |             | Final release candidate for Build 7.2         |
| 0.12.1              | B7.2rc2   |          | 0495         | 2018-11-01 |             | Second release candidate for Build 7.2        |
| 0.12.0              | B7.2rc1   |          | 0493         | 2018-10-09 |             | First release candidate for Build 7.2         |
| 0.11.0              |           |          | 0482         | 2018-09-10 |             | DMS test, no delivery to I&T                  |
| 0.10.0              |           |          | 0477         | 2018-07-31 |             | DMS test, no delivery to I&T                  |
| 0.9.6               | B7.1.3    | 2018_1   | 0468         | 2018-06-08 |             | Final release candidate for Build 7.1.3       |
| 0.9.5               | B7.1.3rc3 |          | 0468         | 2018-06-06 |             | Third release candidate for Build 7.1.3       |
| 0.9.4               | B7.1.3rc2 |          | 0463         | 2018-05-29 |             | Second release candidate for Build 7.1.3      |
| 0.9.3               | B7.1.3rc1 |          | 0457         | 2018-05-11 |             | First release candidate for Build 7.1.3       |
| 0.9.2               |           |          | 0441         | 2018-03-28 |             | DMS test, no delivery to I&T                  |
| 0.9.1               |           |          | 0432         | 2018-02-16 |             | DMS test, no delivery to I&T                  |
| 0.9.0               | B7.1.2    |          | 0422         | 2017-12-22 |             | DMS patch release to I&T 2018-02-15           |
| 0.8.0               | B7.1.1    |          | 0422         | 2017-11-06 |             | DMS patch release to I&T 2018-01-17           |
| 0.8.0               | B7.1      | 2017_1   | 0422         | 2017-11-06 |             | Final release for Build 7.1                   |
| 0.7.7               | B7.0      | 2016_2   | 0303         | 2016-12-13 |             | Final release for Build 7.0                   |
