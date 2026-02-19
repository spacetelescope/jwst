# JWST Calibration Pipeline

[![Build Status](https://github.com/spacetelescope/jwst/actions/workflows/ci.yml/badge.svg)](https://github.com/spacetelescope/jwst/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/spacetelescope/jwst/branch/main/graph/badge.svg?token=Utf5Zs9g7z)](https://codecov.io/gh/spacetelescope/jwst)
[![Documentation Status](https://readthedocs.org/projects/jwst-pipeline/badge/?version=latest)](https://jwst-pipeline.readthedocs.io/en/latest/?badge=latest)
[![Pre-Commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Powered by STScI Badge](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](https://www.stsci.edu)
[![Powered by Astropy Badge](https://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](https://www.astropy.org/)
[![DOI](https://zenodo.org/badge/60551519.svg)](https://zenodo.org/badge/latestdoi/60551519)
[![Python Versions](https://img.shields.io/pypi/pyversions/jwst)](https://pypi.org/project/jwst/)

![STScI Logo](docs/_static/stsci_logo.png)

This package processes uncalibrated data for both imagers and spectrographs onboard the James Webb Space Telescope (JWST), an orbiting infrared observatory stationed at Earth-Sun L<sub>2</sub>. It performs a series of calibration steps that result in standard data products usable for science.
More information on running this pipeline, including explanations of specific stages and how to obtain reference files,
can be found [here](https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline).

> [!IMPORTANT]
> The JWST calibration pipeline currently supports Linux and macOS.
> Native Windows builds are **not** currently supported; [use WSL instead](https://stenv.readthedocs.io/en/latest/windows.html).

## Installation

We recommend using an isolated Python environment to install `jwst`.

> [!TIP]
> Python "environments" are isolated Python installations, confined to a single directory, where you can install packages, dependencies, and tools without cluttering your system Python libraries. Examples of commonly-used environments are `virtualenv`, `conda`, `mamba`, `uv`, etc.

You can create multiple environments to allow for switching between different versions of the `jwst` package (e.g. a released version versus the current development version).

In all cases, installing to an environment is a 3-step process:
* Create a new environment
* Activate that environment (change shell variables in the current session to point to the isolated Python installation)
* Install the desired version of the `jwst` package into that environment

```shell
virtualenv jwst_env/
source jwst_env/bin/activate
pip install jwst==1.20.1
```

> [!WARNING]
> Installing `jwst` version `1.15.1` through `1.16.1` pulls an incompatible version of `gwcs`.
> Remedy this issue by downgrading `gwcs`:
> ```shell
> pip install 'gwcs<0.22'
> ```

> [!NOTE]
> Please contact the [JWST Help Desk](https://jwsthelp.stsci.edu) if you have issues installing `jwst`.

### Installing the latest unreleased development version directly from the source code

> [!IMPORTANT]
> You need a C compiler in order to build the JWST calibration pipeline (and dependencies) from source.

```shell
pip install git+https://github.com/spacetelescope/jwst
```

### Installing a DMS Operational Build

There may be occasions where an exact copy of an operational DMS build is
desired (e.g. for validation testing or debugging issues in operations).
We package releases for DMS builds [as environment snapshots that specify the exact versions of all packages to be installed](https://ssb.stsci.edu/stasis/releases/jwst/).

Note that these environment files are delivered as Conda YAMLs. To install one, you will need to have `conda`, `mamba`, or `micromamba` installed on your computer.

The [Software vs DMS build version map](#software-vs-dms-build-version-map), shown below, displays the `jwst` version for each DMS build.
For example, use `jwst==1.17.1` for DMS build **11.2**.
Also note that Linux and macOS systems require different snapshot files:
```shell
conda env create --file https://ssb.stsci.edu/stasis/releases/jwst/JWSTDP-1.18.1/delivery/latest-py312-macos-arm64.yml
conda activate JWSTDP-1.18.1-1-py312-macos-arm64
```

> [!NOTE]
> Starting with `jwst==1.16.1`, the JWST pipeline uses [`stasis`](https://github.com/spacetelescope/stasis) to package environments and deliver releases. If you need a version of `jwst` prior to `1.16.1`, use a slightly different procedure:
> ```shell
> conda create -n jwstdp-1.16.0 --file https://ssb.stsci.edu/releases/jwstdp/1.16.0/conda_python_macos-stable-deps.txt
> conda activate jwstdp-1.16.0
> pip install -r https://ssb.stsci.edu/releases/jwstdp/1.16.0/reqs_macos-stable-deps.txt
> ```

## Calibration References Data System (CRDS) Setup

The [Calibration References Data System (CRDS)](https://jwst-crds.stsci.edu/static/users_guide/index.html) serves and manages reference files for several telescope calibration pipelines, including `jwst`.

The JWST CRDS server is available at https://jwst-crds.stsci.edu

> [!WARNING]
> The CRDS PUB Server (`https://jwst-crds-pub.stsci.edu`) was decommissioned in March 2023.
> To use historical files from the PUB server, contact the [JWST Help Desk](https://jwsthelp.stsci.edu).

Set `CRDS_SERVER_URL` and `CRDS_PATH` to run the calibration pipeline with access to reference files from CRDS:
```shell
export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
export CRDS_PATH=$HOME/data/crds_cache/
```

The pipeline will automatically download individual reference files and cache them in the `CRDS_PATH` directory.
Expect to use upwards of 200 gigabytes of disk space for reference files.

> [!TIP]
> If you are inside the STScI network, physically or via VPN, you do not need to set the `CRDS_PATH` environment variable (it defaults to shared network storage).

To use a specific CRDS context other than that [automatically associated with a given pipeline version](https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline/crds-migration-to-quarterly-calibration-updates), you may also explicitly set the ``CRDS_CONTEXT`` environment variable:
```shell
export CRDS_CONTEXT=jwst_1179.pmap
```

## Documentation

Package documentation is available at https://jwst-pipeline.readthedocs.io

For more user-focused documentation, see the JDox pages at https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline

> [!TIP]
> Information on the latest build is available on JDox at https://jwst-docs.stsci.edu/jwst-science-calibration-pipeline/jwst-operations-pipeline-build-information

> [!TIP]
> Public API and its deprecation policy is described at https://jwst-pipeline.readthedocs.io/en/stable/jwst/user_documentation/more_information.html

## Contributions and Feedback

We welcome contributions and feedback on the project. Please follow [`CONTRIBUTING.md`](CONTRIBUTING.md) to submit an issue or a pull request.

We strive to provide a welcoming community to all of our users by abiding with [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md)

If you have questions or concerns regarding the software, please open an issue at https://github.com/spacetelescope/jwst/issues or contact the [JWST Help Desk](https://jwsthelp.stsci.edu).


## Software vs DMS build version map

The table below provides information on each release of the `jwst` package and its relationship to software builds used in the STScI JWST DMS operations environment.
The `Released` column gives the date on which the `jwst` tag was released on PyPI, and the `Ops Install` column gives the date on which the build incorporating that release was installed in DMS operations.
Note that the `CRDS_CONTEXT` listed is a minimum context that can be used with that release. A release should work with any contexts between the specified context and the context specified for the next release.

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


## Testing

See [`TESTING.md`](./TESTING.md)
