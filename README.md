# Spectral decomposition of heat current in Python

This repository contains Python code for calculating the spectral decomposition of heat current from the data produced by non-equilibrium molecular dynamics simulation with [LAMMPS](http://lammps.sandia.gov) software. The relevant equations were published in:

1) K. Sääskilahti, J. Oksanen, J. Tulkki, and S. Volz, [Phys. Rev. B 90, 134312 (2014)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.134312)
2) K. Sääskilahti, J. Oksanen, S. Volz, and J. Tulkki, [Phys. Rev. B 91, 115426 (2015)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.245411)

These codes are meant to help anyone interested in implementing the spectral heat current decomposition calculations for their own applications. If you want to use the codes for research purposes, please cite the above-mentioned publications and let me know.

## Caveats

- This is research code and not unit-tested
- The logic how atom IDs are tracked throughtout the process is fragile. For example, the interface atoms are computed independently both when starting the LAMMPS simulation (to know which velocities to dump) and when calculating force constants afterwards. Currently these atom IDs must match.

## Installation

```bash
$ git clone git@github.com:ksaaskil/shc-python-tools.git
$ cd shc-python-tools
# If you need development dependencies:
$ pip install -e '.[dev]'
# If not:
$ pip install -e .
```

## Prerequisites

- Using LAMMPS from Python requires that you have built LAMMPS as a dynamically shared library as instructed in the [LAMMPS manual](http://lammps.sandia.gov/doc/Section_python.html)
- You need to build `compactify_vels.cpp` in `scripts` folder and have that available in your `$PATH`.
- Simulation uses the `sw` pair style, which is included in the `MANYBODY` package.
See [here](https://lammps.sandia.gov/doc/Build_package.html) how to include packages in your
LAMMPS build.

## Contents

The actual library for computing spectral heat current distributions is found
in the [`sdhc`](./sdhc) folder. It contains:

- [`SHCPostProc.py`](./sdhc/SHCPostProc.py): Python class for performing the post-processing
- [`calcFC.py`](./sdhc/calcFC.py): Class for calculating the force constants (note that the definition of the "left" and "right" interfaces must be the same in the NEMD simulation and in the calculation of the force constants)

Additionally, the `scripts` folder contains:

- [`compactify_vels.cpp`](./scripts/compactify_vels.cpp): C++ script for formatting the LAMMPS's dump velocity file into a more easily readable column file. The program can be compiled by running `make` in `scripts` folder (if `g++` is found, otherwise modify `Makefile` such that appropriate compiler is defined in variable `CC`).

In addition, the root directory contains the script [calcSHC.py](./calcSHC.py) demonstrating how the post-processing class is used and how the data could be saved to file.

## Usage

```python
from sdhc import SHCPostProc
import numpy as np

postprocessor = SHCPostProc(*args, **kwargs)
postprocessor.postProcess()

# Save frequencies and smoothened spectral heat currents as NumPy files
np.save('angular_frequencies.npy', postprocessor.oms_fft)
np.save('heat_currents.npy', postprocessor.SHC_smooth)

# Save the frequencies and smoothened spectral heat currents to text file
np.savetxt('frequencies_and_currents.txt', np.column_stack((oms, postprocessor.SHC_smooth)))
```

## Example

Folder [example](./example) contains a self-contained example for calculating the spectral heat current flowing across a slab of amorphous Si. The script to be run is called `silicon_example.py`. It performs the following steps:

1. prepares a box of atoms,
1. call LAMMPS to perform the quenching procedure contained in LAMMPS input file `run_quench.lmp`,
1. call LAMMPS to perform the actual NEMD calculation for a-Si using `run_simulation.lmp`, and
1. perform the post-processing using `sdhc`

## Development

Format files using Black:

```bash
$ black sdhc example
```

Check code style with Flake8:

```bash
$ flake8 .
```

## Installing LAMMPS on macOS (2021)

Download the tarball from [LAMMPS downloads](https://www.lammps.org/download.html)

```bash
$ mkdir ~/lammps
$ cd ~/lammps
$ wget https://download.lammps.org/tars/lammps.tar.gz
$ tar -xvf lammps.tar.gz
$ cd lammps-29Oct20
```

Edit the target `Makefile` such as `src/MAKE/Makefile.serial` to use `clang++`:

```Makefile
# Makefile.serial
CC =            clang++ -std=c++11 -stdlib=libc++
LINK =          clang++
```

Build as shared library and [include any required packages](https://docs.lammps.org/Build_package.html):

```bash
$ cd src
$ make yes-MANYBODY
$ make mode=shared serial
```

Test the executable:

```bash
$ ./lmp_serial -i ../examples/min/in.min
```

Make the executable available in your `PATH`:

```bash
# Assuming you have ~/bin in your PATH
$ ln -sf ${LAMMPS_PATH}/src/lmp_serial ~/bin/lmp_serial
```

Add LAMMPS to `LD_LIBRARY_PATH`, `DYLD_LIBRARY_PATH` and `PYTHONPATH`:

```bash
# .bash_profile
export LAMMPS_PATH=${HOME}/lammps/lammps-30Jul2021
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LAMMPS_PATH}/src
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${LAMMPS_PATH}/src
export PYTHONPATH=${PYTHONPATH}:${LAMMPS_PATH}/python
export PATH=$PATH:${HOME}/git/shc-python-tools/scripts
```

For some reason, setting library paths didn't work as `lammps` Python package was searching for the `liblammps.so` file in the folder of the Python package. So I added a soft link:

```bash
$ ln -sf ${LAMMPS_PATH}/src/liblammps.so ${LAMMPS_PATH}/python/lammps/liblammps.so
```

I also needed to update `python/core.py` as follows:

```python
  if platform.system() == "Darwin":
     # lib_ext = ".dylib"
     lib_ext = ".so"
```

Now test running LAMMPS from Python:

```python
>>> from lammps import lammps
>>> lmp = lammps()
LAMMPS (30 Jul 2021)
```
