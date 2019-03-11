# Spectral heat current Python tools

This repository contains Python tools for calculating the spectral heat current from the data produced by non-equilibrium molecular dynamics simulation with [LAMMPS](http://lammps.sandia.gov) software. The relevant theory was published in:

1) K. Sääskilahti, J. Oksanen, J. Tulkki, and S. Volz, [Phys. Rev. B 90, 134312 (2014)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.134312)
2) K. Sääskilahti, J. Oksanen, S. Volz, and J. Tulkki, [Phys. Rev. B 91, 115426 (2015)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.245411)

These codes are meant to help anyone interested in implementing the spectral heat current calculations for their own applications. If you want to use the codes for research purposes, please cite the above-mentioned publications and let me know.

See detailed API documentation in [readthedocs.io](https://shc-python-tools.readthedocs.io/en/latest/).

## Contents

The actual library for computing spectral heat current distributions is found
in the [sdhc](./sdhc) folder. It contains:
- [SHCPostProc.py](./sdhc/SHCPostProc.py): Python class for performing the post-processing
- [calcFC.py](./sdhc/calcFC.py): Class for calculating the force constants (note that the definition of the "left" and "right" interfaces must be the same in the NEMD simulation and in the calculation of the force constants)
- [compactify_vels.cpp](./sdhc/compactify_vels.cpp): C++ script for formatting the LAMMPS's dump velocity file into a more easily readable column file. The program can be compiled by running `make` in `shc-tools` folder (if `g++` is found, otherwise modify `Makefile` such that appropriate compiler is defined in variable CC).

In addition, the root directory contains the script [calcSHC.py](./calcSHC.py) demonstrating how the post-processing class is used and how the data could be saved to file.

## Usage

Ensure that `sdhc` is found in your `PYTHONPATH` and do

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

See the example below for the arguments required by `SHCPostProc`.

## Example

Folder [example](./example) contains a self-contained example for calculating the spectral heat current flowing across a slab of amorphous Si. The script to be run is called `silicon_example.py`. It performs the following steps:

1. prepare a box of atoms,
1. call LAMMPS to perform the quenching procedure contained in LAMMPS input file `quench_Si.lmp`,
1. call LAMMPS to perform the actual NEMD calculation for a-Si using `amorphous_interface.lmp`, and
1. perform the post-processing using `shc-tools`

### Prerequisites

- Using LAMMPS from Python requires that you have built LAMMPS as a dynamically shared library as instructed in the [LAMMPS manual](http://lammps.sandia.gov/doc/Section_python.html)
- Simulation uses the `sw` pair style, which is included in the `MANYBODY` package.
See [here](https://lammps.sandia.gov/doc/Build_package.html) how to include packages in your
LAMMPS build.
- You need to have [numpy](https://docs.scipy.org/doc/numpy/index.html) installed.

## TODO
- Documentation for calcSHC.py
- Documentation for the example
- Preparation of visualization files
- Documentation for fcCalc.py
