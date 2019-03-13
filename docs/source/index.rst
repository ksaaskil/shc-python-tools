.. sdhc documentation master file, created by
   sphinx-quickstart on Mon Mar 11 15:27:37 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Spectral decomposition of heat current (SDHC) tools
===================================================

This package contains Python tools for calculating the
spectral decomposition of heat current from the data produced by
non-equilibrium molecular dynamics simulation with
`LAMMPS <http://lammps.sandia.gov>`_ software.
The relevant theory was published in:

1) K. Sääskilahti, J. Oksanen, J. Tulkki, and S. Volz,
`Phys. Rev. B 90, 134312 (2014) <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.134312>`_.

2) K. Sääskilahti, J. Oksanen, S. Volz, and J. Tulkki,
`Phys. Rev. B 91, 115426 (2015) <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.91.115426>`_.

The package is meant to help anyone interested in
implementing the spectral heat current calculations
for their own applications. If you want to use the code
for research purposes, please cite the above-mentioned
publications and let me know.

   *Note that I'm currently rewriting both the documentation and code so this is still work in progress!*

Installation
------------

Clone the `GitHub repository <https://github.com/ksaaskil/shc-python-tools>`_
and install the package from the repository root::

   $ pip install -e .

The command installs the library and its dependencies listed in ``setup.py``.

If you wish to use Python package to compute the force constants,
you need to have

I will later also add the package to PyPI so one can simply do ``pip install sdhc``.

Usage
-----

Example::

   from sdhc import SHCPostProc
   import numpy as np

   postprocessor = SHCPostProc(*args, **kwargs)
   postprocessor.postProcess()

   # Save frequencies and smoothened spectral heat currents as NumPy files
   np.save('angular_frequencies.npy', postprocessor.oms_fft)
   np.save('heat_currents.npy', postprocessor.SHC_smooth)

   # Save the frequencies and smoothened spectral heat currents to text file
   np.savetxt('frequencies_and_currents.txt', np.column_stack((oms, postprocessor.SHC_smooth)))


Example for a-Si
----------------

Folder `example <https://github.com/ksaaskil/shc-python-tools/tree/master/example>`_
contains a self-contained example for calculating the
spectral decomposition of heat current flowing across a
slab of amorphous Si.
The script to be run is called ``silicon_example.py``.
It performs the following steps:

1. prepare a box of atoms,
2. call LAMMPS to perform the quenching procedure contained in LAMMPS input file `quench_Si.lmp`,
3. call LAMMPS to perform the actual NEMD calculation for a-Si using `amorphous_interface.lmp`, and
4. perform the post-processing using `sdhc`

Prerequisites
~~~~~~~~~~~~~

- Using LAMMPS from Python requires that you have built LAMMPS as a dynamically shared library as instructed in the `LAMMPS manual <http://lammps.sandia.gov/doc/Section_python.html>`_
- Simulation uses the ``sw`` pair style, which is included in the ``MANYBODY`` package. See `here <https://lammps.sandia.gov/doc/Build_package.html>`_ how to include packages in your LAMMPS build.

API documentation
-----------------

.. automodule:: sdhc
   :members:
   :undoc-members:
