# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2015

from __future__ import division
import numpy as np
from LammpsDataFileWriter import writer


class atombox:
    """
    Used to create a box of atoms at random locations and to write to file.

    Methods:
      __init__: Documentation below
      fillBox(seed): Create the co-ordinate vectors using the seed (default 1234)
      writeToFile(fileToWrite,atomMass): Write the coordinates to a LAMMPS data file called "fileToWrite" and set atom mass to "atomMass".
    """

    def __init__(self, length, width, Natoms):
        """
        Arguments.
          length: Length of the box in chosen units
          width: Width of the box in chosen units
          Natoms: Number of atoms to be stored into the box
        """
        self.length = length
        self.width = width
        self.xs = None
        self.ys = None
        self.zs = None
        self.Natoms = Natoms

    def __enter__(self):
        return self

    def __exit__(self, t1, t2, t3):
        return False

    def fillBox(self, seed=1234):
        np.random.seed(seed)
        self.xs = np.random.rand(self.Natoms) * self.length
        self.ys = np.random.rand(self.Natoms) * self.width
        self.zs = np.random.rand(self.Natoms) * self.width

    def writeToFile(self, filetowrite, atomMass):
        if self.xs is None:
            print("You must fill the box first!")
            return False
        with writer(
            filetowrite=filetowrite,
            xs=self.xs,
            ys=self.ys,
            zs=self.zs,
            xlo=0,
            xhi=self.length,
            ylo=0,
            yhi=self.width,
            zlo=0,
            zhi=self.width,
            masses=[atomMass],
        ) as wr:
            try:
                wr.writeToFile()
                return True
            except:  # Catch any exceptions
                return False


if __name__ == "__main__":

    filetowrite = "test.dat"

    length = 100.0
    width = 50.0
    # Atom mass
    mass = 28.0

    # Calculate number of atoms assuming density
    rho = 2.291  # In g/cm^3
    Natoms = int(
        np.round(rho * 1e-3 / 1e-6 * length * width ** 2 * 1e-30 / (mass * 1.66e-27))
    )

    print("Number of atoms is %d." % (Natoms))

    with atombox(length=length, width=width, Natoms=Natoms) as ab:
        # Generate the coordinates
        ab.fillBox()
        # Write to file
        ab.writeToFile(filetowrite, mass)
