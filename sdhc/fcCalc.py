# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2015
import numpy as np

__all__ = ["fcCalc"]


class fcCalc:
    """
    Legacy class for computing force constants between atoms.
    Uses the `Python library interface <https://lammps.sandia.gov/doc/Python_library.html>`_
    of LAMMPS so you need to have (1) `lammps` in your `PYTHONPATH` and
    (2) `liblammps.so` available for the Python package.

    :param fileprefix: File prefix (TODO What is this)
    :type fileprefix: str
    :param restartfile: LAMMPS restart file (TODO What is this)
    :type restartfile: str
    """

    def __init__(self, fileprefix, restartfile):
        self.fileprefix = fileprefix
        self.restartfile = restartfile
        self.Kij = None
        self.inds_L = None
        self.inds_R = None
        self.ids_L = None
        self.ids_R = None
        self.Natoms = None
        self.lmp = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def preparelammps(
        self, pair_style=None, pair_coeff=None, x_interface=0.5, w_interface=3.0
    ):
        """
        Prepare the LAMMPS object for computing force constants.

        :param pair_style: LAMMPS pair style to set
        :type pair_style: str, optional
        :param pair_coeff: LAMMPS `pair_coeff` to set
        :type pair_coeff: str, optional
        :param x_interface: Position of the interface relative to `boxxlo` and `boxxhi`, defaults to 0.5
        :type x_interface: float, optional
        :param w_interface: Width of the area of atoms to include in the interface, defaults to 3.0
        :type w_interface: float, optional
        :return (ids_L, ids_R): tuple of atom indices on left and right
        """
        from lammps import lammps

        restartfile = self.restartfile
        self.lmp = lammps()

        print(f"Reading restart file: {restartfile}")
        self.lmp.command(f"read_restart {restartfile} remap")

        if pair_style is not None:
            self.lmp.command("pair_style " + pair_style)
        if pair_coeff is not None:
            self.lmp.command("pair_coeff " + pair_coeff)

        self.lmp.command("fix NVE all nve")

        assert self.lmp is not None
        # xlo = self.lmp.extract_box.extract_global("boxxlo", 1)
        # xhi = self.lmp.extract_global("boxxhi", 1)
        boxlo, boxhi, xy, yz, xz, periodicity, box_change = self.lmp.extract_box()

        xlo = boxlo[0]  # type: float
        xhi = boxhi[0]  # type: float
        print("Box is [%f,%f]." % (xlo, xhi))

        assert xlo is not None
        assert xhi is not None
        # The position of the interface, at the middle by default (0.5)
        x_interface = (xlo + xhi) * x_interface

        xmax = x_interface + w_interface
        xmin = x_interface - w_interface

        self.lmp.command("region middle block %f %f INF INF INF INF" % (xmin, xmax))
        self.lmp.command("group interface region middle")

        self.lmp.command("compute fxs interface property/atom fx")
        self.lmp.command("compute fys interface property/atom fy")
        self.lmp.command("compute fzs interface property/atom fz")

        # Coordinates ordered by atom ID
        coords_data = self.lmp.gather_atoms("x", 1, 3)

        assert coords_data is not None

        # Coordinates in a numpy array
        coords = np.array(coords_data[:], dtype=np.dtype("f8"))
        self.natoms = self.lmp.get_natoms()  # extract_global("natoms", 0)

        assert self.natoms is not None
        assert self.natoms > 0

        coords = np.reshape(coords, (self.natoms, 3))

        # X-coordinates in a Numpy array
        xs = coords[:, 0]

        # Atom on the left side?
        mask_left = np.logical_and(xs < x_interface, xs > xmin)
        # Atom on the right side?
        mask_right = np.logical_and(xs > x_interface, xs < xmax)

        # Note that these indices differ from atom IDs by a factor of one
        self.inds_left = np.where(mask_left)[0]
        self.inds_right = np.where(mask_right)[0]

        # All atom indices sorted by atom ID, duplicates removed
        inds_interface = np.unique(np.concatenate((self.inds_left, self.inds_right)))
        # Where are the atoms of the left atom set
        self.ids_L = np.in1d(inds_interface, self.inds_left)
        self.ids_L = np.where(self.ids_L)[0]
        # Atoms of the right set
        self.ids_R = np.in1d(inds_interface, self.inds_right)
        self.ids_R = np.where(self.ids_R)[0]

        return self.ids_L, self.ids_R

    def fcCalc(self, hstep) -> np.ndarray:
        """
        Compute force constants and store to `self.Kij`.

        :param hstep: Step to use in finite differences
        :type hstep: float
        :return: Kij: Force constant matrix
        """
        lmp = self.lmp

        assert lmp is not None
        inds_left = self.inds_left
        inds_right = self.inds_right
        # One-dimensional indices of the atoms on the right side
        inds_right_1d = np.concatenate(
            (3 * inds_right, 3 * inds_right + 1, 3 * inds_right + 2)
        )
        inds_right_1d = np.sort(inds_right_1d)

        Kij = np.zeros((len(inds_left) * 3, len(inds_right) * 3))

        # Loop over the atoms on the left side
        for i1 in range(0, len(inds_left)):
            #        for i1 in range(0,10):
            # Index of the atom on the left
            ind1 = inds_left[i1]
            print("\n Moving atom %i/%i. \n" % (i1 + 1, len(inds_left)))

            # Move atom to directions x, y, and z
            for direction in [0, 1, 2]:
                # Index of the displaced degree of freedom
                index = 3 * ind1 + direction
                # Get the coordinates from LAMMPS
                # ids = lmp.gather_atoms("id", 0, 1)

                xc = lmp.gather_atoms("x", 1, 3)

                assert xc is not None

                # Move the atom
                xc[index] += hstep

                # Communicate to LAMMPS
                # print("Scattering to", xc)
                # This will fail with a warning if `atom_modify map` has not been set!
                lmp.scatter_atoms("x", 1, 3, xc)

                # Run LAMMPS to update the forces
                lmp.command("run 0 post no")
                # Gather the forces
                fc1 = lmp.gather_atoms("f", 1, 3)
                assert fc1 is not None

                fc1 = np.array(fc1, dtype=np.dtype("f8"))

                # Move to negative direction
                xc[index] -= 2 * hstep
                lmp.scatter_atoms("x", 1, 3, xc)
                lmp.command("run 0 post no")

                fc2 = lmp.gather_atoms("f", 1, 3)
                fc2 = np.array(fc2, dtype=np.dtype("f8"))

                # Fill one row of spring constant matrix
                Kij[3 * i1 + direction, :] = (
                    fc1[inds_right_1d] - fc2[inds_right_1d]
                ) / (2.0 * hstep)
                xc[index] += hstep
                lmp.scatter_atoms("x", 1, 3, xc)

        self.Kij = Kij
        return Kij

    def writeToFile(self):
        """
        Write `self.Kij` to files starting with `self.fileprefix`.

        :return: None
        """
        np.save(self.fileprefix + ".Kij.npy", self.Kij)
        np.save(self.fileprefix + ".ids_L.npy", self.ids_L)
        np.save(self.fileprefix + ".ids_R.npy", self.ids_R)


if __name__ == "__main__":
    # import argparse
    # parser=argparse.ArgumentParser()
    # parser.add_argument("filePrefix",help="The prefix of file for which to calculate the force constants")
    # parser.add_argument("hstep",default=0.001,help="The displacement used in the finite-difference evaluation of force constants.")

    # args=parser.parse_args()
    # fileprefix=args.filePrefix
    # hstep=args.hstep
    fileprefix = "270115a"
    restartfile = fileprefix + ".quenched.restart"
    hstep = 0.001

    with fcCalc(fileprefix, restartfile) as fc:
        fc.preparelammps(
            pair_style="sw", pair_coeff="* * Si_vbwm.sw Si", w_interface=3.0
        )
        fc.fcCalc(hstep)
        fc.writeToFile()
