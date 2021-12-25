import logging
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger("fc")


def write_to_file(
    Kij: np.ndarray, ids_L: np.ndarray, ids_R: np.ndarray, fileprefix: str
):
    """
    Write `self.Kij` to files starting with `self.fileprefix`.

    :return: None
    """
    np.save(fileprefix + ".Kij.npy", Kij)
    np.save(fileprefix + ".ids_L.npy", ids_L)
    np.save(fileprefix + ".ids_R.npy", ids_R)


def calculate_force_constants(
    restartfile: Path,
    hstep: float,
    pair_style: str = None,
    pair_coeff: str = None,
    x_interface=0.5,
    w_interface=3.0,
):
    """
    Version 2 of force constant calculation.

    Uses the `Python library interface <https://lammps.sandia.gov/doc/Python_library.html>`_
    of LAMMPS. You need to have (1) `lammps` in your `PYTHONPATH` and
    (2) `liblammps.so` available for the Python package.

    TODO: Split into smaller parts
    """
    from lammps import lammps

    lmp = lammps()

    print(f"Reading restart file: {restartfile}")
    lmp.command(f"read_restart {restartfile} remap")

    if pair_style is not None:
        lmp.command("pair_style " + pair_style)
    if pair_coeff is not None:
        lmp.command("pair_coeff " + pair_coeff)

    lmp.command("fix NVE all nve")

    boxlo, boxhi, xy, yz, xz, periodicity, box_change = lmp.extract_box()

    xlo = boxlo[0]  # type: float
    xhi = boxhi[0]  # type: float
    logger.info("Box is [%f,%f]." % (xlo, xhi))

    assert xlo is not None
    assert xhi is not None

    # The position of the interface, at the middle by default (0.5)
    x_interface = (xlo + xhi) * x_interface

    xmax = x_interface + w_interface
    xmin = x_interface - w_interface

    lmp.command("region middle block %f %f INF INF INF INF" % (xmin, xmax))
    lmp.command("group interface region middle")

    lmp.command("compute fxs interface property/atom fx")
    lmp.command("compute fys interface property/atom fy")
    lmp.command("compute fzs interface property/atom fz")

    # Coordinates ordered by atom ID
    coords_data = lmp.gather_atoms("x", 1, 3)

    assert coords_data is not None

    # Coordinates in a numpy array
    coords = np.array(coords_data[:], dtype=np.dtype("f8"))
    natoms = lmp.get_natoms()  # extract_global("natoms", 0)

    assert natoms is not None
    assert natoms > 0

    coords = np.reshape(coords, (natoms, 3))

    # X-coordinates in a Numpy array
    xs = coords[:, 0]

    # Atom on the left side?
    mask_left = np.logical_and(xs < x_interface, xs > xmin)
    # Atom on the right side?
    mask_right = np.logical_and(xs > x_interface, xs < xmax)

    # Note that these indices differ from atom IDs by a factor of one
    inds_left = np.where(mask_left)[0]
    inds_right = np.where(mask_right)[0]

    # All atom indices sorted by atom ID, duplicates removed
    inds_interface = np.unique(np.concatenate((inds_left, inds_right)))
    # Where are the atoms of the left atom set
    ids_L = np.in1d(inds_interface, inds_left)
    ids_L = np.where(ids_L)[0]
    # Atoms of the right set
    ids_R = np.in1d(inds_interface, inds_right)
    ids_R = np.where(ids_R)[0]

    assert lmp is not None

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
            Kij[3 * i1 + direction, :] = (fc1[inds_right_1d] - fc2[inds_right_1d]) / (
                2.0 * hstep
            )
            xc[index] += hstep
            lmp.scatter_atoms("x", 1, 3, xc)

    return Kij, ids_L, ids_R
