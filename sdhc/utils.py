# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2021
from datetime import date
from pathlib import Path

import numpy as np


def write_lammps_data_file(
    output: Path,
    xs,
    ys,
    zs,
    xlo=None,
    ylo=None,
    zlo=None,
    xhi=None,
    yhi=None,
    zhi=None,
    tags=None,
    masses=None,
):
    """
    Write LAMMPS data file.

    tags (list) : List of atom types (defaults to a list of ones)
    masses (list): List of atom masses, length must coincide with the maximum tag (defaults to one for each type)
    header (str): The header written in the first line of the output file, prints the date on default (# included automatically in the beginning.).
    """
    today = date.today()
    header = (
        "File written on "
        + str(today.year)
        + "-"
        + str(today.month)
        + "-"
        + str(today.day)
    )

    Natoms = len(xs)
    if Natoms != len(ys) or Natoms != len(zs):
        raise AttributeError("xs, ys and zs must have the same size.")

    tags = tags or [1] * Natoms

    if len(tags) != Natoms:
        raise ValueError(
            "The length of tag list must coincide with the length of xs,ys,zs."
        )

    maxtag = max(tags)

    masses = masses or [1] * maxtag

    if len(masses) != maxtag:
        raise AttributeError(
            "The number of given masses must coincide with the maximum tag."
        )

    shift = 0.1
    xlo = xlo or min(xs) - shift
    xhi = xhi or max(xs) + shift
    ylo = ylo or min(ys) - shift
    yhi = yhi or max(ys) + shift
    zlo = zlo or min(zs) - shift
    zhi = zhi or max(zs) + shift

    print("Writing to file %s." % output)
    with output.open("w") as f:

        f.write("# " + header + "\n")
        f.write("%d atoms\n" % Natoms)
        f.write("%d atom types\n" % max(tags))

        f.write("\n")

        f.write("%.2f %.2f xlo xhi\n" % (xlo, xhi))
        f.write("%.2f %.2f ylo yhi\n" % (ylo, yhi))
        f.write("%.2f %.2f zlo zhi\n" % (zlo, zhi))

        f.write("\nMasses\n\n")
        for i in range(0, max(tags)):
            f.write("%d %.3f\n" % (i + 1, masses[i]))

        f.write("\nAtoms\n\n")

        for i in range(0, Natoms):
            f.write("%d %d %.5f %.5f %.5f\n" % (i + 1, tags[i], xs[i], ys[i], zs[i]))

        print(f"Finished writing to file {output}")


def make_atombox(
    length: float, width: float, n_atoms: int, atom_mass: float, output: Path
):
    xs = np.random.rand(n_atoms) * length
    ys = np.random.rand(n_atoms) * width
    zs = np.random.rand(n_atoms) * width

    write_lammps_data_file(
        output=output,
        xs=xs,
        ys=ys,
        zs=zs,
        xlo=0,
        xhi=length,
        ylo=0,
        yhi=width,
        zlo=0,
        zhi=width,
        masses=[atom_mass],
    )
