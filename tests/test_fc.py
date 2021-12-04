from pathlib import Path

import pytest
import numpy as np

from sdhc.fcCalc import fcCalc
from tests.utils import np_load

RESOURCES_PATH = Path("tests").joinpath("resources")


@pytest.mark.lammps
def test_fccalc():

    restart_file = RESOURCES_PATH.joinpath("quenched.restart")

    hstep = 0.001

    fc = fcCalc(fileprefix=None, restartfile=str(restart_file))

    pair_coeff_file = Path("tests").joinpath("resources").joinpath("Si_vbwm.sw")

    ids_L, ids_R = fc.preparelammps(
        pair_style="sw", pair_coeff=f"* * {pair_coeff_file} Si", w_interface=3.0
    )

    # Atoms on the left and right interface
    assert ids_L is not None
    assert ids_R is not None

    # Expected number of atoms on the left interface
    assert len(ids_L) == 60
    ids_L_expected = np_load(file=RESOURCES_PATH.joinpath("force_constants.ids_L.npy"))
    assert np.all(ids_L == ids_L_expected)

    # Expected number of atoms on the right interface
    assert len(ids_R) == 63

    ids_R_expected = np_load(file=RESOURCES_PATH.joinpath("force_constants.ids_R.npy"))
    assert np.all(ids_R == ids_R_expected)

    fc.fcCalc(hstep)

    # Force constant matrix
    Kij = fc.Kij

    assert Kij is not None
    assert Kij.shape == (len(ids_L) * 3, len(ids_R) * 3)

    Kij_expected = np_load(file=RESOURCES_PATH.joinpath("force_constants.Kij.npy"))

    assert np.allclose(Kij, Kij_expected)
