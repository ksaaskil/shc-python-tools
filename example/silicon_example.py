from pathlib import Path
import os

import numpy as np

# from randomAtomBox import atombox
from sdhc import SHCPostProc
from sdhc.utils import make_atombox
from lammps import lammps

SCALE = 10

QUENCH_STEPS_HEATING = 5e5 / SCALE
QUENCH_STEPS_QUENCH = 1e6 / SCALE
QUENCH_STEPS_COOLED = 5e5 / SCALE

SIMU_STEPS_EQUIL = 5e5 / SCALE
SIMU_STEPS_STEADY = 1e6 / SCALE
SIMU_STEPS_SIMULATION = 1e6 / SCALE

SYSTEM_LENGTH = 200
SYSTEM_WIDTH = 20

QUENCH_LMP_PATH = Path(__file__).parent.joinpath("quench.lmp")
SIMULATION_LMP_PATH = Path(__file__).parent.joinpath("simulation.lmp")


def iterateFile(lmp: lammps, filename: Path):
    """
    Do the same as lmp.file(filename) but allow the script to be
    continued after quit.
    """
    with filename.open("r") as f:
        for line in f:
            print(line)
            if "quit" in line and line[0] != "#":
                return
            else:
                lmp.command(line)
    return


def write_initial_positions_file(filename: Path):
    mass = 28.0
    rho = 2.291  # Density in g/cm^3
    n_atoms = int(
        np.round(
            rho
            * 1e-3
            / 1e-6
            * SYSTEM_LENGTH
            * SYSTEM_WIDTH ** 2
            * 1e-30
            / (mass * 1.66e-27)
        )
    )
    make_atombox(
        length=SYSTEM_LENGTH,
        width=SYSTEM_WIDTH,
        n_atoms=n_atoms,
        atom_mass=mass,
        output=filename,
    )


def perform_quench(folder, atom_positions_file, restart_file):
    lmp = lammps()
    file_prefix = os.path.join(folder, "quench")
    lmp.command("variable filename string '" + file_prefix + "'")
    lmp.command("variable datafile string '" + atom_positions_file + "'")
    lmp.command("variable restartfile string '" + restart_file + "'")

    lmp.command("variable steps_heating equal {}".format(QUENCH_STEPS_HEATING))
    lmp.command("variable steps_quench equal {}".format(QUENCH_STEPS_QUENCH))
    lmp.command("variable steps_cooled equal {}".format(QUENCH_STEPS_COOLED))

    iterateFile(lmp, QUENCH_LMP_PATH)
    lmp.close()


def perform_simulation(folder, restart_file):
    file_prefix = os.path.join(folder, "simu")
    lmp = lammps()
    lmp.command("variable filename string '" + file_prefix + "'")
    lmp.command("variable restartfile string '" + restart_file + "'")
    lmp.command("variable steps_equil equal {}".format(SIMU_STEPS_EQUIL))
    lmp.command("variable steps_steady equal {}".format(SIMU_STEPS_STEADY))
    lmp.command("variable steps_simu equal {}".format(SIMU_STEPS_SIMULATION))

    iterateFile(lmp, SIMULATION_LMP_PATH)
    lmp.close()


def compute_sdhc(folder, restart_file):

    compact_velocities_file = os.path.join(folder, "vels.dat.compact")
    atomic_velocities_file = os.path.join(folder, "simu.vels.dat")
    frequency_window_width = 0.5e12

    backup_prefix = os.path.join(folder, "backup")
    force_constant_file_prefix = os.path.join(folder, "force_constants")
    unit_scaling_factor = 1.602e-19 / 1e-20 * 1e4
    md_timestep = 2.5e-15

    postprocessor = SHCPostProc(
        compact_velocities_file,
        force_constant_file_prefix,
        dt_md=md_timestep,
        scaleFactor=unit_scaling_factor,
        LAMMPSDumpFile=atomic_velocities_file,
        widthWin=frequency_window_width,
        NChunks=20,
        chunkSize=50000,
        backupPrefix=backup_prefix,
        LAMMPSRestartFile=restart_file,
        reCalcVels=True,
        reCalcFC=True,
    )

    postprocessor.postProcess()
    return postprocessor


def main(folder: Path = Path("lammps-output")):
    """
    Run SDHC example for a-Si.

    folder: Path to folder where to store output. For example, `090419a`
    """

    folder.mkdir(exist_ok=True)

    atom_positions_file = folder.joinpath("Si.dat")
    write_initial_positions_file(atom_positions_file)

    # Do quenching
    restart_file = str(folder.joinpath("quenched.restart"))
    perform_quench(str(folder), atom_positions_file, restart_file)

    # Gather data from simulation
    perform_simulation(str(folder), restart_file)

    postprocessor = compute_sdhc(str(folder), restart_file)

    # Saving into numpy files
    np.save(folder.joinpath("oms.npy"), postprocessor.oms_fft)
    np.save(folder.joinpath("SHC.npy"), postprocessor.SHC_smooth)

    # Saving the frequencies and heat currents to file
    np.savetxt(
        folder.joinpath("SHC.txt"),
        np.column_stack((postprocessor.oms_fft, postprocessor.SHC_smooth)),
    )


if __name__ == "__main__":
    main()
