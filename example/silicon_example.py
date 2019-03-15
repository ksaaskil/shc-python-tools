from __future__ import print_function
import numpy as np
from randomAtomBox import atombox
from sdhc import SHCPostProc
from lammps import lammps

QUENCH_STEPS_HEATING = 5e5
QUENCH_STEPS_QUENCH = 1e6
QUENCH_STEPS_COOLED = 5e5

SIMU_STEPS_EQUIL = 5e5
SIMU_STEPS_STEADY = 1e6
SIMU_STEPS_SIMULATION = 1e6

SYSTEM_LENGTH = 200
SYSTEM_WIDTH = 20


def iterateFile(lmp, filename):
    """
    Do the same as lmp.file(filename) but allow the script to be
    continued after quit.
    """
    with open(filename, "r") as f:
        for line in f:
            print(line)
            if "quit" in line and line[0] != "#":
                return
            else:
                lmp.command(line)
    return


def write_initial_positions_file(filename):
    mass = 28.0
    rho = 2.291  # Density in g/cm^3
    n_atoms = np.int(np.round(rho * 1e-3 / 1e-6 * SYSTEM_LENGTH * SYSTEM_WIDTH ** 2 * 1e-30 / (mass * 1.66e-27)))
    ab = atombox(SYSTEM_LENGTH, SYSTEM_WIDTH, n_atoms)
    ab.fillBox(seed=1234)
    ab.writeToFile(filename, mass)


def perform_quench(file_prefix, atom_positions_file, restart_file):
    lmp = lammps()
    lmp.command("variable filename string '" + file_prefix + "'")
    lmp.command("variable datafile string '" + atom_positions_file + "'")
    lmp.command("variable restartfile string '" + restart_file + "'")

    lmp.command("variable steps_heating equal {}".format(QUENCH_STEPS_HEATING))
    lmp.command("variable steps_quench equal {}".format(QUENCH_STEPS_QUENCH))
    lmp.command("variable steps_cooled equal {}".format(QUENCH_STEPS_COOLED))

    iterateFile(lmp, "quench_Si.lmp")
    lmp.close()


def perform_simulation(file_prefix, restart_file):
    lmp = lammps()
    lmp.command("variable filename string '" + file_prefix + "'")
    lmp.command("variable restartfile string '" + restart_file + "'")
    lmp.command("variable steps_equil equal {}".format(SIMU_STEPS_EQUIL))
    lmp.command("variable steps_steady equal {}".format(SIMU_STEPS_STEADY))
    lmp.command("variable steps_simu equal {}".format(SIMU_STEPS_SIMULATION))

    iterateFile(lmp, "amorphous_interface.lmp")
    lmp.close()


def compute_sdhc(file_prefix, restart_file):
    compact_velocities_file = file_prefix + '.vels.dat.compact'
    atomic_velocities_file = file_prefix + '.vels.dat'
    frequency_window_width = 0.5e12

    force_constant_file_prefix = file_prefix
    unit_scaling_factor = 1.602e-19 / (1e-20) * 1e4
    md_timestep = 2.5e-15

    postprocessor = SHCPostProc(compact_velocities_file,
                                force_constant_file_prefix,
                                dt_md=md_timestep,
                                scaleFactor=unit_scaling_factor,
                                LAMMPSDumpFile=atomic_velocities_file,
                                widthWin=frequency_window_width,
                                NChunks=20,
                                chunkSize=50000,
                                backupPrefix=file_prefix,
                                LAMMPSRestartFile=restart_file,
                                reCalcVels=True,
                                reCalcFC=True)

    postprocessor.postProcess()
    return postprocessor


def main(file_prefix):
    """
    Run SDHC example for a-Si.

    :param file_prefix: File prefix, e.g. `090419a`
    :type file_prefix: str
    :return: None
    """

    atom_positions_file = file_prefix + '_Si.dat'
    restart_file = file_prefix + '.quenched.restart'

    write_initial_positions_file(atom_positions_file)

    # Do quenching
    perform_quench(file_prefix, atom_positions_file, restart_file)

    # Gather data from simulation
    perform_simulation(file_prefix, restart_file)

    postprocessor = compute_sdhc(file_prefix, restart_file)

    # Pickling the post-processing object into file
    import cPickle as pickle
    with open(file_prefix + '_PP.pckl', 'w') as f:
        pickle.dump(postprocessor, f)

    # Saving into numpy files 
    np.save(file_prefix + '_oms.npy', postprocessor.oms_fft)
    np.save(file_prefix + '_SHC.npy', postprocessor.SHC_smooth)

    # Saving the frequencies and heat currents to file
    np.savetxt(file_prefix + '_SHC.txt', np.column_stack((postprocessor.oms_fft, postprocessor.SHC_smooth)))


if __name__ == '__main__':
    main(file_prefix='094015a')
