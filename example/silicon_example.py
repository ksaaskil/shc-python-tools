from __future__ import print_function
import numpy as np
from randomAtomBox import atombox
from sdhc import SHCPostProc
from lammps import lammps


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


QUENCH_STEPS_HEATING = 5e5
QUENCH_STEPS_QUENCH = 1e6
QUENCH_STEPS_COOLED = 5e5

SIMU_STEPS_EQUIL = 5e5
SIMU_STEPS_STEADY = 1e6
SIMU_STEPS_SIMULATION = 1e6

SYSTEM_LENGTH = 200
SYSTEM_WIDTH = 20


def main(filePrefix):
    """
    Run SDHC example for a-Si.

    :param filePrefix: File prefix, e.g. `090419a`
    :type filePrefix: str
    :return: None
    """

    dataFile = filePrefix + '_Si.dat'
    restartFile = filePrefix + '.quenched.restart'

    mass = 28.0
    rho = 2.291  # Density in g/cm^3

    Natoms = np.int(np.round(rho * 1e-3 / 1e-6 * SYSTEM_LENGTH * SYSTEM_WIDTH ** 2 * 1e-30 / (mass * 1.66e-27)))
    # Create the box of silicon atoms, write to datafile
    ab = atombox(SYSTEM_LENGTH, SYSTEM_WIDTH, Natoms)
    ab.fillBox(seed=1234)
    ab.writeToFile(dataFile, mass)
    del ab

    # Minimize the atom positions
    lmp = lammps()
    lmp.command("variable filename string '" + filePrefix + "'")
    lmp.command("variable datafile string '" + dataFile + "'")
    lmp.command("variable restartfile string '" + restartFile + "'")

    lmp.command("variable steps_heating equal {}".format(QUENCH_STEPS_HEATING))
    lmp.command("variable steps_quench equal {}".format(QUENCH_STEPS_QUENCH))
    lmp.command("variable steps_cooled equal {}".format(QUENCH_STEPS_COOLED))

    iterateFile(lmp, "quench_Si.lmp")
    lmp.close()

    lmp = lammps()
    lmp.command("variable filename string '" + filePrefix + "'")
    lmp.command("variable restartfile string '" + restartFile + "'")
    lmp.command("variable steps_equil equal {}".format(SIMU_STEPS_EQUIL))
    lmp.command("variable steps_steady equal {}".format(SIMU_STEPS_STEADY))
    lmp.command("variable steps_simu equal {}".format(SIMU_STEPS_SIMULATION))

    iterateFile(lmp, "amorphous_interface.lmp")

    fileCompactVels = filePrefix + '.vels.dat.compact'
    fileVels = filePrefix + '.vels.dat'
    widthWin = 0.5e12

    KijFilePrefix = filePrefix
    scaleFactor = 1.602e-19 / (1e-20) * 1e4
    dt_md = 2.5e-15

    pP = SHCPostProc(fileCompactVels,
                     KijFilePrefix,
                     dt_md=dt_md,
                     scaleFactor=scaleFactor,
                     LAMMPSDumpFile=fileVels,
                     widthWin=widthWin,
                     NChunks=20,
                     chunkSize=50000,
                     backupPrefix=filePrefix,
                     LAMMPSRestartFile=restartFile,
                     reCalcVels=True,
                     reCalcFC=True)

    pP.postProcess()

    # Pickling the post-processing object into file
    import cPickle as pickle
    with open(filePrefix + '_PP.pckl', 'w') as f:
        pickle.dump(pP, f)

    # Saving into numpy files 
    np.save(filePrefix + '_oms.npy', pP.oms_fft)
    np.save(filePrefix + '_SHC.npy', pP.SHC_smooth)

    # Saving the frequencies and heat currents to file
    np.savetxt(filePrefix + '_SHC.txt', np.column_stack((pP.oms_fft, pP.SHC_smooth)))


if __name__ == '__main__':
    main(filePrefix='094015a')
