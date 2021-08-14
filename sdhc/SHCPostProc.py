from dataclasses import dataclass
import os
import subprocess

import numpy as np

from sdhc.config import logger
from sdhc.fcCalc import fcCalc


__all__ = ["SHCPostProc"]


@dataclass(frozen=True)
class SdhcResult:
    smooth: np.ndarray
    smooth2: np.ndarray
    average: np.ndarray
    error: np.ndarray
    oms_fft: np.ndarray


class SHCPostProc:
    """
    Compute the spectral decomposition of heat current from the data
    produced using LAMMPS MD simulation.

    The velocities are read from the "compact" file produced with the C++-code
    compactify_vels.cpp from a LAMMPS dump file. If the dump file does not exist,
    it is produced by calling the binary ``compactify_vels``,
    which must be found in the environment's ``$PATH``.

    Minimal usage in Python::

        pP = SHCPostProc(compactVelocityFile, KijFilePrefix)
        pP.postProcess() # Calculate the heat current spectrum

    Public attributes::

        - SHC_smooth (numpy float array): The chunk-averaged, smoothened spectral heat current
        - SHC_smooth2 (numpy float array): Square of the chunk-averaged, smoothened spectral heat current, used for estimating the error from the between-chunk variance
        - SHC_average (numpy float array): The chunk-averaged spectral heat current without smoothing
        - SHC_error (numpy float array): The estimated error from the between-chunk variance, None if only one chunk evaluated
        - oms_fft (numpy float array): The angular frequency grid (in the units of Hz if dt_md is given in the units of seconds in the initialization)

    - compactVelocityFile (str): Path to file from where the velocities are read. Produced using the binary compactify_vels if the file does not exist. In this case, you must also supply the keyword argument LAMMPSDumpFile containing the velocities produced using LAMMPS.
    - KijFilePrefix (str): The prefix used in trying to find the force constant matrix file KijFilePrefix.Kij.npy. If the file does not exist, the force constant calculator fcCalc is called using the keyword argument LAMMPSRestartFile (which must be supplied in this case).
    - reCalcVels (str): Force recreation of the compact velocity file if True, defaults to False
    - reCalcFC (str): Re-calculate force constants if True, defaults to False
    - dt_md (float): Timestep used in the NEMD simulation (seconds), used for inferring the sampling timestep and the frequency grid (default 1.0)
    - scaleFactor (float): Multiply the spectral heat current by this factor to convert to correct units (default 1.0)
    - LAMMPSDumpFile (str): Use this velocity dump file produced by LAMMPS for post-processing, needed only if the compact velocity file cannot be found (default None)
    - LAMMPSRestartFile (str): Use this restart file for calculating the force constants if the force constant file cannot be found (default None)
    - widthWin (float): Use this width for the smoothing window (Hz) (default 1.0)
    - chunkSize (int): Used chunk size for reading the velocities, affects the frequency grid (default 50000). Performing FFT is faster if chunkSize is a power of 2.
    - NChunks (int): The number of chunks to be read, this should be set to a sufficiently large value if the whole velocity file should be read (default 20)
    - backupPrefix (str): Prefix for pickling the post-processing object to a file after the read of each chunk (default None)
    - hstep (float): The displacement used in calculating the force constants by the finite-displacement method (default 0.001)
    """

    def __init__(
        self,
        compactVelocityFile: str,
        KijFilePrefix: str,
        reCalcVels=False,
        reCalcFC=False,
        dt_md=1.0,
        scaleFactor=1.0,
        LAMMPSDumpFile: str = None,
        LAMMPSRestartFile: str = None,
        widthWin: float = 1.0,
        chunkSize: int = 50000,
        backupPrefix: str = None,
        hstep: float = 0.001,
    ):
        # Attributes set by positional arguments
        self.compactVelocityFile = compactVelocityFile
        self.KijFilePrefix = KijFilePrefix

        self.reCalcVels = reCalcVels
        self.reCalcFC = reCalcFC

        # Other attributes
        self.Kij = None
        self.ids_L = None
        self.ids_R = None
        self.NL = None
        self.NR = None
        self.SHC_smooth = None
        self.SHC_smooth2 = None
        self.SHC_average = None
        self.SHC_error = None
        self.oms_fft = None

        self.dt_md = dt_md
        self.scaleFactor = scaleFactor
        self.LAMMPSDumpFile = LAMMPSDumpFile
        self.LAMMPSRestartFile = LAMMPSRestartFile
        self.backupPrefix = backupPrefix
        self.widthWin = widthWin
        self.chunkSize = chunkSize
        self.NChunks = 20
        self.hstep = hstep

        if self.reCalcVels or not os.path.isfile(
            self.compactVelocityFile
        ):  # Check if the velocity file exists
            # Check that the LAMMPS Dump file exists
            if self.LAMMPSDumpFile is None or not os.path.isfile(self.LAMMPSDumpFile):
                raise ValueError(
                    f"You must give the LAMMPS velocity dump file as an argument to create the file {self.compactVelocityFile}"
                )
            # print self.compactVelocityFile + " does not exist, creating by reading from file " + self.LAMMPSDumpFile
            # Run the C++ script
            self._compactVels(self.LAMMPSDumpFile, self.compactVelocityFile)

        else:
            logger.info(
                self.compactVelocityFile
                + " exists, using the file for post-processing."
            )

        # Check the force constant file
        if self.reCalcFC or not os.path.isfile(
            self.KijFilePrefix + ".Kij.npy"
        ):  # Check if the force constant file exists
            logger.info(f"Creating force constants with prefix: {self.KijFilePrefix}")
            if self.LAMMPSRestartFile is None:
                raise ValueError(
                    f"You must give the LAMMPSRestartFile as an argument so that the file {self.KijFilePrefix}.Kij.npy can be created"
                )
            self._calcFC(self.KijFilePrefix, self.LAMMPSRestartFile)
        else:  # Load the force constants from file
            logger.info(f"Loading force constants from: {self.KijFilePrefix}")
            self._loadFC(self.KijFilePrefix)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def _calcFC(self, fileprefix, restartfile):

        with fcCalc(fileprefix, restartfile) as fc:
            fc.preparelammps(
                pair_style="sw", pair_coeff="* * Si_vbwm.sw Si", w_interface=3.0
            )
            fc.fcCalc(self.hstep)
            fc.writeToFile()
            self.Kij = fc.Kij  # Reference
            assert self.Kij is not None
            logger.info(
                "Size of the Kij file is (3*%d)x(3*%d)."
                % (np.size(self.Kij, 0) / 3, np.size(self.Kij, 1) / 3)
            )
            self.ids_L = fc.ids_L  # Reference
            self.ids_R = fc.ids_R  # Reference

            assert self.ids_L is not None
            assert self.ids_R is not None

            self.NL = len(self.ids_L)
            logger.info(
                "Number of atoms on the left interface: len(ids_L)=%d" % self.NL
            )
            self.NR = len(self.ids_R)
            logger.info(
                "Number of atoms on the right interface: len(ids_R)=%d" % self.NR
            )

    def _loadFC(self, KijFilePrefix):
        kij_file = f"{KijFilePrefix}.Kij.npy"
        logger.info(f"Loading the force constants from {kij_file}")
        self.Kij = np.load(kij_file)
        logger.info(
            "Size of the Kij file is (3*%d)x(3*%d)."
            % (np.size(self.Kij, 0) / 3, np.size(self.Kij, 1) / 3)
        )
        left_ids_file = f"{KijFilePrefix}.ids_L.npy"
        logger.info("Loading left interfacial atom indices from {left_ids_file}")
        self.ids_L = np.load(left_ids_file)
        right_ids_file = f"{KijFilePrefix}.ids_R.npy"
        logger.info(f"Loading right interfacial atom indices from {right_ids_file}")
        self.ids_R = np.load(right_ids_file)
        self.NL = len(self.ids_L)
        logger.info("Number of atoms on the left: len(ids_L)=%d" % self.NL)
        self.NR = len(self.ids_R)
        logger.info("Number of atoms on the right: len(ids_R)=%d" % self.NR)
        if (np.size(self.Kij, 0) / 3 != self.NL) or (
            np.size(self.Kij, 1) / 3 != self.NR
        ):
            raise ValueError("Sizes in Kij and ids_L/R do not match!")

    def _compactVels(self, fileVels, finalFileVels):
        command = ["compactify_vels", fileVels, finalFileVels]
        logger.info(f"Running command: {' '.join(command)}")
        subprocess.run(command)

    def _smoothen(self, df, func, widthWin):
        Nwindow = int(np.ceil(widthWin / df))
        daniellWindow = np.ones(Nwindow) / Nwindow
        # daniellWindow/=np.sum(daniellWindow)
        # Smooth the value
        smooth = np.convolve(func, daniellWindow, "same")
        return smooth

    def postProcess(self):
        """
        Calculate the spectral decomposition and store to ``self.SHC_smooth``.

        """

        logger.info(f"Reading the compact velocity file: {self.compactVelocityFile}")
        f = open(self.compactVelocityFile, "r")
        s = f.readline().split()
        NAtoms = int(s[1])

        assert self.NL is not None
        assert self.NR is not None

        if NAtoms != (self.NL + self.NR):
            raise ValueError(
                f"Mismatch in the numbers of atoms in the read velocity file and the used force constant file: velocity file has {NAtoms} and force constants has: {self.NL}x{self.NR} "
            )

        s = f.readline()
        # logger.info(s)
        s = s.split()
        self.sampleTimestep = int(s[1]) * self.dt_md

        s = f.readline()  # Atom ids:
        # logger.info(s)
        # Read the atom ids
        _ = np.fromfile(f, dtype=int, count=NAtoms, sep=" ")

        s = f.readline()  # ------
        # logger.info(s)

        # Total number of degrees of freedom
        NDOF = 3 * (self.NL + self.NR)

        self.oms_fft = (
            np.fft.rfftfreq(self.chunkSize, d=self.sampleTimestep) * 2 * np.pi
        )
        Nfreqs = np.size(self.oms_fft)
        # Initialize the spectral heat current arrays
        self.SHC_smooth = np.zeros(Nfreqs)
        self.SHC_smooth2 = np.zeros(Nfreqs)
        self.SHC_average = np.zeros(Nfreqs)

        exitFlag = False

        assert self.ids_L is not None
        assert self.ids_R is not None
        assert self.Kij is not None

        for k in np.arange(self.NChunks):  # Start the iteration over chunks
            #        for k in range(0,2): # Start the iteration over chunks
            logger.info("Chunk %d/%d." % (k + 1, self.NChunks))
            # Read a chunk of velocitites
            velArray = np.fromfile(
                f, dtype=np.dtype("f8"), count=self.chunkSize * NDOF, sep=" "
            )
            # Prepare for exit if the read size does not match the chunk size
            if np.size(velArray) == 0:
                logger.info("Finished the file, exiting.")
                self.NChunks = k - 1
                break
            elif np.size(velArray) != self.chunkSize * NDOF:
                # Reaching the end of file
                self.chunkSize = int(np.size(velArray) / NDOF)
                if k > 0:  # Not the first chunk
                    self.NChunks = k - 1
                    break
                else:
                    exitFlag = True
                    self.oms_fft = (
                        np.fft.rfftfreq(self.chunkSize, d=self.sampleTimestep)
                        * 2
                        * np.pi
                    )
                    Nfreqs = np.size(self.oms_fft)
                    logger.info(
                        "Changing chunk size to "
                        + str(int(np.size(velArray) / NDOF))
                        + "!"
                    )

            # Reshape the array so that each row corresponds to different degree of freedom (e.g. particle 1, direction x etc.)
            velArray = np.reshape(velArray, (NDOF, self.chunkSize), order="F")

            # FFT with respect to the second axis (NOTE THE USE OF RFFT)
            velFFT = np.fft.rfft(velArray, axis=1)
            velFFT *= self.sampleTimestep

            velsL = np.zeros((3 * self.NL, Nfreqs), dtype=np.complex128)
            velsR = np.zeros((3 * self.NR, Nfreqs), dtype=np.complex128)

            velsL[0::3, :] = velFFT[3 * self.ids_L, :]
            velsL[1::3, :] = velFFT[3 * self.ids_L + 1, :]
            velsL[2::3, :] = velFFT[3 * self.ids_L + 2, :]

            velsR[0::3, :] = velFFT[3 * self.ids_R, :]
            velsR[1::3, :] = velFFT[3 * self.ids_R + 1, :]
            velsR[2::3, :] = velFFT[3 * self.ids_R + 2, :]

            # Spectral heat current for the specific chunk
            SHC = np.zeros(Nfreqs)

            for ki in range(1, Nfreqs):  # Skip the first one with zero frequency
                SHC[ki] = (
                    -2.0
                    * np.imag(
                        np.dot(velsL[:, ki], np.dot(-self.Kij, np.conj(velsR[:, ki])))
                    )
                    / self.oms_fft[ki]
                )

            # Normalize correctly
            SHC /= self.chunkSize * self.sampleTimestep

            # Change units
            SHC *= self.scaleFactor

            # daniellWindow=np.ones(np.ceil(self.widthWin*2*np.pi/(self.oms_fft[1]-self.oms_fft[0])))
            # daniellWindow/=np.sum(daniellWindow)

            SHC_orig = SHC.copy()
            # Smooth the value
            # SHC=np.convolve(SHC,daniellWindow,'same')
            df = (self.oms_fft[1] - self.oms_fft[0]) / (2 * np.pi)
            SHC = self._smoothen(df, SHC, self.widthWin)

            if (
                not exitFlag
            ):  # If Nfreqs has changed, the running averaging cannot be performed
                self.SHC_smooth = (k * self.SHC_smooth + SHC) / (k + 1.0)
                # The square
                self.SHC_smooth2 = (k * self.SHC_smooth2 + SHC ** 2) / (k + 1.0)
                # The non-smoothened average
                self.SHC_average = (k * self.SHC_average + SHC_orig) / (k + 1.0)
                if self.backupPrefix is not None:
                    np.save(self.backupPrefix + "_backup_oms.npy", self.oms_fft)
                    np.save(self.backupPrefix + "_backup_SHC.npy", self.SHC_smooth)
            elif (
                exitFlag and k == 0
            ):  # First chunk and new chunk size, needs re-initializing the vectors as Nfreqs may have changed
                self.SHC_smooth = SHC
                self.SHC_smooth2 = SHC ** 2
                self.SHC_average = SHC_orig
                self.NChunks = 1
                break
            else:  # This should never be reached
                assert (
                    False
                ), "SHCPostProc should not reach here (exitFlag=True and k>0)."

        # Calculate the error estimate at each frequency from the between-chunk variances
        if self.NChunks > 1:
            logger.info("Calculating error estimates.")
            samplevar = (self.NChunks / (self.NChunks - 1.0)) * (
                self.SHC_smooth2 - self.SHC_smooth ** 2
            )
            self.SHC_error = np.sqrt(samplevar) / np.sqrt(self.NChunks)
        else:
            self.SHC_error = None

        f.close()

        logger.info("Finished post-processing.")
