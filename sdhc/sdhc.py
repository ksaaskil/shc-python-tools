"V2 version of heat current calculation"
from dataclasses import dataclass
from pathlib import Path
import typing

import numpy as np

from sdhc.config import logger


def _smoothen(df, func, widthWin):
    Nwindow = int(np.ceil(widthWin / df))
    daniellWindow = np.ones(Nwindow) / Nwindow
    # daniellWindow/=np.sum(daniellWindow)
    # Smooth the value
    smooth = np.convolve(func, daniellWindow, "same")
    return smooth


@dataclass(frozen=True)
class SdhcResult:
    oms_fft: np.ndarray
    SHC_smooth: np.ndarray
    SHC_smooth2: np.ndarray
    SHC_average: np.ndarray
    SHC_error: typing.Optional[np.ndarray]


def calculate_sdhc(
    compact_vels_file: Path,
    Kij: np.ndarray,
    ids_L: np.ndarray,
    ids_R: np.ndarray,
    dt_md: float = 1.0,
    chunkSize: int = 50000,
    NChunks: int = 20,
    scaleFactor: float = 1.0,
    widthWin: float = 1.0,
) -> SdhcResult:
    """
    Calculate the spectral decomposition.
    """

    NL = len(ids_L)
    NR = len(ids_R)

    with open(compact_vels_file, "r") as f:
        s = f.readline().split()
        NAtoms = int(s[1])

        if NAtoms != (NL + NR):
            raise ValueError(
                f"""
Mismatch in the numbers of atoms in the read velocity file and the used force constant file:
velocity file has {NAtoms} and force constants has: {NL}x{NR}
"""
            )

        s = f.readline()
        # logger.info(s)
        s = s.split()
        sampleTimestep = int(s[1]) * dt_md

        s = f.readline()  # Atom ids:
        # logger.info(s)
        # Read the atom ids
        _ = np.fromfile(f, dtype=int, count=NAtoms, sep=" ")

        s = f.readline()  # ------
        # logger.info(s)

        # Total number of degrees of freedom
        NDOF = 3 * (NL + NR)

        oms_fft = np.fft.rfftfreq(chunkSize, d=sampleTimestep) * 2 * np.pi
        Nfreqs = np.size(oms_fft)
        # Initialize the spectral heat current arrays
        SHC_smooth = np.zeros(Nfreqs)
        SHC_smooth2 = np.zeros(Nfreqs)
        SHC_average = np.zeros(Nfreqs)
        SHC_error = None

        exitFlag = False

        for k in np.arange(NChunks):  # Start the iteration over chunks
            #        for k in range(0,2): # Start the iteration over chunks
            logger.info("Chunk %d/%d." % (k + 1, NChunks))
            # Read a chunk of velocitites
            velArray = np.fromfile(
                f, dtype=np.dtype("f8"), count=chunkSize * NDOF, sep=" "
            )
            # Prepare for exit if the read size does not match the chunk size
            if np.size(velArray) == 0:
                logger.info("Finished the file, exiting.")
                NChunks = k - 1
                break
            elif np.size(velArray) != chunkSize * NDOF:
                # Reaching the end of file
                chunkSize = int(np.size(velArray) / NDOF)
                if k > 0:  # Not the first chunk
                    NChunks = k - 1
                    break
                else:
                    exitFlag = True
                    oms_fft = np.fft.rfftfreq(chunkSize, d=sampleTimestep) * 2 * np.pi
                    Nfreqs = np.size(oms_fft)
                    logger.info(
                        "Changing chunk size to "
                        + str(int(np.size(velArray) / NDOF))
                        + "!"
                    )

            # Reshape the array so that each row corresponds to different degree of freedom (e.g. particle 1, direction x etc.)
            velArray = np.reshape(velArray, (NDOF, chunkSize), order="F")

            # FFT with respect to the second axis (NOTE THE USE OF RFFT)
            velFFT = np.fft.rfft(velArray, axis=1)
            velFFT *= sampleTimestep

            velsL = np.zeros((3 * NL, Nfreqs), dtype=np.complex128)
            velsR = np.zeros((3 * NR, Nfreqs), dtype=np.complex128)

            velsL[0::3, :] = velFFT[3 * ids_L, :]
            velsL[1::3, :] = velFFT[3 * ids_L + 1, :]
            velsL[2::3, :] = velFFT[3 * ids_L + 2, :]

            velsR[0::3, :] = velFFT[3 * ids_R, :]
            velsR[1::3, :] = velFFT[3 * ids_R + 1, :]
            velsR[2::3, :] = velFFT[3 * ids_R + 2, :]

            # Spectral heat current for the specific chunk
            SHC = np.zeros(Nfreqs)

            for ki in range(1, Nfreqs):  # Skip the first one with zero frequency
                SHC[ki] = (
                    -2.0
                    * np.imag(np.dot(velsL[:, ki], np.dot(-Kij, np.conj(velsR[:, ki]))))
                    / oms_fft[ki]
                )

            # Normalize correctly
            SHC /= chunkSize * sampleTimestep

            # Change units
            SHC *= scaleFactor

            # daniellWindow=np.ones(np.ceil(self.widthWin*2*np.pi/(self.oms_fft[1]-self.oms_fft[0])))
            # daniellWindow/=np.sum(daniellWindow)

            SHC_orig = SHC.copy()
            # Smooth the value
            # SHC=np.convolve(SHC,daniellWindow,'same')
            df = (oms_fft[1] - oms_fft[0]) / (2 * np.pi)
            SHC = _smoothen(df, SHC, widthWin)

            if (
                not exitFlag
            ):  # If Nfreqs has changed, the running averaging cannot be performed
                SHC_smooth = (k * SHC_smooth + SHC) / (k + 1.0)
                # The square
                SHC_smooth2 = (k * SHC_smooth2 + SHC ** 2) / (k + 1.0)
                # The non-smoothened average
                SHC_average = (k * SHC_average + SHC_orig) / (k + 1.0)
                # if self.backupPrefix is not None:
                #     np.save(backupPrefix + "_backup_oms.npy", self.oms_fft)
                #     np.save(self.backupPrefix + "_backup_SHC.npy", self.SHC_smooth)
            elif (
                exitFlag and k == 0
            ):  # First chunk and new chunk size, needs re-initializing the vectors as Nfreqs may have changed
                SHC_smooth = SHC
                SHC_smooth2 = SHC ** 2
                SHC_average = SHC_orig
                NChunks = 1
                break
            else:  # This should never be reached
                raise Exception(
                    "SHCPostProc should not reach here (exitFlag=True and k>0)."
                )

        # Calculate the error estimate at each frequency from the between-chunk variances
        if NChunks > 1:
            logger.info("Calculating error estimates...")
            samplevar = (NChunks / (NChunks - 1.0)) * (SHC_smooth2 - SHC_smooth ** 2)
            SHC_error = np.sqrt(samplevar) / np.sqrt(NChunks)
        else:
            logger.info(
                "Skipping calculating error estimates as the number of chunks was one"
            )
            SHC_error = None

        logger.info("Finished post-processing.")

    result = SdhcResult(
        oms_fft=oms_fft,
        SHC_smooth=SHC_smooth,
        SHC_smooth2=SHC_smooth2,
        SHC_average=SHC_average,
        SHC_error=SHC_error,
    )
    return result
