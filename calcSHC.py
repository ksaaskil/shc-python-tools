#!/usr/bin/python
# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2015
from __future__ import division
from subprocess import call

import numpy as np

from sdhc.SHCPostProc import SHCPostProc

fileprefix = "270315a"
outputFolder = "DATA/" + fileprefix + "_tar"

# Create the data folder

command = ["mkdir", "-p", outputFolder]
print(" ".join(command))
call(command)

dt_md = 2.5e-15  # Timestep used in MD, affects the frequency grid
widthWin = 0.5e12  # Width of the Daniell smoothing window in Hz
# The velocity dump file from LAMMPS
fileVels = fileprefix + ".vels.dat"
# The compactly formatted velocity file, produced using a C++ script if not found
fileCompactVels = fileprefix + ".vels.dat.compact"
# Post-processor searches/saves file "KijFilePrefix.Kij.npy"
KijFilePrefix = "270115a"
# Use this restart file if the force constant file cannot be found
restartFile = KijFilePrefix + ".quenched.restart"
# Correct the units, this assumes the unit of eV/(A^2)*(A/ps)^2 for the v_iK_{ij}v_j product (LAMMPS metal units)
scaleFactor = 1.602e-19 / (1e-20) * 1e4

# Prepare the post-processor
pP = SHCPostProc(
    fileCompactVels,
    KijFilePrefix,
    dt_md=dt_md,
    scaleFactor=scaleFactor,
    LAMMPSDumpFile=fileVels,
    widthWin=widthWin,
    LAMMPSRestartFile="270115a.quenched.restart",
    NChunks=200,
    chunkSize=50000,
    backupPrefix=fileprefix,
    reCalcVels=False,
    reCalcFC=False,
)
# Post-process
pP.postProcess()  # All variables will be contained in the object pP

# Saving into numpy files
np.save(outputFolder + "/" + fileprefix + "_oms.npy", pP.oms_fft)
np.save(outputFolder + "/" + fileprefix + "_SHC.npy", pP.SHC_smooth)

# Saving to file
print("Saving to file " + outputFolder + "/" + fileprefix + "_SHC.txt")
np.savetxt(
    outputFolder + "/" + fileprefix + "_SHC.txt",
    np.column_stack((pP.oms_fft, pP.SHC_smooth)),
)

command = ["tar", "-czvf", fileprefix + "_tar.tgz", outputFolder]
print(" ".join(command))
call(command)

# Plotting if available
# import matplotlib.pylab as plt
# plt.plot(pP.oms_fft/(2*np.pi*1.0e12),pP.SHC_smooth)
# plt.xlabel('Frequency (THz)')
# plt.ylabel('Spectral current')
# plt.savefig(fileprefix+'_SHC.eps')
