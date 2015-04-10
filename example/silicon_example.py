# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2015

import numpy as np
from randomAtomBox import atombox
from lammps import lammps
from SHCPostProc import SHCPostProc 

filePrefix='090415a'
dataFile=filePrefix+'_Si.dat'
restartFile=filePrefix+'.quenched.restart'
length=200
width=20
mass=28.0
rho=2.291 # Density in g/cm^3

Natoms=np.int(np.round(rho*1e-3/1e-6*length*width**2*1e-30/(mass*1.66e-27)))
# Create the box of silicon atoms, write to datafile
ab=atombox(length,width,Natoms)
ab.fillBox(seed=1234)
ab.writeToFile(dataFile,mass)
del ab

# Minimize the atom positions
lmp=lammps()
lmp.command("variable filename string '"+filePrefix+"'")
lmp.command("variable datafile string '"+dataFile+"'")
lmp.command("variable restartfile string '"+restartFile+"'")

def iterateFile(lmp,filename):
    '''
    Do the same as lmp(filename) but allow the script to be continued after quit.
    '''
    with open(filename,"r") as f:
        for line in f:
            print line
            if "quit" in line and line[0]!="#":
                return True
            else:
                lmp.command(line)
        return True

# lmp.file("quench_Si.lmp") # If quit is found, python quits
iterateFile(lmp,"quench_Si.lmp")
lmp.close()
# Create a new LAMMPS object
lmp=lammps()
lmp.command("variable filename string '"+filePrefix+"'")
lmp.command("variable restartfile string '"+restartFile+"'")
iterateFile(lmp,"amorphous_interface.lmp")
# lmp.file("amorphous_interface.lmp")

fileCompactVels=filePrefix+'.vels.dat.compact'
fileVels=filePrefix+'.vels.dat'
widthWin=0.5e12

KijFilePrefix=filePrefix
scaleFactor=1.602e-19/(1e-20)*1e4
dt_md=2.5e-15

pP=SHCPostProc.SHCPostProc(fileCompactVels,KijFilePrefix,
                           dt_md=dt_md,scaleFactor=scaleFactor,
                           LAMMPSDumpFile=fileVels,
                           widthWin=widthWin,
                           NChunks=20,chunkSize=50000,
                           backupPrefix=filePrefix,
                           LAMMPSRestartFile=restartFile,
                           reCalcVels=True,
                           reCalcFC=True)

pP.postProcess()

# Pickling the post-processing object into file
import cPickle as pickle
with open(filePrefix+'_PP.pckl','w') as f:
    pickle.dump(pP,f)
  
# Saving into numpy files 
np.save(filePrefix+'_oms.npy',pP.oms_fft)
np.save(filePrefix+'_SHC.npy',pP.SHC_smooth)

# Saving the frequencies and heat currents to file
np.savetxt(fileprefix+'_SHC.txt',np.column_stack((oms,pP.SHC_smooth)))
