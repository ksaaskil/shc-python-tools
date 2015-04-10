# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2015

from __future__ import division
import numpy as np

# Class for post-processing
# Positional arguments:
# Compact velocity file produced with compactifyVels.cpp, written if does not exist
# KijFilePrefix

class SHCPostProc:
    '''
    Post-process the data produced using LAMMPS Molecular Dynamics simulation to calculate the spectral heat current.

    The velocities are read from the "compact" file produced with the C++-code compactify_vels.cpp from a LAMMPS dump file. If the file does not exist, it is produced by calling the binary "compactify_vels", which must be found in the environment's $PATH.

    Minimal usage:
      pP=SHCPostProc(compactVelocityFile,KijFilePrefix) # See the documentation for arguments below
      pP.postProcess() # Calculate the heat current spectrum

    Public attributes:
      SHC_smooth (numpy float array): The smoothened spectral heat current
      oms_fft (numpy float array): The angular frequencies
    '''

    def __init__(self,compactVelocityFile,KijFilePrefix,reCalcVels=False,reCalcFC=False,**args):
        '''
        Positional arguments:
          compactVelocityFile (str): The file where the velocities are read. Produced using the binary compactify_vels if the file does not exist. In this case, you must also supply the keyword argument LAMMPSDumpFile containing the velocities produced using LAMMPS.
          KijFilePrefix (str): The prefix used in trying to find the force constant matrix file KijFilePrefix.Kij.npy. If the file does not exist, the force constant calculator fcCalc is called using the keyword argument LAMMPSRestartFile (which must be supplied in this case). 

        Keyword arguments:

        '''
        self.compactVelocityFile=compactVelocityFile
        self.KijFilePrefix=KijFilePrefix
        self.Kij=None
        self.ids_L=None
        self.ids_R=None
        self.NL=None
        self.NR=None
        self.SHC_smooth=None
        self.oms_fft=None
        self.dt_md=1.0 # Default
        self.scaleFactor=1.0 # Default
        self.LAMMPSDumpFile=None
        self.LAMMPSRestartFile=None
        self.widthWin=1.0 # Default
        self.chunkSize=50000
        self.NChunks=20
        self.backupPrefix=None
        self.hstep=0.001
        self.reCalcVels=reCalcVels
        self.reCalcFC=reCalcFC

        for key,value in args.items():
            if not hasattr(self,key):
                raise ValueError, "Invalid argument " + key + " to PostProc!"
            print "Using the value "+key+"="+str(value)+"."
            setattr(self,key,value)

        import os
        if self.reCalcVels or not os.path.isfile(self.compactVelocityFile): # Check if the velocity file exists
            # Check that the LAMMPS Dump file exists
            if self.LAMMPSDumpFile is None or not os.path.isfile(self.LAMMPSDumpFile):
                raise ValueError, "You must give the LAMMPS velocity dump file as an argument to create the file "+self.compactVelocityFile+"!"
            #print self.compactVelocityFile + " does not exist, creating by reading from file " + self.LAMMPSDumpFile
            # Run the C++ script
            self._compactVels(self.LAMMPSDumpFile,self.compactVelocityFile)

        else:
            print self.compactVelocityFile + " exists, using the file for post-processing."

        # Check the force constant file
        if self.reCalcFC or not os.path.isfile(self.KijFilePrefix+'.Kij.npy'): # Check if the force constant file exists
            print "Creating file "+self.KijFilePrefix+"."
            if self.LAMMPSRestartFile is None:
                raise ValueError, "You must give the LAMMPSRestartFile as an argument so that the file "+self.KijFilePrefix+".Kij.npy can be created!"
            self._calcFC(self.KijFilePrefix,self.LAMMPSRestartFile)
        else: # Load the force constants from file
            self._loadFC(self.KijFilePrefix)


    def __enter__(self):
        return self

    def __exit__(self,t1,t2,t3):
        return False

    def _calcFC(self,fileprefix,restartfile):
        from fcCalc import fcCalc
        with fcCalc(fileprefix,restartfile) as fc:
           fc.preparelammps(pair_style='sw',pair_coeff='* * Si_vbwm.sw Si',w_interface=3.0)
           # hstep=0.001
           fc.fcCalc(self.hstep)
           fc.writeToFile()
           self.Kij=fc.Kij
           print "Size of the Kij file is (3*%d)x(3*%d)." % (np.size(self.Kij,0)/3,np.size(self.Kij,1)/3)
           self.ids_L=fc.ids_L
           self.ids_R=fc.ids_R
           self.NL=len(self.ids_L)
           print "len(ids_L)=%d" % (self.NL)
           self.NR=len(self.ids_R)
           print "len(ids_R)=%d" % (self.NR)

    def _loadFC(self,KijFilePrefix):
        print "Loading the force constants from "+ KijFilePrefix+'.Kij.npy'
        self.Kij=np.load(KijFilePrefix+'.Kij.npy')
        print "Size of the Kij file is (3*%d)x(3*%d)." % (np.size(self.Kij,0)/3,np.size(self.Kij,1)/3)
        print "Loading left interfacial atom indices from "+ KijFilePrefix+'.ids_L.npy'
        self.ids_L=np.load(KijFilePrefix+'.ids_L.npy')
        print "Loading right interfacial atom indices from "+ KijFilePrefix+'.ids_R.npy'
        self.ids_R=np.load(KijFilePrefix+'.ids_R.npy')
        self.NL=len(self.ids_L)
        print "len(ids_L)=%d" % (self.NL)
        self.NR=len(self.ids_R)
        print "len(ids_R)=%d" % (self.NR)
        if (np.size(self.Kij,0)/3!=self.NL) or (np.size(self.Kij,1)/3!=self.NR):
            raise ValueError, "Sizes in Kij and ids_L/R do not match!"

    def _compactVels(self,fileVels,finalFileVels):
        from subprocess import call
        command=["compactify_vels",fileVels,finalFileVels]
        print "Running "+" ".join(command)
        call(command)

    def postProcess(self):
        # dt_md=self.dt_md
        # widthWin=self.widthWin
        fileCompactVels=self.compactVelocityFile


        print "Reading the compact velocity file "+fileCompactVels+"."
        f=open(fileCompactVels,'r')
        s=f.readline().split()
        NAtoms=int(s[1])
        #print NAtoms, self.NL, self.NR
        if NAtoms!=self.NL+self.NR:
            raise ValueError, 'Mismatch in the numbers of atoms in the read velocity file and the used force constant file!'
        
        s=f.readline()
        print s
        s=s.split()
        sampleTimestep=int(s[1])*self.dt_md

        s=f.readline() # Atom ids:
        print s
        # Read the atom ids
        indArray=np.fromfile(f,dtype=int,count=NAtoms,sep=" ")
        # print indArray
        
        s=f.readline() # ------
        print s
  
        # Total number of degrees of freedom
        NDOF=3*(self.NL+self.NR)

        self.oms_fft=np.fft.rfftfreq(self.chunkSize,d=sampleTimestep)*2*np.pi
        Nfreqs=np.size(self.oms_fft)
        self.SHC_smooth=np.zeros(Nfreqs)

        exitFlag=False

        for k in np.arange(self.NChunks): # Start the iteration over chunks
#        for k in range(0,2): # Start the iteration over chunks
            print "Chunk %d/%d." % (k+1,self.NChunks)
            # Read the velocitites
            velArray=np.fromfile(f,dtype=np.dtype('f8'),count=self.chunkSize*NDOF,sep=" ")
            # Prepare for exit if the read size does not match the chunk size
            if np.size(velArray)==0:
                print "Finished the file, exiting."
                break
            if np.size(velArray)!=self.chunkSize*NDOF:
                # Reaching the end of file           
                self.chunkSize=int(np.size(velArray)/NDOF)             
                if k>0:
                    break
                else:
                    exitFlag=True
                    self.oms_fft=np.fft.rfftfreq(self.chunkSize,d=sampleTimestep)*2*np.pi
                    Nfreqs=np.size(self.oms_fft)
                    print "Changing chunk size to "+str(int(np.size(velArray)/NDOF))+"!"
                    
 
            # Reshape the array so that each row corresponds to different degree of freedom (e.g. particle 1, direction x etc.)
            velArray=np.reshape(velArray,(NDOF,self.chunkSize),'F')
            
            # FFT with respect to the second axis (NOTE THE USE OF RFFT)
            velFFT=np.fft.rfft(velArray,axis=1)
            # print type(velFFT[0,0])
            
            velsL=np.zeros((3*self.NL,Nfreqs),dtype=np.complex128)
            velsR=np.zeros((3*self.NR,Nfreqs),dtype=np.complex128)
            
            velsL[0::3,:]=velFFT[3*self.ids_L,:]
            velsL[1::3,:]=velFFT[3*self.ids_L+1,:]
            velsL[2::3,:]=velFFT[3*self.ids_L+2,:]

            velsR[0::3,:]=velFFT[3*self.ids_R,:]
            velsR[1::3,:]=velFFT[3*self.ids_R+1,:]
            velsR[2::3,:]=velFFT[3*self.ids_R+2,:]
            
            # Spectral heat current for the specific chunk
            SHC=np.zeros(Nfreqs)

            for ki in range(1,Nfreqs): # Skip the first one with zero frequency
                SHC[ki]=-2.0*np.imag(np.dot(velsL[:,ki],np.dot(-self.Kij,np.conj(velsR[:,ki]))))/self.oms_fft[ki]
            # SHC[0]=0.0 

            # Normalize correctly
            SHC/=(self.chunkSize*sampleTimestep)

            # Change units             
            SHC*=self.scaleFactor
            
            daniellWindow=np.ones(np.ceil(self.widthWin*2*np.pi/(self.oms_fft[1]-self.oms_fft[0])))
            daniellWindow/=np.sum(daniellWindow)

            # Smooth the value
            SHC=np.convolve(SHC,daniellWindow,'same')

            if not exitFlag:
                self.SHC_smooth=(k*self.SHC_smooth+SHC)/(k+1.0)
                if self.backupPrefix is not None:
                    np.save(self.backupPrefix+'_backup_oms.npy',self.oms_fft)
                    np.save(self.backupPrefix+'_backup_SHC.npy',self.SHC_smooth)
            elif exitFlag and k==0: # First chunk
                self.SHC_smooth=SHC
                break
            else: # Small chunk but not the first one, exit
                break

if __name__=="__main__":
    
    fileprefix='270315a'
    dt_md=2.5e-15 # Timestep used in MD, affects the frequency grid
    widthWin=0.5e12 # Width of the Daniell smoothing window in Hz
    # The velocity dump file from LAMMPS
    fileVels=fileprefix+'.vels.dat' 
    # The compactly formatted velocity file, produced using a C++ script if not found
    fileCompactVels=fileprefix+'.vels.dat.compact' 
    # Searches/saves file "fileprefix.Kij.npy"
    KijFilePrefix='270115a' 
    # Uses this restart file if the force constant file cannot be found
    restartFile=KijFilePrefix+'.quenched.restart' 
    # Correct the units, this assumes the unit of eV/(A^2)*(A/ps)^2 for the v_iK_{ij}v_j product (LAMMPS metal units)
    scaleFactor=1.602e-19/(1e-20)*1e4

    # Prepare the post-processor
    pP=SHCPostProc(fileCompactVels,KijFilePrefix,dt_md=dt_md,scaleFactor=scaleFactor,LAMMPSDumpFile=fileVels,widthWin=widthWin,NChunks=20,chunkSize=50000,backupPrefix=fileprefix)
    # Post-process
    pP.postProcess() # All variables will be contained in the object pP

    import cPickle as pickle
    with open(fileprefix+'_PP.pckl','w') as f:
        pickle.dump(pP,f)
    
    np.save(fileprefix+'_oms.npy',pP.oms_fft)
    np.save(fileprefix+'_SHC.npy',pP.SHC_smooth)

    # Plotting if available
    # import matplotlib.pylab as plt
    # plt.plot(pP.oms_fft/(2*np.pi*1.0e12),pP.SHC_smooth)
    # plt.xlabel('Frequency (THz)')
    # plt.ylabel('Spectral current')
    # plt.savefig(fileprefix+'_SHC.eps')
    
      
        
 
