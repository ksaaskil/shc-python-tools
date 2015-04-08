
import numpy as np

class SHCPostProc:

    def __init__(self,fileprefix,KijFilePrefix):
        self.dt_md=dt_md
        self.fileprefix=fileprefix
        self.Kij=None
        self.ids_L=None
        self.ids_R=None
        self.NL=None
        self.NR=None
        self.SHC_smooth=None
        self.oms_fft=None

    def __enter__(self):
        return self

    def __exit__(self,t1,t2,t3):
        return False

    def calcFC(self,fileprefix,restartfile):
        from fcCalc import fcCalc
        with fcCalc(fileprefix,restartfile) as fc:
           fc.preparelammps(pair_style='sw',pair_coeff='* * Si_vbwm.sw Si',w_interface=3.0)
           hstep=0.001
           fc.fcCalc(hstep)
           fc.writeToFile()
           self.Kij=fc.Kij
           self.ids_L=fc.ids_L
           self.ids_R=fc.ids_R
           self.NL=len(self.ids_L)
           print "len(ids_L)=%d" % (self.NL)
           self.NR=len(self.ids_R)
           print "len(ids_R)=%d" % (self.NR)

    def loadFC(self,KijFilePrefix):
        self.Kij=np.load(KijFilePrefix+'.Kij.npy')
        self.ids_L=np.load(KijFilePrefix+'.ids_L.npy')
        self.ids_R=np.load(KijFilePrefix+'.ids_R.npy')
        self.NL=len(self.ids_L)
        print "len(ids_L)=%d" % (self.NL)
        self.NR=len(self.ids_R)
        print "len(ids_R)=%d" % (self.NR)

    def compactVels(self,fileVels,finalFileVels):
        from subprocess import call
        command=["./compactify_vels",fileVels,finalFileVels]
        print "Running "+" ".join(command)
        call(command)

    def postProcess(self,fileCompactVels,dt_md,widthWin):
        f=open(fileCompactVels,'r')
        s=f.readline().split()
        NAtoms=int(s[1])
        print NAtoms, self.NL, self.NR
        assert NAtoms==self.NL+self.NR, 'Mismatch in the numbers of atoms in the read velocity file and the used force constant file!'
        s=f.readline().split()
        sampleTimestep=int(s[1])*dt_md

        f.readline() # Atom ids:

        indArray=np.fromfile(f,dtype=int,count=NAtoms,sep=" ")
        #print indArray
        f.readline() # ------

        NChunks=100
        chunkSize=10000

        NDOF=3*(self.NL+self.NR)

        self.SHC_smooth=np.zeros(chunkSize)
        oms_fft=np.arange(0,chunkSize)/(chunkSize*sampleTimestep)*2*np.pi
        self.oms_fft=oms_fft
        exitFlag=False
        for k in np.arange(NChunks): # Start the iteration over chunks
#        for k in range(0,2): # Start the iteration over chunks
            print "Chunk %d/%d." % (k+1,NChunks)
            # Read the velocitites
            velArray=np.fromfile(f,dtype=np.dtype('f8'),count=chunkSize*NDOF,sep=" ")
            # Prepare for exit if the read size does not match the chunk size
            if not np.size(velArray)==chunkSize*NDOF:
                # Reaching the end of file
                print "Changing chunk size to "+str(np.size(velArray)/NDOF)+"!"
                chunkSize=np.size(velArray)/NDOF
                exitFlag=True

            # Reshape the array   
            velArray=np.reshape(velArray,(NDOF,chunkSize))
            
            # FFT with respect to the second axis
            velFFT=np.fft.fft(velArray,axis=-1)
            
            velsL=np.zeros((3*self.NL,chunkSize),dtype=np.complex64)
            velsR=np.zeros((3*self.NR,chunkSize),dtype=np.complex64)
            
            velsL[0::3,:]=velFFT[3*self.ids_L,:]
            velsL[1::3,:]=velFFT[3*self.ids_L+1,:]
            velsL[2::3,:]=velFFT[3*self.ids_L+2,:]

            velsR[0::3,:]=velFFT[3*self.ids_R,:]
            velsR[1::3,:]=velFFT[3*self.ids_R+1,:]
            velsR[2::3,:]=velFFT[3*self.ids_R+2,:]
            
            # Spectral heat current for the specific chunk
            SHC=np.zeros(chunkSize)

            for ki in range(1,chunkSize): # Skip the first one with zero frequency
                SHC[ki]=-2.0*np.imag(np.dot(velsL[:,ki],np.dot(self.Kij,np.conj(velsR[:,ki]))))/oms_fft[ki]
            SHC[0]=0.0 

            # Normalize correctly
            SHC/=(chunkSize*sampleTimestep)
            
            daniellWindow=np.ones(np.ceil(widthWin*2*np.pi/(oms_fft[2]-oms_fft[1])))
            daniellWindow/=np.sum(daniellWindow)

            # Smooth the value
            SHC=np.convolve(SHC,daniellWindow,'same')

            if not exitFlag:
                self.SHC_smooth=(k*self.SHC_smooth+SHC)/(k+1.0)
            elif exitFlag and k==0: # First chunk
                self.SHC_smooth=SHC
                break
            else: # Small chunk but not the first one, exit
                break

if __name__=="__main__":
    
    fileprefix='270315a'
    dt_md=0.25e-15 # Timestep used in MD, affects the frequency grid
    widthWin=1e12 # Width of the Daniell smoothing window in Hz
    fileVels=fileprefix+'.vels.dat' # The velocity dump file from LAMMPS
    fileCompactVels=fileprefix+'.vels.dat.compact' # The compactly formatted velocity file, produced using a C++ script if not found
    KijFilePrefix='270115a' # Searches/saves file "fileprefix.Kij.npy"
    restartFile=KijFilePrefix+'.quenched.restart' # Uses this restart file if the force constant file cannot be found

    with SHCPostProc(fileprefix,KijFilePrefix) as pP:
        import os
        if not os.path.isfile(fileCompactVels): # Check if the velocity file exists
            # Run the C++ script
            pP.compactVels(fileVels,fileCompactVels)

        if not os.path.isfile(KijFilePrefix+'.Kij.npy'): # Check if the force constant file exists
            assert restartFile is not None, "You must give the restart file if the KijFile does not exist!"
            pP.calcFC(KijFilePrefix,restartFile)
        else: # Load the force constants from file
            pP.loadFC(KijFilePrefix)

        pP.postProcess(fileCompactVels,dt_md,widthWin)
        np.save(fileprefix+'_oms.npy',pP.oms_fft)
        np.save(fileprefix+'_SHC.npy',pP.SHC_smooth)
 
