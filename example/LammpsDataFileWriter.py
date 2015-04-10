# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2015

class writer(object):
    """
    Write a text file suitable for reading into LAMMPS using the read_data command.
    """

    def __init__(self,filetowrite,xs,ys,zs,**kwargs):
        """
        Positional arguments:
          fileToWrite (str): The file to which the output is written.
          xs (float array): The x-coordinates of atoms.
          ys (float array): The y-coordinates of atoms.
          zs (float array): The z-coordinates of atoms.

        Keyword arguments:
          xlo,ylo,zlo (float): Box lower bounds (default minimum coordinates minus a shift of 0.1)
          xhi,yhi,zhi (float): Box upper bounds (defaults to maximum coordinates plus a shift of 0.1)
          tags (list) : List of atom types (defaults to a list of ones)
          masses (list): List of atom masses, length must coincide with the maximum tag (defaults to one for each type)
          header (str): The header written in the first line of the output file, prints the date on default (# included automatically in the beginning.).       
        """

        self.file=filetowrite
        self.xs=xs
        self.ys=ys
        self.zs=zs

        self.xlo=None
        self.ylo=None
        self.zlo=None
        self.xhi=None
        self.yhi=None
        self.zhi=None
        self.tags=None
        self.masses=None
        self.Natoms=None
        self.atom_style=None

        # Initialize the header to current date
        from datetime import date
        today=date.today()
        self.header='File written on '+str(today.year)+'-'+str(today.month)+'-'+str(today.day)

        for key,value in kwargs.items():
            if not hasattr(self,key):
                raise ValueError, "Unknown argument to the writer."
            setattr(self,key,value)

        self.Natoms=len(xs)
        if self.Natoms!=len(ys) or self.Natoms!=len(zs):
            raise AttributeError, "xs, ys and zs must have the same size."

        if self.tags is None:
            self.tags=[1]*self.Natoms
            
        if len(self.tags)!=self.Natoms:
            raise ValueError, "The length of tag list must coincide with the length of xs,ys,zs."

        maxtag=max(self.tags)

        if self.masses is None:
            self.masses=[1]*maxtag

        if len(self.masses)!=maxtag:
            raise AttributeError, "The number of given masses must coincide with the maximum tag."

        shift=0.1
        if self.xlo is None:
            self.xlo=min(xs)-shift

        if self.xhi is None:
            self.xhi=max(xs)+shift

        if self.ylo is None:
            self.ylo=min(ys)-shift

        if self.yhi is None:
            self.yhi=max(ys)+shift

        if self.zlo is None:
            self.zlo=min(zs)-shift

        if self.zhi is None:
            self.zhi=max(zs)+shift

         
            

    def __enter__(self):
        return self

    def __exit__(self,t1,t2,t3):
        return False

    def writeToFile(self):
        print "Writing to file %s." % self.file
        with open(self.file,'w') as f:

            f.write('# '+self.header+'\n')
            f.write('%d atoms\n' % self.Natoms)
            tags=self.tags
            f.write('%d atom types\n' % max(tags))

            f.write('\n')
            
            f.write('%.2f %.2f xlo xhi\n' % (self.xlo,self.xhi))
            f.write('%.2f %.2f ylo yhi\n' % (self.ylo,self.yhi))
            f.write('%.2f %.2f zlo zhi\n' % (self.zlo,self.zhi))

            f.write('\nMasses\n\n')
            for i in range(0,max(tags)):
                f.write('%d %.3f\n' % (i+1,self.masses[i]))

            f.write('\nAtoms\n\n')

            for i in range(0,self.Natoms):
                f.write('%d %d %.5f %.5f %.5f\n' % (i+1,self.tags[i],self.xs[i],self.ys[i],self.zs[i]))

            print "Finished writing to file "+self.file+"."
