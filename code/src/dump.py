class dump:
    """
    Read LAMMPS dump file.

    Usage:
        d = dump("dump.atom")                                 # Read one dump file.
        alltime = d.time()                                    # Get all timesteps in the dump file.
        d.read(N)                                             # Read the Nth timestep (Nmin = 0, Nmax = number of snapshots - 1).
        d.wrap()                                              # Calculate wrapped unscaled coordinate.
        d.unwrap()                                            # Calculate unwrapped unscaled coordinate.
        X, Y, Z = d.map("x","y","z")                          # Find the index for required properties.
        idmapping = d.idtoindex()                             # Map AtomID to Atom index for current step.
        timestep,natoms,box,props,atoms = d.snapshot()        # Get all data from current snapshot.
        d.printatom(AtomID)                                   # Print the information of one atom with required AtomID at current step.
        d.changeatoms(new_atoms,new_props)                    # Change and replace old atoms data by new atoms data (new_atoms) and new property list (new_props).
        d.write("new_dump.atom",properties,"w")               # Create a new file "new_dump.atom" and write the current snapshot to it, including properties in `properties`.
        d.write("new_dump.atom",properties,"a")               # Append the current snapshot to file "new_dump.atom", including properties in `properties`. File "new_dump.atom" will be created if it does not exist.
    """
    def __init__(self, filename):
        """
        dump(filname)
            Read one dump file.

        Parameters:
            filename: string
                Name of dump file.

        Return:
            None.
        """
        self._dump = _ReadDump(filename)
        self._step = _STEP()
        self._idmapping = _VecInt()

    def time(self):
        """
        time()
            Get all timesteps in the dump file.

        Parameters:
            None.

        Return:
            timesteps: python list (int)
                All the timesteps in the dump file.
        """
        tmp_VecInt = _VecInt()
        ret = self._dump.GetTimesteps(tmp_VecInt)
        return list(tmp_VecInt)

    def read(self, N):
        """
        read(N)
            Read the Nth timestep (Nmin = 0, Nmax = number of snapshots - 1).

        Parameters:
            N: int
                Index of timestep.

        Return
            None.

        Example:
            d = dump("dump.atom")
            alltime = d.time()
            for i in range(len(alltime)):
                d.read(i)
                d.wrap()
                timestep,natoms,box,props,atoms = d.snapshot()
                # do something
        """
        ret = self._dump.Read_Frame(N, self._step)
        if ret == 0:
            ### build the AtomID --> Atom index mapping
            self._dump.IDMapping(self._idmapping, self._step)
            print("Current step's timestep = %d" % self._step.Timestep)

    def wrap(self):
        """
        wrap()
            Calculate wrapped unscaled coordinate.
            This will change the internal coordinates data, and change property names for coordinates to "x","y",and "z".
            Changed data can be accessed by snapshot()

        Parameters:
            None.

        Return:
            None.
        """
        ret = self._dump.CalcWrappedCoord(self._step)
        if ret == 0:
            print("Timestep = %d is wrapped." % self._step.Timestep)
        else:
            print("Failed to wrap timestep = %d." % self._step.Timestep)

    def unwrap(self):
        """
        unwrap()
            Calculate unwrapped unscaled coordinate.
            This will change the internal coordinates data, and change property names for coordinates to "xu","yu",and "zu".
            Changed data can be accessed by snapshot()

        Parameters:
            None.

        Return:
            None.
        """
        ret = self._dump.CalcUnwrappedCoord(self._step)
        if ret == 0:
            print("Timestep = %d is unwrapped." % self._step.Timestep)
        else:
            print("Failed to unwrap timestep = %d." % self._step.Timestep)

    def map(self,*props):
        """
        map(prop1,prop2,...)
            Find the index for required properties.

        Parameters:
            prop1, prop2, ...: string
                Property names.

        Return:
            index1, index2, ...: int
                Property indices. If required property is not found, index will be assigned to -1

        Example:
            import numpy as np
            d = dump("dump.atom")
            d.read(0)
            d.wrap()
            X,Y,Z,IX,IY,IZ = d.map("x","y","z","ix","iy","iz")
            timestep,natoms,box,props,atoms = d.snapshot()
            atoms = np.transpose(atoms)
            print(atoms[X],atom[IX])    # will print wrapped unscaled x coordinate and box image ix of all atoms
        """
        if len(props) > 1:
            tmp_list = [-1 for i in range(len(props))]
            for i in range(0,len(props)):
                ret, tmp_list[i] = self._dump.PropertyMapping(props[i],self._step)
                if tmp_list[i] == -1:
                    print("Warning: dump map() cannot find required property %s" % props[i])
            return tmp_list
        elif len(props) == 1:
            ret, tmp = self._dump.PropertyMapping(props[0],self._step)
            if tmp == -1:
                print("Warning: dump map() cannot find required property %s" % props)
            return tmp
        else:
            print("Invalid syntax for map()")

    def idtoindex(self):
        """
        idtoindex():
            Map AtomID to Atom index for current step.

        Parameters:
            None.

        Return:
            idmapping: numpy 1d array (int)
                Mapping of AtomID --> Atom index.
                i.e., idmapping[AtomID] = Atom index
                If the AtomID is not found, Atom index will be set to (number of atoms in current step + 1)

        Example:
            import numpy as np
            d = dump("dump.atom")
            d.read(0)
            d.wrap()
            timestep,natoms,box,props,atoms = d.snapshot()
            idmapping = d.idtoindex()
            print(atoms[idmapping[10]])    # will print properties of the atom whose AtomID = 10
        """
        tmp_idmapping = _ConvertToNumpy_VecInt(self._idmapping)
        ### since numpy take negative index (count from the end), we should change -1 (indicate AtomID not found) to (number of atoms in current step + 1), so numpy can throw an error
        mask = (tmp_idmapping == -1)
        tmp_idmapping = tmp_idmapping + mask * (self._step.NumOfAtoms+ 1)
        return tmp_idmapping

    def snapshot(self):
        """
        snapshot()
            Get all data from current snapshot.

        Parameters:
            None.

        Return:
            timestep: int
                Current timestep.
            natoms: int
                Number of atoms in current snapshot.
            box: python list (double)
                A list of box edge of current snapshot: [xlo, ylo, zlo, xhi, yhi, zhi] (non triclinic box)
                or [xlo_bound, ylo_bound, zlo_bound, xhi_bound, yhi_bound, zhi_bound, xy, xz, yz] (triclinic box)
            props: python list (string)
                A list of property names in current snapshot.
            atoms: numpy 2d array (double)
                A numpy 2d array in shape (natoms, len(props)). Rows are atoms, columns are properties. Not sorted.

        Example:
            import numpy as np
            d = dump("dump.atom")
            d.read(0)
            d.wrap()
            X,Y,Z,IX,IY,IZ = d.map("x","y","z","ix","iy","iz")
            timestep,natoms,box,props,atoms = d.snapshot()
            print(atoms[0])    # will print all properties of one atoms
            atoms = np.transpose(atoms)
            print(atoms[X],atom[IX])    # will print wrapped unscaled x coordinate and box image ix of all atoms
        """
        timestep = self._step.Timestep
        natoms = self._step.NumOfAtoms
        if self._step.Box.triclinic:
            box = list(_VecDouble([self._step.Box.xlo_bound,self._step.Box.ylo_bound,self._step.Box.zlo_bound,
                                   self._step.Box.xhi_bound,self._step.Box.yhi_bound,self._step.Box.zhi_bound,
                                   self._step.Box.xy,self._step.Box.xz,self._step.Box.yz]))
        else:
            box = list(_VecDouble([self._step.Box.xlo,self._step.Box.ylo,self._step.Box.zlo,self._step.Box.xhi,self._step.Box.yhi,self._step.Box.zhi]))
        props = list(self._step.Properties)
        atoms = _ConvertToNumpy_VecVecDouble(self._step.Atoms)
        if self._step.wrapped_flag == 2:
            print("Warning: The coordinate is not really wrapped!")
        return timestep,natoms,box,props,atoms

    def printatom(self,AtomID):
        """
        printatom(AtomID):
            Print the information of one atom with required AtomID at current step.

        Patameters:
            AtomID: int
                Required AtomID.

        Return:
            None.
        """
        print("timestep = %d" % self._step.Timestep)
        if AtomID >= self._idmapping.size():
            print("Required AtomID = %d not found." % AtomID)
        elif self._idmapping[AtomID] == -1:
            print("Required AtomID = %d not found." % AtomID)
        else:
            for i in range(self._step.Properties.size()):
                print("'%s' : %g" % (self._step.Properties[i], self._step.Atoms[self._idmapping[AtomID]][i]))

    def changeatoms(self,new_atoms,new_props):
        """
        changeatoms(new_atoms)
            Change atoms data. Overwrite current atoms data by a new numpy 2d array (new_atoms), which
            include different properties (new_props) for each atom.

        Parameters:
            new_atoms: numpy 2d array (double)
                A numpy 2d array with atoms data.
            new_props: python list (string)
                A list of property names in new_atoms.

        Return:
            None.
        """
        if new_atoms.shape[1] != len(new_props):
            print("The size of 2d array with new atoms data does not match the number of provided new properties.")
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Failed to change atoms data.")
        else:
            _ConvertToC_2dDouble(self._step.Atoms, new_atoms)
            self._step.NumOfAtoms = new_atoms.shape[0]
            self._step.Properties.resize(len(new_props))
            self._step.ID = self._step.TYPE = self._step.MOL = self._step.X = self._step.Y = self._step.Z = self._step.IX = self._step.IY = self._step.IZ \
            = self._step.VX = self._step.VY = self._step.VZ = self._step.FX = self._step.FY = self._step.FZ = self._step.Q = -1
            for i in range(len(new_props)):
                self._step.Properties[i] = new_props[i]
                if new_props[i] == "id":
                    self._step.ID = i
                elif new_props[i] == "type":
                    self._step.TYPE = i
                elif new_props[i] == "mol":
                    self._step.MOL = i
                elif new_props[i] == "x" or new_props[i] == "xu" or new_props[i] == "xs" or new_props[i] == "xsu":
                    self._step.X = i
                elif new_props[i] == "y" or new_props[i] == "yu" or new_props[i] == "ys" or new_props[i] == "ysu":
                    self._step.Y = i
                elif new_props[i] == "z" or new_props[i] == "zu" or new_props[i] == "zs" or new_props[i] == "zsu":
                    self._step.Z = i
                elif new_props[i] == "ix":
                    self._step.IX = i
                elif new_props[i] == "iy":
                    self._step.IY = i
                elif new_props[i] == "iz":
                    self._step.IZ = i
                elif new_props[i] == "vx":
                    self._step.VX = i
                elif new_props[i] == "vy":
                    self._step.VY = i
                elif new_props[i] == "vz":
                    self._step.VZ = i
                elif new_props[i] == "fx":
                    self._step.FX = i
                elif new_props[i] == "fy":
                    self._step.FY = i
                elif new_props[i] == "fz":
                    self._step.FZ = i
                elif new_props[i] == "q":
                    self._step.Q = i
            if 'x' in new_props and 'y' in new_props and 'z' in new_props:
                self._step.wrapped_flag = 1
                self._step.scaled_flag = 0
            elif 'xu' in new_props and 'yu' in new_props and 'zu' in new_props:
                self._step.wrapped_flag = 0
                self._step.scaled_flag = 0
            elif 'xs' in new_props and 'ys' in new_props and 'zs' in new_props:
                self._step.wrapped_flag = 1
                self._step.scaled_flag = 1
            elif 'xsu' in new_props and 'ysu' in new_props and 'zsu' in new_props:
                self._step.wrapped_flag = 0
                self._step.scaled_flag = 1
            else:
                print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Inconsistant coordinate type.")
            ### build the AtomID --> Atom index mapping
            self._dump.IDMapping(self._idmapping, self._step)

    def write(self,filename,properties,mode="a"):
        """
        write(filename, properties, mode):
            Write the current snapshot (read by calling member function `read`) to the file `filename`, including properties listed in `properties`.

        Parameters:
            filename: string
                Name of the output LAMMPS dump file
            properties: python list (string)
                A list of properties to be written to the output file.
            mode: string
                Writing mode. Should be either "w" (for creating a new file or overwriting the file that already exists) or "a" (for appending to the file).
                If "a" is given but the file `filename` does not exist, a new file `filename` will be created.

        Return:
            None.
        """
        writer = _WriteDump(filename,properties,mode)
        ret = writer.Write(self._step)
        if ret == 0:
            print("timestep = %d is written to file \"%s\"" % (self._step.Timestep, filename))
        else:
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Failed to write dump into file.")
