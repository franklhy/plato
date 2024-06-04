class data:
    """
    Read and write LAMMPS data file. Trajectory from LAMMPS dump file can also be loaded.

    Usage:
        d = data()                                         # Declare an instance of data class
        d.print_info(True)                                 # whether to print information when operating on files and data. False by default.
        d.read("data.atom","molecular")                    # Read the data file "data.atom" according to the format of the atom type "molecular"
        d.read("data.atom")                                # Read the data file "data.atom". The format will be determined according to the data file.
        d.set_traj_file("dump.atom")                       # Set a dump file as the trajactory while keeping the topology (bond,angle,...) information
        d.set_traj_files("dump1", "dump2", ...)            # Set multiple dump files as the trajactory
        alltime = d.traj_time()                            # Return all timesteps in the trajactory dump file
        d.load_traj(i)                                     # Read the ith frame from the trajactory dump file
        d.wrap()                                           # Calculate wrapped unscaled coordinate.
        d.unwrap()                                         # Calculate unwrapped unscaled coordinate.
        X, Y, Z = d.map("x","y","z")                       # Find the index for required properties.
        idmapping_atom = d.idtoindex_atom()                # Map AtomID to Atom index for current step.
        timestep,natoms,box,props,atoms = d.snapshot()     # Get all data (except topology data) from current snapshot.
        masses = d.masses()                                # Get the mass of each atom type.
        bonds,angles,dihedrals,impropers = d.topology()    # Get all topology data from current snapshot.
        d.change_timestep(new_timestep)                    # Change and replace old timestep.
        d.change_box(new_box,triclinic)                    # Change and replace old simulation box dimensions by new simulation box (new_box)
        d.change_masses(new_masses)                        # Change and replace old masses by the new ones (new_masses).
        d.change_atoms(new_atoms,new_props)                # Change and replace old atoms data by new atoms data (new_atoms) and new property list (new_props).
        d.change_bonds(new_bonds)                          # Change and replace old bonds data by new bonds data (new_bonds).
        d.change_angles(new_angles)                        # Change and replace old angles data by new angles data (new_angles).
        d.change_dihedrals(new_dihedrals)                  # Change and replace old dihedrals data by new dihedrals data (new_dihedrals).
        d.change_impropers(new_impropers)                  # Change and replace old impropers data by new impropers data (new_impropers).
        d.reset_atomid()                                   # Reset the atom ID of atoms, so they are numbered contiguously.
        d.flip()                                           # Flip the simulation box is the tilt factor is larger than 0.5
        d.write_data("new_data.atom","molecular")          # Write the current step into LAMMPS data file "new_data.atom" in "molecular" format.
        d.write_dump("new_dump.atom", props, "a")          # Write the current step into LAMMPS dump file "new_dump.atom", including properties in "props", with appending mode ("a")
    """
    def __init__(self):
        """
        data()

        Parameters:
            None

        Return:
            None.
        """
        self._print_flag = False
        self._data = _ReadData()
        self._dumps = None
        self._step = _STEP()
        self._timesteps = None
        self._frame_dict = None    # key: frame index; value: (dumpfile index, frame index of the dumpfile, timestep)
        self._idmapping_atom = _VecInt()
        self._writer = _WriteData()

    def print_info(self, print_flag=False):
        """
        print_info(print_flag)
            Whether to print information when operating on files and data.

        Parameters:
            print_flag: bool
                Print flag. Default = False.
        """
        self._print_flag = print_flag

    def read(self,filename,atomstyle=""):
        """
        read(filename,atomstyle="")
            Read the data file (filename) according to the format of an atom type (atomstyle)

        Parameters:
            filename: string
                Name of LAMMPS data file.
            atomstyle: string
                Name of atom style ("atomic", "molecular", "full").
                It can also be an empty string "". In this case, the Atoms section of data file should start with something like "Atom # molecular",
                so the code can recognize the atom style of the data file.
                By default, atomstyle = ""

        Return
            None.

        Example:
            d = data()
            d.read("data.atom")
            d.unwrap()
            timestep,natoms,box,props,atoms = d.snapshot()
            # do something
        """
        self._step = _STEP()
        ret = self._data.Read_File(filename, self._step, atomstyle)
        if ret == 0:
            ### build the AtomID --> Atom index mapping
            self._data.AtomIDMapping(self._idmapping_atom, self._step)
            print("Current step's timestep = %d" % self._step.Timestep)
        else:
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Fail to read.")

    def set_traj_files(self, *dumpfiles, overwrite=True, require_same_dump_freq=True):
        """
        set_traj_files(file1, file2, file3, ...)
            Set multiple dump files as the trajactory while keeping the topology (bond,angle,...) information read from data file.

        Parameters:
            dumpfile1, dumpfile2, ...: string
                Name of LAMMPS dump file(s).

            overwrite: bool (default=True)
                When multiple dump files are read, some dump files might contain snapshots with the same timestep as the snapshots from 
                previous dump files. This could happen if simulations are restarted from a specific snapshot from previous dump file.
                It decides whether to overwrite the old snapshots from previous dump files with the new ones if they share the same timestep.

            require_same_dump_freq: bool (default=True)
                When multiple dump files are read, they might not be generated with the same dump frequency. This could cause trouble when
                calculating time dependent properties. When set to True, it will check if all dump files share the same dump frequency.

        Return:
            None.
        """
        def check_atom_id_and_type(tmp_dump, dumpfile):
            """
            check atom id and type, they should be identical to those in data file (stored in self._step)
            """
            tmp_step = _STEP()
            ret = tmp_dump.Read_Frame(0, tmp_step)
            if ret:
                raise RuntimeError("Fail to read the 1st frame of %s." % dumpfile)
            tmp_atoms = _ConvertToNumpy_VecVecDouble(tmp_step.Atoms)
            ret, tmpID = tmp_dump.PropertyMapping("id", tmp_step)
            if ret:
                raise RuntimeError("Fail to map property \"id\" from trajectory dump file %s." % dumpfile)
            ret, tmpTYPE = tmp_dump.PropertyMapping("type", tmp_step)
            if ret:
                raise RuntimeError("Fail to map property \"type\" from trajectory dump file %s." % dumpfile)
            tmp_orderid = np.argsort(tmp_atoms[:,tmpID])
            tmpaid = tmp_atoms[tmp_orderid][:,tmpID]
            tmpatype = tmp_atoms[tmp_orderid][:,tmpTYPE]

            atoms = _ConvertToNumpy_VecVecDouble(self._step.Atoms)
            ret, ID = self._data.PropertyMapping("id", self._step)
            if ret:
                raise RuntimeError("Fail to map property \"id\" from data file.")
            ret, TYPE = self._data.PropertyMapping("type", self._step)
            if ret:
                raise RuntimeError("Fail to map property \"type\" from data file.")
            orderid = np.argsort(atoms[:,ID])
            aid = atoms[orderid][:,ID]
            atype = atoms[orderid][:,TYPE]

            if not np.array_equal(aid, tmpaid) or not np.array_equal(atype, tmpatype):
                raise RuntimeError("The trajectory dump file %s has different atoms then data file." % dumpfile)

        tmp_VecInt = _VecInt()
        if len(dumpfiles) > 0:
            self._dumps = [_ReadDump(dumpfile) for dumpfile in dumpfiles]

            # the first dump file
            check_atom_id_and_type(self._dumps[0], dumpfiles[0])
            _ = self._dumps[0].GetTimesteps(tmp_VecInt)
            self._timesteps = np.array(list(tmp_VecInt))
            self._frame_dict = {ii: (0, ii, self._timesteps[ii]) for ii in range(len(self._timesteps))}

            # the remaining dump files
            for i in range(1,len(dumpfiles)):
                check_atom_id_and_type(self._dumps[i], dumpfiles[i])
                _ = self._dumps[i].GetTimesteps(tmp_VecInt)
                tmp_timesteps = np.array(list(tmp_VecInt))
                if require_same_dump_freq and tmp_timesteps[1] - tmp_timesteps[0] != self._timesteps[1] - self._timesteps[0]:
                    raise RuntimeError("Dump files don't share the same dump frequency.")
                ### if the first timestep is 0, it assumes the first timestep is the same as the last timestep from the previous dump file
                if tmp_timesteps[0] == 0:
                    tmp_timesteps += self._timesteps[-1]
                for ii in range(len(tmp_timesteps)):
                    if tmp_timesteps[ii] in self._timesteps:
                        if overwrite:
                            ind = np.where(self._timesteps == tmp_timesteps[ii])[0][0]
                            print("timestep %d: replace dump file no.%d frame no.%d -> dump file no.%d frame no.%d" % 
                                  (tmp_timesteps[ii], self._frame_dict[ind][0]+1, self._frame_dict[ind][1]+1, i+1, ii+1))
                            self._frame_dict[ind] = (i, ii, tmp_timesteps[ii])
                    else:
                        self._timesteps = np.append(self._timesteps, tmp_timesteps[ii])
                        self._frame_dict[len(self._timesteps)-1] = (i, ii, tmp_timesteps[ii])
        else:
            raise RuntimeError("no dump file is provided.")

    def set_traj_file(self, dumpfile):
        """
        set_traj_file(dumpfile):
            Set a dump file (dumpfile) as the trajactory while keeping the topology (bond,angle,...) information read from data file.

        Parameters:
            dumpfile: string
                Name of LAMMPS dump file.

        Returun:
            None.
        """
        self.set_traj_files(dumpfile)

    def traj_time(self):
        """
        traj_time()
            Return all timesteps in the trajactory dump file.

        Parameters:
            None.

        Return:
            timesteps: python list (int)
                All timesteps in the trajactory dump file.
        """
        if self._dumps is None:
            raise RuntimeError("Please set trajectory dump file by member function set_traj_file(dumpfile) first.")
        return list(self._timesteps)

    def load_traj(self, N):
        """
        Read the Nth frame from the trajactory dump file(s).

        Parameters:
            N, int
                Index of a timestep.

        Return:
            None.

        Example:
            d = data()
            d.read("atom.dat")
            d.set_traj_file("atom.dump")
            alltime = d.traj_time()
            for i in range(len(alltime)):
                d.load_traj(i)
        """
        if self._dumps is None:
            raise RuntimeError("Please set trajectory dump file by member function set_traj_file(dumpfile) first.")
        dumpid, frameid, timestep = self._frame_dict[N]
        ret = self._dumps[dumpid].Read_Frame(frameid, self._step)
        if ret == 0:
            self._step.Timestep = int(timestep)
            ### build the AtomID --> Atom index mapping
            self._dumps[dumpid].IDMapping(self._idmapping_atom, self._step)
            if self._print_flag:
                print("Current step's timestep = %d" % self._step.Timestep)
        else:
            raise RuntimeError("Fail to read dump file.")

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
        ret = self._data.CalcWrappedCoord(self._step)
        if ret == 0:
            if self._print_flag:
                print("Timestep = %d is wrapped." % self._step.Timestep)
        else:
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Failed to wrap timestep = %d." % self._step.Timestep)

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
        ret = self._data.CalcUnwrappedCoord(self._step)
        if ret == 0:
            if self._print_flag:
                print("Timestep = %d is unwrapped." % self._step.Timestep)
        else:
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Failed to unwrap timestep = %d." % self._step.Timestep)

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
            d = data()
            d.read("data.atom","molecular")
            d.wrap()
            X,Y,Z,IX,IY,IZ = d.map("x","y","z","ix","iy","iz")
            timestep,natoms,box,props,atoms = d.snapshot()
            atoms = np.transpose(atoms)
            print(atoms[X],atom[IX])    # will print wrapped unscaled x coordinate and box image ix of all atoms
        """
        if len(props) > 1:
            tmp_list = [-1 for i in range(len(props))]
            for i in range(0,len(props)):
                ret, tmp_list[i] = self._data.PropertyMapping(props[i],self._step)
                if tmp_list[i] == -1:
                    print("\x1b[31m" + "Warning: " + "\x1b[0m" + "data map() cannot find required property %s" % props[i])
            return tmp_list
        elif len(props) == 1:
            ret, tmp = self._data.PropertyMapping(props[0],self._step)
            if tmp == -1:
                print("\x1b[31m" + "Warning: " + "\x1b[0m" + "data map() cannot find required property %s" % props)
            return tmp
        else:
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Invalid syntax for map()")

    def idtoindex_atom(self):
        """
        idtoindex_atom():
            Map AtomID to Atom index for current step.

        Parameters:
            None.

        Return:
            idmapping_atom: numpy 1d array (int)
                Mapping of AtomID --> Atom index.
                i.e., idmapping_atom[AtomID] = Atom index
                If the AtomID is not found, Atom index will be set to (number of atoms in current step + 1)

        Example:
            import numpy as np
            d = data()
            d.read("data.atom", "molecular")
            d.unwrap()
            timestep,natoms,box,props,atoms = d.snapshot()
            idmapping_atom = d.idtoindex_atom()
            print(atoms[idmapping_atom[10]])    # will print properties of the atom whose AtomID = 10

        Example:
            ### Order atoms by their ID:
            d = data()
            d.read("data.atom", "molecular")
            d.unwrap()
            timestep,natoms,box,props,atoms = d.snapshot()
            idmapping_atom = d.idtoindex_atom()
            order_by_id = idmapping_atom[1:]
            atoms = atoms[order_by_id]
            ### now rows of atoms are ordered by atom ID 
        """
        tmp_idmapping = _ConvertToNumpy_VecInt(self._idmapping_atom)
        ### since numpy take negative index (count from the end), we should change -1 (indicate AtomID not found) to tmp_idmapping.size+1, so numpy can throw an error
        mask = (tmp_idmapping == -1)
        tmp_idmapping = tmp_idmapping + mask * (self._step.NumOfAtoms + 1)
        return tmp_idmapping

    def snapshot(self):
        """
        snapshot()
            Get all data (except topology data) from current snapshot.

        Parameters:
            None.

        Return:
            timestep: int
                Current timestep.
            natoms: int
                Number of atoms in current snapshot.
            box: python list (double)
                A list of box edge of current snapshot: [xlo, ylo, zlo, xhi, yhi, zhi] (non triclinic box)
                or [xlo, ylo, zlo, xhi, yhi, zhi, xy, xz, yz] (triclinc box)
            props: python list (string)
                A list of property names in current snapshot.
            atoms: numpy 2d array (double)
                A numpy 2d array in shape (natoms, len(props)). Rows are atoms, columns are properties. Not sorted.

        Example:
            import numpy as np
            d = data()
            d.read("data.atom","molecular")
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
            box = list(_VecDouble([self._step.Box.xlo,self._step.Box.ylo,self._step.Box.zlo,
                                   self._step.Box.xhi,self._step.Box.yhi,self._step.Box.zhi,
                                   self._step.Box.xy,self._step.Box.xz,self._step.Box.yz]))
        else:
            box = list(_VecDouble([self._step.Box.xlo,self._step.Box.ylo,self._step.Box.zlo,self._step.Box.xhi,self._step.Box.yhi,self._step.Box.zhi]))
        props = list(self._step.Properties)
        atoms = _ConvertToNumpy_VecVecDouble(self._step.Atoms)
        if self._step.wrapped_flag == 2:
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "The coordinate is not really wrapped!")
        return timestep,natoms,box,props,atoms

    def masses(self):
        """
        masses()
            Get the mass of each atom type.

        Parameters:
            None.

        Return:
            masses: python list (double)
                A python list with (Ntype + 1) elements, where Ntype is the number of atom types. masses[i] is the mass of type i. masses[0] = 0 is redundant.
                If masses is not included in the data file being read, then the return will be a python list with only one element masses[0] = 0.
        """
        masses = list(self._step.Masses)
        if len(masses) == 0:
            print("No masses information in the data file.")
            masses = [0]
        return masses

    def topology(self):
        """
        topology()
            Get all topology data from current snapshot.

        Parameters:
            None.

        Return:
            bonds: numpy 2d array (int)
                A numpy 2d array in shape (number of bonds, 4). Rows are bonds, columns are BondID, BondType, BondAtom1ID, BondAtom2ID
            angles: numpy 2d array (int)
                A numpy 2d array in shape (number of angles, 5). Rows are angles, columns are AngleID, AngleType, Atom1~3's ID
            dihedrals: numpy 2d array (int)
                A numpy 2d array in shape (number of dihedrals, 6). Rows are dihedrals, columns are DihedralID, DihedralType, Atom1~4's ID
            impropers: numpy 2d array (int)
                A numpy 2d array in shape (number of impropers, 6). Rows are impropers, columns are ImproperID, ImproperType, Atom1~4's ID

        Example:
            d = data()
            d.read("data.atom", "molecular")
            d.unwrap()
            bonds,angles,dihedrals,impropers = d.topology()
        """ 
        bonds = _ConvertToNumpy_VecVecInt(self._step.Bonds)
        angles = _ConvertToNumpy_VecVecInt(self._step.Angles)
        dihedrals = _ConvertToNumpy_VecVecInt(self._step.Dihedrals)
        impropers = _ConvertToNumpy_VecVecInt(self._step.Impropers)
        return bonds,angles,dihedrals,impropers

    def change_timestep(self, new_timestep):
        """
        change_timestep(new_timestep):
            Change timestep. Overwrite the current timestep.

        Parameters:
            new_timestep: int
                New timestep

        Return:
            None.
        """
        self._step.Timestep = new_timestep

    def change_box(self, new_box, triclinic=False):
        """
        change_box(new_box):
            Change simulation box. Overwrite current dimensions of simulation box by a python list with
            6 elements for non triclinic box ([xlo, ylo, zlo, xhi, yhi, zhi]) or 9 elements for triclinic
            box ([xlo, ylo, zlo, xhi, yhi, zhi, xy, xz, yz]).

        Parameters:
            new_box: python list (double) or numpy 1d array (double)
                A list or numpy 1d array recording [xlo, ylo, zlo, xhi, yhi, zhi] or 
                [xlo, ylo, zlo, xhi, yhi, zhi, xy, xz, yz]
            triclinic: bool (default=False)
                Whether the new box is triclinic

        Return:
            None.
        """
        if not triclinic and len(new_box) != 6:
            print("There should be 6 elements for new_box if not triclinic, which are [xlo, ylo, zlo, xhi, yhi, zhi].")
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Failed to change simluation box.")
        elif triclinic and len(new_box) != 9:
            print("There should be 6 elements for new_box if triclinic, which are [xlo, ylo, zlo, xhi, yhi, zhi, xy, xz, yz].")
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Failed to change simluation box.")
        else:
            self._step.Box.xlo = new_box[0]
            self._step.Box.ylo = new_box[1]
            self._step.Box.zlo = new_box[2]
            self._step.Box.xhi = new_box[3]
            self._step.Box.yhi = new_box[4]
            self._step.Box.zhi = new_box[5]
            self._step.Box.xy = 0
            self._step.Box.xz = 0
            self._step.Box.yz = 0
            self._step.Box.triclinic = False
            if triclinic:
                self._step.Box.xy = new_box[6]
                self._step.Box.xz = new_box[7]
                self._step.Box.yz = new_box[8]
                self._step.Box.triclinic = True
            self._step.Box.set_box(False)

    def change_masses(self,new_masses):
        """
        change_masses(new_masses):
            Change masses of each atom type. Overwrite current masses by new_massses.

        Parameters:
            new_masses: python list (double) or numpy 1d array (double)
                Mass of each atom type. The mass of ith atom type is new_masses[i]. new_masses[0] = 0 is a redundant element (place holder).
                The number of elements (including the 1st redundant element new_masses[0]) should be the same as 1 + the number of atom types
                in the current snapshot.
        
        Return:
            None.
        """
        if len(new_masses) - 1 != self._step.NumOfAtomTypes:
            raise ValueError("The number of atom types in new_masses is not the same as that of the current snapshot.")
        _ConvertToC_Double(self._step.Masses, new_masses)

    def change_atoms(self,new_atoms,new_props):
        """
        change_atoms(new_atoms)
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
        elif 'type' not in new_props:
            print("Missing the necessary property \"type\".")
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Failed to change atoms data.")
        else:
            _ConvertToC_2dDouble(self._step.Atoms, new_atoms)
            self._step.NumOfAtoms = new_atoms.shape[0]
            TYPE = new_props.index("type")
            self._step.NumOfAtomTypes = int(np.max(new_atoms[:,TYPE]))
            masses = list(self._step.Masses)
            if len(masses) - 1 > self._step.NumOfAtomTypes:
                new_masses = list(self._step.Masses)
                new_masses = new_masses[:self._step.NumOfAtomTypes + 1]
                _ConvertToC_Double(self._step.Masses, new_masses)
                print("Due to the change in number of atom types, masses section has also been modified by removing redundant types.")
            elif len(masses) - 1 < self._step.NumOfAtomTypes:
                new_masses = masses + [0.0] * (self._step.NumOfAtomTypes - len(masses) + 1)
                _ConvertToC_Double(self._step.Masses, new_masses)
                print("Due to the change in number of atom types, masses section has also been modified by appending extra types with mass=0.")
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
            self._data.AtomIDMapping(self._idmapping_atom, self._step)
            
    def change_bonds(self,new_bonds):
        """
        change_bonds(new_bonds)
            Change bonds data. Overwrite current bonds data by a new numpy 2d array (new_bonds).
            If remove all bonds, then the input new_bonds should be np.array([]).reshape(0,0)

        Parameters:
            new_bonds: numpy 2d array (int)
                A numpy 2d array with bonds data.

        Return:
            None.
        """
        if not (isinstance(new_bonds, np.ndarray) and new_bonds.ndim == 2):
            raise RuntimeError("The input new_bonds should be an 2d numpy array.")
        if len(new_bonds) == 0:
            self._step.NumOfBonds = 0
            self._step.NumOfBondTypes = 0
            _ConvertToC_2dInt(self._step.Bonds, new_bonds.astype(np.int32))
        elif new_bonds.shape[1] != 4:
            print("Wrong size of 2d array with new bonds data.")
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Failed to change bonds data.")
        else:
            _ConvertToC_2dInt(self._step.Bonds, new_bonds.astype(np.int32))
            self._step.NumOfBonds = new_bonds.shape[0]
            self._step.NumOfBondTypes = int(np.max(new_bonds[:,1]))

    def change_angles(self,new_angles):
        """
        change_angles(new_angles)
            Change angles data. Overwrite current angles data by a new numpy 2d array (new_angles).
            If remove all angles, then the input new_angles should be np.array([]).reshape(0,0)

        Parameters:
            new_angles: numpy 2d array (int)
                A numpy 2d array with angles data.

        Return:
            None.
        """
        if not (isinstance(new_angles, np.ndarray) and new_angles.ndim == 2):
            raise RuntimeError("The input new_angles should be an 2d numpy array.")
        if len(new_angles) == 0:
            self._step.NumOfAngles = 0
            self._step.NumOfAngleTypes = 0
            _ConvertToC_2dInt(self._step.Angles, new_angles.astype(np.int32))
        elif new_angles.shape[1] != 5:
            print("Wrong size of 2d array with new angles data.")
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Failed to change angles data.")
        else:
            _ConvertToC_2dInt(self._step.Angles, new_angles.astype(np.int32))
            self._step.NumOfAngles = new_angles.shape[0]
            self._step.NumOfAngleTypes = int(np.max(new_angles[:,1]))

    def change_dihedrals(self,new_dihedrals):
        """
        change_dihedrals(new_dihedrals)
            Change dihedrals data. Overwrite current dihedrals data by a new numpy 2d array (new_dihedrals).
            If remove all dihedrals, then the input new_dihedrals should be np.array([]).reshape(0,0)

        Parameters:
            new_dihedrals: numpy 2d array (int)
                A numpy 2d array with dihedrals data.

        Return:
            None.
        """
        if not (isinstance(new_dihedrals, np.ndarray) and new_dihedrals.ndim == 2):
            raise RuntimeError("The input new_dihedrals should be an 2d numpy array.")
        if len(new_dihedrals) == 0:
            self._step.NumOfDihedrals = 0
            self._step.NumOfDihedralTypes = 0
            _ConvertToC_2dInt(self._step.Dihedrals, new_dihedrals.astype(np.int32))
        elif new_dihedrals.shape[1] != 6:
            print("Wrong size of 2d array with new dihedrals data.")
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Failed to change dihedrals data.")
        else:
            _ConvertToC_2dInt(self._step.Dihedrals, new_dihedrals.astype(np.int32))
            self._step.NumOfDihedrals = new_dihedrals.shape[0]
            self._step.NumOfDihedralTypes = int(np.max(new_dihedrals[:,1]))

    def change_impropers(self,new_impropers):
        """
        change_impropers(new_impropers)
            Change impropers data. Overwrite current impropers data by a new numpy 2d array (new_impropers).
            If remove all impropers, then the input new_impropers should be np.array([]).reshape(0,0)

        Parameters:
            new_impropers: numpy 2d array (int)
                A numpy 2d array with impropers data.

        Return:
            None.
        """
        if not (isinstance(new_impropers, np.ndarray) and new_impropers.ndim == 2):
            raise RuntimeError("The input new_impropers should be an 2d numpy array.")                    
        if len(new_impropers) == 0:
            self._step.NumOfImpropers = 0
            self._step.NumOfImproperTypes = 0
            _ConvertToC_2dInt(self._step.Impropers, new_impropers.astype(np.int32))
        elif new_impropers.shape[1] != 6:
            print("Wrong size of 2d array with new impropers data.")
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Failed to change impropers data.")
        else:
            _ConvertToC_2dInt(self._step.Impropers, new_impropers.astype(np.int32))
            self._step.NumOfImpropers = new_impropers.shape[0]
            self._step.NumOfImproperTypes = int(np.max(new_impropers[:,1]))

    def reset_atomid(self):
        """
        Reset the atom ID of atoms, so they are numbered contiguously.
        This will also change the bonds, angles, dihedrals and imporpers section accordingly.
        """
        atoms = _ConvertToNumpy_VecVecDouble(self._step.Atoms)
        bonds = _ConvertToNumpy_VecVecInt(self._step.Bonds)
        angles = _ConvertToNumpy_VecVecInt(self._step.Angles)
        dihedrals = _ConvertToNumpy_VecVecInt(self._step.Dihedrals)
        impropers = _ConvertToNumpy_VecVecInt(self._step.Impropers)
        _, ID = self._data.PropertyMapping("id", self._step)
        sortedid = np.sort(atoms[:,ID]).astype(int)
        old_to_new = {sortedid[i]: i + 1 for i in range(len(sortedid))}
        for i in range(len(atoms)):
            tmp = int(atoms[i,ID])
            atoms[i,ID] = old_to_new[tmp]
        for i in range(len(bonds)):
            for j in range(2):
                tmp = bonds[i,j+2]
                bonds[i,j+2] = old_to_new[tmp]
        for i in range(len(angles)):
            for j in range(3):
                tmp = angles[i,j+2]
                angles[i,j+2] = old_to_new[tmp]
        for i in range(len(dihedrals)):
            for j in range(4):
                tmp = dihedrals[i,j+2]
                dihedrals[i,j+2] = old_to_new[tmp]
        for i in range(len(impropers)):
            for j in range(4):
                tmp = impropers[i,j+2]
                impropers[i,j+2] = old_to_new[tmp]
        _ConvertToC_2dDouble(self._step.Atoms, atoms)
        self._data.AtomIDMapping(self._idmapping_atom, self._step)
        _ConvertToC_2dInt(self._step.Bonds, bonds.astype(np.int32))
        _ConvertToC_2dInt(self._step.Angles, angles.astype(np.int32))
        _ConvertToC_2dInt(self._step.Dihedrals, dihedrals.astype(np.int32))
        _ConvertToC_2dInt(self._step.Impropers, impropers.astype(np.int32))

    def flip(self):    
        """
        Flip the simulation box if the tilt factor is larger than 0.5
        """
        self._step.flip()


    def write_data(self,filename,atomstyle):
        """
        write(filename,atomstyle)
            Write the current step into a LAMMPS data file (filename) in a certain format (atomstyle).

        Parameters:
            filename: string
                Name of output LAMMPS data file.
            atomstyle: string
                Name of atom style ("atomic", "molecular", "sphere", or "full").

        Return:
            None.
        """
        ret = self._writer.Write(filename, self._step, atomstyle)
        if ret == 0:
            if self._print_flag:
                print("timestep = %d is written to file \"%s\"" % (self._step.Timestep, filename))
        else:
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Failed to write data into file.")

    def write_dump(self,filename,properties,mode="a"):
        """
        write_dump(filename, properties, mode):
            Write the current snapshot to a LAMMPS dump file `filename`, including properties listed in `properties`, with designated mode.

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
        self._step.Box.xboundary = "pp"
        self._step.Box.yboundary = "pp"
        self._step.Box.zboundary = "pp"
        ret = writer.Write(self._step)
        if ret == 0:
            if self._print_flag:
                print("timestep = %d is written to file \"%s\"" % (self._step.Timestep, filename))
        else:
            print("\x1b[31m" + "Warning: " + "\x1b[0m" + "Failed to write dump into file.")

    ############################
    ## backward compatibility ##
    ############################
    def changetimestep(self, new_timestep):
        self.change_timestep(new_timestep)

    def changebox(self, new_box, triclinic=False):
        self.change_box(new_box, triclinic=False)

    def changemasses(self,new_masses):
        self.change_masses(new_masses)

    def changeatoms(self,new_atoms,new_props):
        self.change_atoms(new_atoms,new_props)

    def changebonds(self,new_bonds):
        self.change_bonds(new_bonds)

    def changeangles(self,new_angles):
        self.change_angles(new_angles)

    def changedihedrals(self,new_dihedrals):
        self.change_dihedrals(new_dihedrals)

    def changeimpropers(self,new_impropers):
        self.change_impropers(new_impropers)
