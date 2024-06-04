class neigh:
    """
    Analyze neighbor enviroment for a simulation snapshot.

    Usage:
        ne = neigh()                                                # Declare an instance of neigh class.
        ne.set_snapshot(d)                                          # Set the simulation snapshot for analysis. "d" can be an instance of data class or dump class.
        ne.set_topology(d, inter_mol=True, intra_mol=True, exclude_bond=False, exclude_angle=False, exclude_dihedral=False, exclude_improper=False)
                                                                    # Set the simluation topology for analysis. "d" should be an instance of data class.
        neigh_list, dist_list = ne.query(cutoff, centerAtomIndex, force_prepare=True)
                                                                    # Find the neighbor atoms within a cutoff distance from a given atom.
        neigh_list, dist_list = ne.neighbor_analysis(cutoff, center_mask=None, neighbor_mask=None)
                                                                    # Build the neighbor list.
        cluster_id = ne.cluster_analysis(cutoff, cluster_mask=None)
                                                                    # Cluster analysis based on cutoff distance criterion
        rdf = ne.radial_distribution_function(cutoff, nbins, center_mask, neighbor_mask)
                                                                    # Analyze the (partial) radial distribution function
        rdf_dict = ne.radial_distribution_function_by_type(cutoff, nbins)
                                                                    # Analyze the partial radial distribution function for all atom types
    """
    def __init__(self):
        """
        neigh()

        Parameters:
            None.

        Return:
            None.
        """
        self._neigh = _Neigh()
        self._prepare_flag = False    ### whether cell list is prepared by c++ member function "Prepare" of class "Neigh"

        self._step = None
        self._snapshot_set_flag = False
        self._default_atoms_mask = None
        self._C_default_atoms_mask = _VecInt()
        
        self._step_topology = None
        self._topology_set_flag = False
        self._inter_mol = True
        self._intra_mol = True
        self._exclude_bond = False
        self._exclude_angle = False
        self._exclude_dihedral = False
        self._exclude_improper = False

    def set_snapshot(self, snapshot):
        """
        set_snapshot(snapshot):
            Set the simulation snapshot for analysis.

        Parameters:
            snapshot: a instance of class "data" or "dump"
                A instance of class "data" or "dump" which contains the simulation snapshot information for analysis.
                Read the snapshot for analysis and wrap the coordinate before calling this function.

        Return:
            None.

        Example:
            d = data()
            d.read("data.atom")
            d.wrap()
            ne = neigh()
            ne.set_snapshot(d)

            or

            d = dump("dump.atom")
            d.read(i)
            d.wrap()
            ne = neigh()
            ne.set_snapshot(d) # this set the ith frame in the dump file
        """
        self._step = _STEP(snapshot._step)
        if self._step.Timestep == -1:
            raise RuntimeError("Please read one snapshot first.")
        self._snapshot_set_flag = True
        self._default_atoms_mask = np.ones(self._step.NumOfAtoms).astype(np.int32)
        _ConvertToC_Int(self._C_default_atoms_mask, self._default_atoms_mask)
        self._prepare_flag = False

    def set_topology(self, snapshot, inter_mol=True, intra_mol=True, exclude_bond=False, exclude_angle=False, exclude_dihedral=False, exclude_improper=False):
        """
            inter_mol: bool, optional
                Whether to account for inter-molecular atom pairs.
                Default is True. 
            intra_mol: bool, optional
                Whether to account for intra-molecular atom pairs.
                Default is True. 
            exclude_bond: bool, optional
                Whether to exclude atom pairs if they interact with each other through bond potential (i.e. connected by bond).
                Default is False.
            exclude_angle: bool, optional
                Whether to exclude atom pairs if they interact with each other through angle potential.
                Default is False.
            exclude_dihedral: bool, optional
                Whether to exclude atom pairs if they interact with each other through dihedral potential.
                Default is False.
            exclude_improper: bool, optional
                Whether to exclude atom pairs if they interact with each other through improper potential.
                Default is False.

        NOTE:
            Whether to set the topology? 
            If user knows that the topology (bonds, angles, dihedrals and impropers) are not changed compared to
            the snapshot set previously, this can be set to False. This could happen when analyzing a series of snapshot from a 
            LAMMPS dump file while the topology is determined by a LAMMPS data file. However, when analyzing the first snapshot 
            of a series of snapshots (e.g, snapshots in dump file), user need to set this parameter to True to set the topology
            and make sure the snapshot used as the input of the member function "set_snapshot" is a instance of data class (so it
            contains topology information).

        """
        self._step_topology = _STEP(snapshot._step)
        ret = self._neigh.SetTopology(self._step_topology, inter_mol, intra_mol, exclude_bond, exclude_angle, exclude_dihedral, exclude_improper)
        if ret != 0:
            raise RuntimeError("Failed to set topology for rdf calculation.")
        self._inter_mol = inter_mol
        self._intra_mol = intra_mol
        self._exclude_bond = exclude_bond
        self._exclude_angle = exclude_angle
        self._exclude_dihedral = exclude_dihedral
        self._exclude_improper = exclude_improper
        self._topology_set_flag = True

    def query(self, cutoff, centerAtomIndex, force_prepare=True):
        """
        query(cutoff, centerAtomIndex, force_prepare=False):
            Find the neighbor atoms within a cutoff distance from a given atom.

        Parameters:
            cutoff: float
                Cutoff distance from the center atom within which other atoms are considered as neighbor atoms.
            centerAtomIndex: int
                The Atom Index of the atom of which the neighbor will be given.
            force_prepare: boolean, optional
                Force the neighbor analyzer to prepare cell list with the "cutoff" setting.
                By default (True), it will force the preparation. If the cutoff distance is not changed between
                two calls of this function, then this parameter can be set to False.

        Return:
            neighbors: python list (int)
                A list of Atom Index of atoms that locate within cutoff distance from the atom
                whose Atom Index is "centerAtomIndex". The list is not ordered.
            neighbors_distance: python list (double)
                Distance to the neighbor atoms. The ith element in the list corresponds to the distance to
                the ith atom in the list "neighbors".

        Example:
            d = dump("dump.atom")
            d.read(0)
            d.wrap()
            timestep,natoms,box,props,atoms = d.snapshot()
            ID = d.map("id")
            idmapping = d.idtoindex()
            ne = neigh()
            ne.set_snapshot(d) 
            cutoff = 2.5
            centerID = 1
            neighbors, distances = ne.query(cutoff, idmapping[centerID])
            # sorting neighbor atoms according to atom ID
            neighbors = atoms[neighbors,ID]
            sort_id = np.argsort(neighbors)
            neighbors = neighbors[sort_id]
            distances = distances[sort_id]
            print("Within the cutoff distance %g, IDs of neighbor atoms of the center atom ID=%d are:" % (cutoff, centerID))
            print(neighbors.astype(int)) 
        """
        if not self._snapshot_set_flag:
            raise RuntimeError("Please set the snapshot first. Use the member function set_snapshot.")
        if self._step.wrapped_flag != 1:
            raise RuntimeError("Please wrap the coordinate first.")

        if self._prepare_flag == False or force_prepare:
            ret = self._neigh.Prepare(self._step, cutoff, self._C_default_atoms_mask)
            if ret != 0:
                raise RuntimeError("Failed to prepare cell list.")
            self._prepare_flag = True

        tmp_vec_int = _VecInt()
        tmp_vec_double = _VecDouble()
        ret = self._neigh.Query(self._step, cutoff, self._C_default_atoms_mask, int(centerAtomIndex), tmp_vec_int, tmp_vec_double)
        if ret != 0:
            raise RuntimeError("Failed to query.")
        neighbors = _ConvertToNumpy_VecInt(tmp_vec_int)
        distances = _ConvertToNumpy_VecDouble(tmp_vec_double)
        distances = distances**0.5
        return list(neighbors), list(distances)

    def neighbor_analysis(self, cutoff, center_mask=None, neighbor_mask=None):
        """
        neighbor_analysis(cutoff, center_mask=None, neighbor_mask=None):
            Build the neighbor list for all center atoms. If user customizes center atoms and neighbor atoms,
            atoms which are selected by neighbor_mask will be listed if they are in the "cutoff" distance from
            the atom which is selected by center_mask.

        Parameters:
            cutoff: float
                Cutoff distance from the center atom within which other atoms are considered as neighbor atoms
            center_mask: numpy 1d array (int or boolean), optional
                A mask used to select center atoms. 
                If center_mask is a numpy int array, then 1 stands for True, and 0 stands for False. The length 
                of this array should be the same as the number of atoms in current snapshot.
                By default (None), all atoms in current snapshot can be center atoms.
            neighbor_mask: numpy 1d array (int or boolean), optional
                A mask used to select neighbor atoms. 
                If neighbor_mask is a numpy int array, then 1 stands for True, and 0 stands for False. The length 
                of this array should be the same as the number of atoms in current snapshot.
                By default (None), all atoms in current snapshot can be neighbor atoms.

        Return:
            neighbor_list: python list of tuples
                Atom index of neighbor atoms for each atom in current snapshot. i.e. neighbor_list[index] is 
                a list of atom index of neighbor atoms for a center atom whose atom index is "index".
                The number of rows (tuples) in this list is the number of atoms in current snapshot. If one atom has
                no neighbor atom, or it is not selected as center atom, the corresponding row (tuple) of this list
                will be empty.
            distance_list: python list of tuples
                The distance from the center atom to the corresponding neighbor atom in "neighbor_list".
                distance_list[i][j] is the distance between two atoms whose atom indices (not atom id) are i and 
                neighbor_list[i][j].

        Example:
            d = dump("dump.atom")
            d.read(0)
            d.wrap()
            timestep,natoms,box,props,atoms = d.snapshot()
            ID, TYPE = d.map("id", "type")
            idmapping = d.idtoindex()
            ne = neigh()
            ne.set_snapshot(d) 
            cutoff = 2.5
            center_mask = atoms[:,TYPE] == 1
            neighbor_mask = atoms[:,TYPE] == 2
            neighbor_list, distance_list = ne.neighbor_analysis(cutoff, center_mask, neighbor_mask)
            ### print out neighbor information
            atomindex = 0
            #atomid = 1
            #atomindex = idmapping[1]
            if len(neighbor_list[atomindex]) == 0:
                print("No neighbor for atom id = %d." % (atoms[atomindex,ID]))
            else:
                print("For atom with id = %d, atom id of neighbor atoms and distance to them are:" % atoms[atomindex,ID])
                for i in range(len(neighbor_list[atomindex])):
                    print("Neighbor atom id = %d, distance = %g" % (atoms[neighbor_list[atomindex][i],ID], distance_list[atomindex][i]))
        """
        if not self._snapshot_set_flag:
            raise RuntimeError("Please set the snapshot first. Use the member function set_snapshot.")
        if self._step.wrapped_flag != 1:
            raise RuntimeError("Please wrap the coordinate first.")

        if center_mask is None:
            center_atoms_mask = np.ones(self._step.NumOfAtoms).astype(np.int32)
        elif len(center_mask) != self._step.NumOfAtoms:
            raise ValueError("Size of center_mask should be the same as the number of atoms in current snapshot.")
        else:
            center_atoms_mask = center_mask.astype(np.int32)

        if neighbor_mask is None:
            neighbor_atoms_mask = np.ones(self._step.NumOfAtoms).astype(np.int32)
        elif len(neighbor_mask) != self._step.NumOfAtoms:
            raise ValueError("Size of neighbor_mask should be the same as the number of atoms in current snapshot.")
        else:
            neighbor_atoms_mask = neighbor_mask.astype(np.int32)

        Cmask_center = _VecInt()
        Cmask_neighbor = _VecInt()
        _ConvertToC_Int(Cmask_center, center_atoms_mask)
        _ConvertToC_Int(Cmask_neighbor, neighbor_atoms_mask)

        tmp_vec_vec_int = _VecVecInt()
        tmp_vec_vec_double = _VecVecDouble()
        ret = self._neigh.BuildNeighborList(self._step, cutoff, Cmask_center, Cmask_neighbor, tmp_vec_vec_int, tmp_vec_vec_double)
        if ret != 0:
            raise RuntimeError("Failed to build neighbor list.")
        _Sqrt_2dDouble(tmp_vec_vec_double)
        neighbor_list = list(tmp_vec_vec_int)
        distance_list = list(tmp_vec_vec_double)

        ### The c++ function "BuildNeighborList" calls the c++ function "Prepare" with "Cmask_neighbor"
        ### as a parameter. Since python function "query" uses c++ function "Query" to find the neighbors
        ### of a specific atom, the neighbor finder should be re-prepare with "self._C_default_atoms_mask"
        ### in order for the python function "query" to give the corret neighbor list. Refer to the 
        ### c++ header file Neigh.h for detailed reasons.
        self._prepare_flag = False

        return neighbor_list, distance_list

    def cluster_analysis(self, cutoff, cluster_mask=None):
        """
        cluster_analysis(cutoff, cluster_mask=None):
            Group atoms into clusters based on distance between atoms. Within one cluster, each atom
            is within the "cutoff" distance from at least one atom in such cluster.

        Parameters:
            cutoff: float
                Cutoff distance criterion for cluster analysis.
            cluster_mask: numpy 1d array (int or boolean), optional
                A mask used to select atoms for cluster analysis.
                If cluster_mask is a numpy int array, then 1 stands for True, and 0 stands for False. The length 
                of this array should be the same as the number of atoms in current snapshot.
                By default (None), cluster analysis is implemented for all atoms in current snapshot.

        Return:
            cluster_id: numpy 1d array (int)
                The cluster ID of each atom in current snapshot. The cluster ID is sorted in descending order,
                i.e. cluster ID=1 is the largest cluster. For atoms that cluster analysis is not applied to, 
                their cluster ID=0. cluster_id[i] is the cluster ID of the atom whose atom index is "i".

        Example:
            d = dump("dump.atom")
            d.read(0)
            d.wrap()
            timestep,natoms,box,props,atoms = d.snapshot()
            TYPE = d.map("type")
            idmapping = d.idtoindex()
            ne = neigh()
            ne.set_snapshot(d) 
            cutoff = 2.5
            cluster_mask = atoms[:,TYPE] == 1
            cluster_id = ne.cluster_analysis(cutoff, cluster_mask)
        """
        if not self._snapshot_set_flag:
            raise RuntimeError("Please set the snapshot first. Use the member function set_snapshot.")
        if self._step.wrapped_flag != 1:
            raise RuntimeError("Please wrap the coordinate first.")

        if cluster_mask is None:
            cluster_atoms_mask = np.ones(self._step.NumOfAtoms).astype(np.int32)
        elif len(cluster_mask) != self._step.NumOfAtoms:
            raise ValueError("Size of cluster_mask should be the same as the number of atoms in current snapshot.")
        else:
            cluster_atoms_mask = cluster_mask.astype(np.int32)

        Cmask = _VecInt()
        _ConvertToC_Int(Cmask, cluster_atoms_mask)

        tmp_vec_int = _VecInt()
        ret = self._neigh.Cluster(self._step, cutoff, Cmask, tmp_vec_int)
        if ret != 0:
            raise RuntimeError("Failed to analyze clusters.")
        cluster_id = _ConvertToNumpy_VecInt(tmp_vec_int)

        self._prepare_flag = False    ### see function "neighbor_analysis" for the reason that set the flag to False

        return cluster_id

    def radial_distribution_function(self, cutoff, nbins, center_mask=None, neighbor_mask=None):
        """
        radial_distribution_function(cutoff, nbins, center_mask=None, neighbor_mask=None) 
            Calculate the radial distribution function (rdf). User can also calculate the partial rdf
            by selecting center and neighbor atoms through setting "center_mask" and "neighbor_mask".

        Parameters:
            cutoff: double
                Cutoff distance for rdf analysis.
            nbins: int
                Number of bins in the rdf
            center_mask: numpy 1d array (int or boolean), optional
                A mask used to select center atoms.
                If center_mask is a numpy int array, then 1 stands for True, and 0 stands for False. The length
                of this array should be the same as the number of atoms in current snapshot.
                By default (None), all atoms in current snapshot can be center atoms.
            neighbor_mask: numpy 1d array (int or boolean), optional
                A mask used to select neighbor atoms.
                If neighbor_mask is a numpy int array, then 1 stands for True, and 0 stands for False. The length
                of this array should be the same as the number of atoms in current snapshot.
                By default (None), all atoms in current snapshot can be neighbor atoms.

        Return:
            rdf: numpy 2d array (double)
                Radial distribution function g(r). The first column is r, which chosen to be the center of the bin.
                The second column is g.

        Example:
            import matplotlib.pyplot as plt
            ne = neigh()
            d = data()
            d.read("sys.data")
            TYPE = d.map("type")
            timestep,natoms,box,props,atoms = d.snapshot()
            ne.set_snapshot(d)
            ### exclude 1-2, 1-3, 1-4 non bonded pairs
            ne.set_topology(d, inter_mol=True, intra_mol=True, exclude_bond=True, exclude_angle=True, exclude_dihedral=True, exclude_improper=False)
            ### set masks
            cmask = atoms[:,TYPE] == 1
            nmask = atoms[:,TYPE] == 2
            ### calcualte the rdf of atoms with atom type 2 around atoms with atom type 1
            rdf = ne.radial_distribution_function(10, 1000, center_mask=cmask, neighbor_mask=nmask)
            plt.plot(rdf[:,0], rdf[:,1])        ### plot g(r) vs r
        """
        if not self._snapshot_set_flag:
            raise RuntimeError("Please set the snapshot first. Use the member function set_snapshot.")
        if self._step.wrapped_flag != 1:
            raise RuntimeError("Please wrap the coordinate first.")
        if not self._topology_set_flag:
            raise RuntimeError("Please set topology first. Use the member function set_topology.")


        if center_mask is None:
            center_atoms_mask = np.ones(self._step.NumOfAtoms).astype(np.int32)
        elif len(center_mask) != self._step.NumOfAtoms:
            raise ValueError("Size of center_mask should be the same as the number of atoms in current snapshot.")
        else:
            center_atoms_mask = center_mask.astype(np.int32)

        if neighbor_mask is None:
            neighbor_atoms_mask = np.ones(self._step.NumOfAtoms).astype(np.int32)
        elif len(neighbor_mask) != self._step.NumOfAtoms:
            raise ValueError("Size of neighbor_mask should be the same as the number of atoms in current snapshot.")
        else:
            neighbor_atoms_mask = neighbor_mask.astype(np.int32)

        Cmask_center = _VecInt()
        Cmask_neighbor = _VecInt()
        _ConvertToC_Int(Cmask_center, center_atoms_mask)
        _ConvertToC_Int(Cmask_neighbor, neighbor_atoms_mask)

        tmp_vec_double_r = _VecDouble()
        tmp_vec_double_g = _VecDouble()


        ret = self._neigh.RDF(self._step, cutoff, nbins, Cmask_center, Cmask_neighbor, tmp_vec_double_r, tmp_vec_double_g)
        if ret != 0:
            raise RuntimeError("Failed to calculate rdf.")

        r = _ConvertToNumpy_VecDouble(tmp_vec_double_r)
        g = _ConvertToNumpy_VecDouble(tmp_vec_double_g)

        self._prepare_flag = False    ### see function "neighbor_analysis" for the reason that set the flag to False

        return np.column_stack((r, g))

    def radial_distribution_function_by_type(self, cutoff, nbins):
        """
        radial_distribution_function_by_type(cutoff, nbins):
            Calculate the partial radial distribution function (rdf) for all atom types.

        Parameters:
            cutoff: double
                Cutoff distance for rdf analysis.
            nbins: int
                Number of bins in the rdf

        Return:
            rdf_dict: python dictionary, key: tuple of 2 int, value: numpy 2d array (double)
                A dictionary of partial radial distribution function. The key of this dictionary is the pair of atom types. The value
                is a 2d numpy array, recording the partial rdf of such pair of atom types. The first column is r, which chosen to be 
                the center of the bin. The second column is g.
                NOTE: In the key (i,j), where i and j are atom types, i should not be larger than j, i.e. i <= j.
        
        Example:
            import matplotlib.pyplot as plt
            ne = neigh()
            d = data()
            d.read("sys.data")
            ne.set_snapshot(d)
            ### exclude 1-2, 1-3, 1-4 non bonded pairs
            ne.set_topology(d, inter_mol=True, intra_mol=True, exclude_bond=True, exclude_angle=True, exclude_dihedral=True, exclude_improper=False)
            rdf_dict = ne.radial_distribution_function_by_type(10, 1000)
            rdf_1_2 = rdf_dict[(1,2)]                   ### rdf of atoms with atom type 2 around atoms with atom type 1
            plt.plot(rdf_1_2[:,0], rdf_1_2[:,1])        ### plot g_1,2(r) vs r
        """
        if not self._snapshot_set_flag:
            raise RuntimeError("Please set the snapshot first. Use the member function set_snapshot.")
        if self._step.wrapped_flag != 1:
            raise RuntimeError("Please wrap the coordinate first.")
        if not self._topology_set_flag:
            raise RuntimeError("Please set topology first. Use the member function set_topology.")

        tmp_vec_vec_int = _VecVecInt()
        tmp_vec_vec_double_r = _VecVecDouble()
        tmp_vec_vec_double_g = _VecVecDouble()
        ret = self._neigh.RDFByTypes(self._step, cutoff, nbins, tmp_vec_vec_int, tmp_vec_vec_double_r, tmp_vec_vec_double_g)
        if ret != 0:
            raise RuntimeError("Failed to calculate partial rdf by types.")
        tmp_2dlist = list(tmp_vec_vec_int)
        pair_index = {(i, j): tmp_2dlist[i][j] for i in range(1,len(tmp_2dlist)) for j in range(i,len(tmp_2dlist[i]))}
        r = _ConvertToNumpy_VecVecDouble(tmp_vec_vec_double_r)
        g = _ConvertToNumpy_VecVecDouble(tmp_vec_vec_double_g)
        rdf_dict = {}
        for key in pair_index.keys():
            rdf_dict[key] = np.column_stack((r[pair_index[key]], g[pair_index[key]]))

        return rdf_dict
