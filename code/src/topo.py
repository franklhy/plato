class topo:
    """
    Analyze the topology of a simulation snapshot.

    Usage:
        tp = topo()                                    # Declare an instance of topo class.
        tp.set_snapshot(d)                             # Set the simulation snapshot for analysis. "d" should be an instance of data class.
        clusterid = tp.cluster_analysis()              # Analyze the connectivity of atoms and assign cluster id for each atoms
        connect_list = tp.connectivity(bond=True, angle=False, dihedral=False, improper=False)
                                                       # Build the connectivity list for each atom
    """

    def __init__(self):
        """
        topo()

        Parameters:
            None.

        Return:
            None.
        """
        self._topo = _Topo()
        self._step = None
        self._snapshot_set_flag = False

    def set_snapshot(self, snapshot):
        """
        set_snapshot(snapshot):
            Set the simulation snapshot for analysis.

        Parameters:
            snapshot: a instance of class "data" or "dump"
                A instance of class "data" which contains the topology information (bonds, angles, ...) for analysis.
                Read the snapshot for analysis before calling this function.

        Return:
            None.

        Example:
            d = data()
            d.read("data.atom")
            tp = topo()
            tp.set_snapshot(d)
        """
        self._step = _STEP(snapshot._step)
        if self._step.Timestep == -1:
            raise RuntimeError("Please read one snapshot first.")
        self._snapshot_set_flag = True

    def cluster_analysis(self):
        """
        cluster_analysis():
            Analyze the connectivity of atoms and assign cluster id for each atoms.
            Each cluster is a group of atoms connected together through bonds. Clusters are sorted according to the size in a descending order.

        Parameters:
            None.

        Return:
            cluster_id: numpy 1d array (int)
                The cluster ID of each atom in current snapshot. The cluster ID is sorted in descending order, i.e. cluster ID=1 is the largest cluster. 
                cluster_id[i] is the cluster ID of the atom whose atom index (not atom id) is "i".
        """
        if not self._snapshot_set_flag:
            raise RuntimeError("Please set the snapshot first.")

        tmp_vec_int = _VecInt()
        ret = self._topo.ClusterAnalysis(self._step, tmp_vec_int)
        if ret != 0:
            raise RuntimeError("Failed to analyze clusters.")

        cluster_id = _ConvertToNumpy_VecInt(tmp_vec_int)
        return cluster_id

    def connectivity(self, bond=True, angle=False, dihedral=False, improper=False):
        """
        connectivity(bond=True, angle=False, dihedral=False, improper=False):
            Build the connectivity list for each atom. User can chose the criteria for "connect": whether one atom is
            "connected" to other atoms by bond, angle, dihedral or/and improper interactions. For example, in the default
            case, if there exists a bond interaction between two atoms, these two atoms are considered "connected". If
            angle=True, then if there exists a angle interaction among three atoms, these three atoms are considered
            "connected" to each other.
            If two atoms, with atom index (not atom id) i and j respectively, are "connected", then i will appear in
            connect_list[j] and j will appear in connect_list[i]. Indices in connect_list[i] is not ordered.

        Parameters:
            bond: boolean
                Choose connected atoms based on bond interactions. Default is True.
            angle: boolean
                Choose connected atoms based on angle interactions. Default is False.
            dihedral: boolean
                Choose connected atoms based on dihedral interactions. Default is False.
            improper: boolean
                Choose connected atoms based on improper interactions. Default is False.

        Return:
            connect_list: python list of tuples
                Atom index of connected atoms for each atom in current snapshot. i.e. connect_list[index] is a tuple of
                atom index of atoms directly connected through bond to the center atom whose atom index is "index".
                Note: atom index is not atom id.
        """
        if not self._snapshot_set_flag:
            raise RuntimeError("Please set the snapshot first.")

        tmp_vec_vec_int = _VecVecInt()
        ret = self._topo.Connectivity(self._step, bond, angle, dihedral, improper, tmp_vec_vec_int)
        if ret != 0:
            raise RuntimeError("Failed to analyze clusters.")

        connect_list = list(tmp_vec_vec_int)
        
        return connect_list
