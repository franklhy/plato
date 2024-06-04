class hbond:
    '''
    Analyse the hydrogen bonds (D-H...A) from lammps data file.

                  H ... A (e.g., oxygen, nitrogen, ...)
                 /
                D (e.g., oxygen, nitrogen)

    Criteria: (1) Distance between Donor (D) and Acceptor (A) is smaller than max_DA_dist.
              (2) The angle ADH is smaller than max_ADH_angle. Or the angle DHA is larget than min_DHA_angle.
                  Which angle criterion to be used depends on user's choice

    The output format of hydrogen bond existence mapping is inspired by GROMACS (gmx hbond -hbn -hbm).
    '''
    #from plato import neigh  # don't need to import neigh, since class neigh will be included in the same python package (i.e. plato.py)
    from tqdm import tqdm
    def __init__(self, snapshot=None, D_type=[], H_type=[], A_type=[], max_DA_dist=3.0, max_ADH_angle=30, min_DHA_angle=150, angle_choice="ADH"):
        '''
        Parameters:
            snapshot: a instance of class "data"
                A instance of class "data" which contains the simulation snapshot information for analysis.
            D_type: python list
                A list of atom type which should be considered as donor (D).
            H_type: python list
                A list of atom type which should be considered as hydrogen (H).
            A_type: python list
                A list of atom type which should be considered as acceptor (A).
            max_DA_dist: float
                The maximum distance between a pair of donor (D) and acceptor (A).
            max_ADH_angle: float
                The maximum angle of acceptor(A)-donor(D)-hydrogen(H). In degree.
            min_DHA_angle: float
                The maximum angle of acceptor(A)-donor(D)-hydrogen(H). In degree.
            angle_choice: "ADH", "DHA", or "none"
                The choice of angle criteria.
                If "ADH", then the ADH angle should be smaller than max_ADH_angle.
                If "DHA", then the DHA angle should be larger than min_DHA_angle.
                If "none" (not the python object None, but a string), then angle criteria is not applied,
                only distance between the donor (D) and acceptor (A) is considered.
        '''
        self.snapshot = snapshot
        self.D_type = D_type
        self.H_type = H_type
        self.A_type = A_type
        self.max_DA_dist = max_DA_dist
        self.max_ADH_angle = max_ADH_angle
        self.min_DHA_angle = min_DHA_angle
        self.angle_choice = angle_choice
        if angle_choice not in ["ADH", "DHA", "none"]:
            raise RuntimeError("Invalide angle criteria choice. Please choose from [\"ADH\", \"DHA\", \"none\"].")

        self.HBond_DHA = {}
        self.timesteps = []
        self.HBond_count = []
        self.HBond_DHA_list = None
        self.HBond_existence_map = None


    def load_hbond(self, name=""):
        '''
        Load hydrogen bond existence mapping and index results analyzed before.

        Parameters:
            name: string
                Prefix of results file name.
        '''
        self.HBond_existence_map = np.loadtxt(name + "hbond_map.txt")

        counts = np.loadtxt(name + "hbond_num.txt")
        self.timesteps = counts[:,0].tolist()
        self.HBond_count = counts[:,1].astype(int).tolist()

        self.HBond_DHA_list = []
        hbindex = np.loadtxt(name + "hbond_index.txt")
        for hb in hbindex:
            self.HBond_DHA_list.append(tuple(hb.astype(int).tolist()))


    def run(self, start=None, stop=None, every=None, write_datafile=False):
        '''
        Run the analysis to find all hydrogen bonds

        Parameters:
            start: int
                Starting frame
            stop: int
                Stopping frame (included)
            every: int
                Analyze every this much frame
            write_datafile: bool
                If True, a LAMMPS data file will be written, with a timestep tag in the file name, which
                includes hydrogen bonds in "Bonds" sections, and the bond type of hydrogen bonds will be
                the maximum.
        '''
        if self.snapshot is None:
            return None

        self.HBond_DHA = {}
        self.timesteps = []
        self.HBond_count = []

        alltime = self.snapshot.traj_time()
        if start is None:
            start = 0
        if stop is None:
            stop = len(alltime)
        if every is None:
            every = 1

        frames = np.arange(len(alltime))
        for i in frames[start:stop:every]:
            timestep, HBond_DHA_i = self._analyze(i, write_datafile)
            self.timesteps.append(timestep)
            self.HBond_count.append(len(HBond_DHA_i))

            ### update hydrogen bond dictionary self.HBond_DHA
            for hb in HBond_DHA_i:
                key = tuple(hb.tolist())
                if key in self.HBond_DHA.keys():
                    self.HBond_DHA[key].append(int(i))
                else:
                    self.HBond_DHA[key]= [int(i)]


    def generate_hbond_index_and_map(self, name=None):
        '''
        Similar to GROMACS "gmx hbond -hbn -hbm", generate indexes of hydrogen bond and their existence map.
        The index will be saved in "XXX_hbond_index.txt" as a (#H-bonds,3) array, each line are atom id of D, H and A
        The existence map will be saved in "XXX_hbond_map.txt" as a (#H-bonds, #frames) array, each line is the condition of the hydrogen bond triplet
        of the same line in the index file ("XXX_hbond_index.txt"), 0/1 means the hydrogen bond is broken/formed.
        The timesteps will be saved in "XXX_hbond_num.txt" as a (2, #frames) array, each line is (timestep, #H-bonds at this timestep).

        Parameter:
            name: string
                File name prefix.
        '''
        if name is None:
            name = ""
        else:
            name = name + "_"
        np.savetxt(name + "hbond_num.txt", np.column_stack((self.timesteps, self.HBond_count)))

        self.HBond_DHA_list = sorted(self.HBond_DHA.keys(), key=lambda x: (x[0], x[1], x[2]))
        np.savetxt(name + "hbond_index.txt", np.array(self.HBond_DHA_list), fmt="%d\t%d\t%d", header="D\tH\tA")

        self.HBond_existence_map = np.zeros((len(self.HBond_DHA_list), len(self.timesteps)))
        for i, hb in enumerate(self.HBond_DHA_list):
            for t in self.HBond_DHA[hb]:
                self.HBond_existence_map[i][t] = 1
        np.savetxt(name + "hbond_map.txt", self.HBond_existence_map, fmt="%d")


    def autocorrelation_continuous(self, windows=None):
        '''
        Calculate the continuous autocorrelation function of hydrogen bonds.

        Parameters:
            windows: int
                Determines the maximum time (tau_max) between the start (t0) and the end (t0 + tau_max).
                tau_max = windows * dt, where dt is the time difference between two consecutive snapshots.
        '''
        time_span = self.HBond_existence_map.shape[1]
        if windows is None:
            windows = time_span // 2
        if windows >= time_span:
            raise RuntimeError("Windows should be smaller.")

        acf = np.zeros(windows)
        tau = np.arange(windows) * (self.timesteps[1] - self.timesteps[0])
        for taui in self.tqdm(range(windows)):
            for t0 in range(time_span - windows):
                acf[taui] += np.sum(np.all(self.HBond_existence_map[:,t0:t0+taui+1], axis=1)) / np.sum(self.HBond_existence_map[:,t0])
            acf[taui] /= time_span - windows

        return np.column_stack((tau, acf))


    def autocorrelation_intermittent(self, windows=None, intermittent=0):
        '''
        Calculate the intermittent autocorrelation function of hydrogen bonds.

        Parameters:
            windows: int
                The maximum number of frames between the start (t0) and the end (t0 + tau_max).
                tau_max = windows * dt, where dt is the time difference between two consecutive snapshots.
            intermittent: int
                The maximum number of frames for which a hydrogen bond is allowed to break.
                intermittent = 0 means on
        '''
        time_span = self.HBond_existence_map.shape[1]
        if windows is None:
            windows = time_span // 2
        if windows >= time_span:
            raise RuntimeError("Windows should be smaller.")

        hbmap_intermittent = self._intermittent_correction(intermittent)
        acf = np.zeros(windows)
        tau = np.arange(windows) * (self.timesteps[1] - self.timesteps[0])
        for taui in self.tqdm(range(windows)):
            for t0 in range(time_span - windows):
                acf[taui] += np.sum(hbmap_intermittent[:,t0] * hbmap_intermittent[:,t0+taui]) / np.sum(hbmap_intermittent[:,t0])
            acf[taui] /= time_span - windows

        return np.column_stack((tau, acf))



    def _analyze(self, frame, write_datafile):
        frame = int(frame)
        self.snapshot.load_traj(frame)
        self.snapshot.wrap()
        timestep,natoms,box,props,atoms = self.snapshot.snapshot()

        idmapping_atom = self.snapshot.idtoindex_atom()
        ID, TYPE, X, Y, Z = self.snapshot.map("id", "type", "x", "y", "z")
        lx = box[3] - box[0]
        ly = box[4] - box[1]
        lz = box[5] - box[2]

        ### find neighbors D-A pairs
        ne = neigh()
        ne.set_snapshot(self.snapshot)
        #ne.set_topology(snapshot, exclude_angle=True, exclude_dihedral=True)
        ne.set_topology(self.snapshot, exclude_bond=True)
        maskD = np.in1d(atoms[:,TYPE], self.D_type)
        maskA = np.in1d(atoms[:,TYPE], self.A_type)
        neigh_list, _ = ne.neighbor_analysis(self.max_DA_dist, center_mask=maskD, neighbor_mask=maskA)

        if len(self.H_type) != 0:
            ### Find all D with D-H bonds
            bonds,angles,dihedrals,impropers = self.snapshot.topology()
            maskDH1 = np.in1d(atoms[idmapping_atom[bonds[:,2]],TYPE], self.D_type)
            maskDH2 = np.in1d(atoms[idmapping_atom[bonds[:,3]],TYPE], self.H_type)
            maskDH = maskDH1 * maskDH2
            D_indexDH = idmapping_atom[bonds[maskDH][:,2]]
            H_indexDH = idmapping_atom[bonds[maskDH][:,3]]
            DH = np.column_stack((D_indexDH, H_indexDH))
            maskHD1 = np.in1d(atoms[idmapping_atom[bonds[:,2]],TYPE], self.H_type)
            maskHD2 = np.in1d(atoms[idmapping_atom[bonds[:,3]],TYPE], self.D_type)
            maskHD = maskHD1 * maskHD2
            D_indexHD = idmapping_atom[bonds[maskHD][:,3]]
            H_indexHD = idmapping_atom[bonds[maskHD][:,2]]
            HD = np.column_stack((D_indexHD, H_indexHD))
            ### All DH bonds and all D atom index
            DH = np.concatenate((DH, HD))
            D_index = np.unique(DH[:,0])

            DHA = []
            for D in D_index:
                for A in neigh_list[D]:
                    for H in DH[DH[:,0] == D][:,1]:
                        DHA.append([D, H, A])
            DHA = np.array(DHA)   # all possible DHA triplets, as candidates for hydrogen bonds

            if self.angle_choice != "none":
            ### when angle criteria is applied
                ### Vector DA
                DA = atoms[DHA[:,2]][:,(X,Y,Z)] - atoms[DHA[:,0]][:,(X,Y,Z)]
                ### Wrap vector DA
                DA[:,0] -= np.heaviside(np.absolute(DA[:,0]) - 0.5*lx, 0) * np.sign(DA[:,0]) * lx
                DA[:,1] -= np.heaviside(np.absolute(DA[:,1]) - 0.5*lx, 0) * np.sign(DA[:,1]) * ly
                DA[:,2] -= np.heaviside(np.absolute(DA[:,2]) - 0.5*lx, 0) * np.sign(DA[:,2]) * lz
                ### Vector DH
                DH = atoms[DHA[:,1]][:,(X,Y,Z)] - atoms[DHA[:,0]][:,(X,Y,Z)]
                ### Wrap vector DA
                DH[:,0] -= np.heaviside(np.absolute(DH[:,0]) - 0.5*lx, 0) * np.sign(DH[:,0]) * lx
                DH[:,1] -= np.heaviside(np.absolute(DH[:,1]) - 0.5*lx, 0) * np.sign(DH[:,1]) * ly
                DH[:,2] -= np.heaviside(np.absolute(DH[:,2]) - 0.5*lx, 0) * np.sign(DH[:,2]) * lz
                ### Vector HA
                HA = DA - DH
                ### DA distance
                DA_dist = np.sum(DA**2, axis=1)**0.5
                ### DH distance
                DH_dist = np.sum(DH**2, axis=1)**0.5
                ### HA distance
                HA_dist = np.sum(HA**2, axis=1)**0.5
                if self.angle_choice == "ADH":
                    ### Cosine of angle ADH
                    cosADH = np.sum(DA*DH, axis=1) / DA_dist / DH_dist
                    mask = cosADH > np.cos(self.max_ADH_angle/180*np.pi)
                elif self.angle_choice == "DHA":
                    ### Cosine of angle DHA
                    cosDHA = np.sum(-DH*HA, axis=1) / DH_dist / HA_dist
                    mask = cosDHA < np.cos(self.min_DHA_angle/180*np.pi)
                else:
                    raise RuntimeError("Invalid angle criteria choice.")
                ### Hydrogen bonds
                ### don't need to check DA distance since acceptor atoms (A) are choosen
                ### based on cutoff=max_DA_dist distance when build the neigh_list
                HBond_DHA = DHA[mask]
            else:
            ### when no angle criteria is applied
                HBond_DHA = DHA

        else:
        ### When no H is selected. This usually happens in two cases:
        ###   (1) finding hydrogen bond in coarse-grained simulation, where H is not explicited modeled.
        ###   (2) finding pairs of beads in close contact. In this case, donors (D) and acceptors (A) are not limited to hydrogen bonds,
        ###       but they can also be cations and anions of ion-pairs.
        ### In both cases, the criteria is only distance between donors (D) and acceptors (A)
            DHA = []
            D_index = np.where(maskD)[0]
            for D in D_index:
                for A in neigh_list[D]:
                    DHA.append([D, D, A]) ### since no H is selected, use D as a place holder
            HBond_DHA = np.array(DHA)

        indextoid = lambda x: atoms[x,ID]
        vfunc = np.vectorize(indextoid)
        HBond_DHA = vfunc(HBond_DHA)

        # Generate hydrogen bonds and save in LAMMPS data file
        if write_datafile:
            nHBond = len(HBond_DHA)
            A_id = HBond_DHA[:,0]
            H_id = HBond_DHA[:,2]
            hbonds = np.column_stack((np.arange(len(bonds)+1,len(bonds)+1+nHBond), np.ones(nHBond)*np.max(bonds[:,1]+1), A_id, H_id))
            hbonds = hbonds.astype(np.int32)
            new_bonds = np.concatenate((bonds, hbonds))
            self.snapshot.changebonds(new_bonds)
            self.snapshot.write_data("with_HBonds_%d.data" % timestep, "full")
            self.snapshot.changebonds(bonds)

        return timestep, HBond_DHA


    def _intermittent_correction(self, intermittent=0):
        '''
        For existence map of each hydrogen bond, e.g., [0, 0, 1, 1, 0, 0, 1, 0, 1, ...], if the number of consecutive 0 is less than (intermittent+1),
        then change these 0s to 1s.

        For example, if the existence map is [0,1,1,0,1,0,0,0,1,0,1] and intermittent = 1, then the existence map will be changed to [0,1,1,1,1,0,0,0,1,1,1]

        Parameters:
            intermittent: int
                The maximum number of frame a bond is allowed to break before reform. During the break, the hydrogen bond is still considered as not broken.
        '''
        hbmap_intermittent = np.copy(self.HBond_existence_map)

        delta = intermittent + 1
        if delta == 1:
            return hbmap_intermittent

        for i in range(len(hbmap_intermittent)):
            time_seq = np.where(hbmap_intermittent[i])[0]
            diff = time_seq[1:] - time_seq[:-1]
            mask = (diff <= delta) * (diff > 1)
            start = time_seq[:-1][mask]
            conut = diff[mask]
            for si in range(len(start)):
                hbmap_intermittent[i,start[si]+1:start[si]+count[si]] = 1

        return hbmap_intermittent


'''
datafile = "../equil.dat"
dumpfile = "../dump.atom"
H_type = [7, 15]
A_type = [5, 10]
D_type = [6, 12]

snapshot = data()
snapshot.read(datafile)
snapshot.set_traj_file(dumpfile)

HBA = HydrogenBond(snapshot, D_type=D_type, H_type=H_type, A_type=A_type, max_DA_dist=2.97, max_ADH_angle=50)
HBA.run()
HBA.generate_hbond_index_and_map()
'''

