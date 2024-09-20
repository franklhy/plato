/***********************************************************/
/*** Refer to header file for instructions of this class ***/
/***********************************************************/
/*** Version: 07/22/2022 (1.6 Version)  ***/
/*** By: Heyi Liang                     ***/
/******************************************/

#include "Neigh.h"

Neigh::Neigh()
{
    cutoff_ = 0.0;
    large_cutoff_ = false;
    numbinx_ = numbiny_ = numbinz_ = NumOfBins_ = 0;
    binsizex_ = binsizey_ = binsizez_ = 0.0;
    vector< vector<int> >().swap(bins_atom_index_);
    vector< vector<int> >().swap(bins_neigh_bins_);
}

Neigh::~Neigh()
{

}

void Neigh::print_version()
{
    printf("Running Neigh.cpp, version %s\n", NEIGH_VERSION);
}

bool Neigh::FindIn(const int target, const vector<int> &array)
{
    for(int i = 0; i < (int)array.size(); i++)
        if ( target == array[i] )
            return true;
    return false;
}

bool Neigh::ExcludeNeighByTopology(const STEP &cur_step, int center_index, int neigh_index) {
    int cid = cur_step.Atoms[center_index][cur_step.ID];
    int nid = cur_step.Atoms[neigh_index][cur_step.ID];
    bool exclude = FindIn(nid, connect_list_[cid]);
    bool same_mol = cur_step.Atoms[center_index][cur_step.MOL] == cur_step.Atoms[neigh_index][cur_step.MOL];

    if (!exclude && ((same_mol && include_intra_molecule_) || (!same_mol && include_inter_molecule_)) ) {
        return false;
    } else {
        /*** TEST ***/
        /***
        if (exclude)
            {printf("For atom id %d, neighbor atom id %d is excluded.\n", (int)cur_step.Atoms[center_index][cur_step.ID], (int)cur_step.Atoms[neigh_index][cur_step.ID]);}
        if (same_mol && !include_intra_molecule_)
            {printf("For atom id %d and neighbor atom id %d are intramol.\n", (int)cur_step.Atoms[center_index][cur_step.ID], (int)cur_step.Atoms[neigh_index][cur_step.ID]);}
        if (!same_mol && !include_inter_molecule_)
            {printf("For atom id %d and neighbor atom id %d are intermol.\n", (int)cur_step.Atoms[center_index][cur_step.ID], (int)cur_step.Atoms[neigh_index][cur_step.ID]);}
        ***/
        return true;
    }
}

int Neigh::SetTopology(const STEP &cur_step, const bool include_inter_molecule, const bool include_intra_molecule, \
const bool exclude_bond, const bool exclude_angle, const bool exclude_dihedral, const bool exclude_improper)
{
    // check
    if ( exclude_bond && cur_step.NumOfBonds == 0 )
    {
        fprintf(stderr, "Excluding bond is required, but no bond is recorded in current step.\n");
        return 1;
    }
    if ( exclude_angle && cur_step.NumOfAngles == 0 )
    {
        fprintf(stderr, "Excluding angle is required, but no angle is recorded in current step.\n");
        return 1;
    }
    if ( exclude_dihedral && cur_step.NumOfDihedrals == 0 )
    {
        fprintf(stderr, "Excluding dihedral is required, but no dihedral is recorded in current step.\n");
        return 1;
    }
    if ( exclude_improper && cur_step.NumOfImpropers == 0 )
    {
        fprintf(stderr, "Excluding improper is required, but no improper is recorded in current step.\n");
        return 1;
    }

    include_inter_molecule_ = include_inter_molecule;
    include_intra_molecule_ = include_intra_molecule;

    // build connectivity list (using atom ID, not atom index)
    Topo topo;
    vector< vector<int> > tmp_connect_list;
    if ( topo.Connectivity(cur_step, exclude_bond, exclude_angle, exclude_dihedral, exclude_improper, tmp_connect_list) )
        return 1;
    int max_id = 0;
    for (int i = 0; i < int(cur_step.Atoms.size()); i++)
        if (cur_step.Atoms[i][cur_step.ID] > max_id)
            max_id = cur_step.Atoms[i][cur_step.ID];
    int cid, nid;
    vector< vector<int> >().swap(connect_list_);
    connect_list_.resize(max_id + 1);
    for (int i = 0; i < int(tmp_connect_list.size()); i++)
    {
        cid = cur_step.Atoms[i][cur_step.ID];
        connect_list_[cid].resize(tmp_connect_list[i].size(), 0);
        for (int j = 0; j < int(tmp_connect_list[i].size()); j++)
        {
            nid = cur_step.Atoms[tmp_connect_list[i][j]][cur_step.ID];
            connect_list_[cid][j] = nid;
        }
    }
    return 0;
}

int Neigh::Prepare(const STEP &cur_step, const double p_cutoff, const vector<int> &neighbor_atoms_flag)    // single processor is enough
{
    if ( cur_step.Atoms.size() != neighbor_atoms_flag.size() )
    {
        fprintf(stderr, "Input vector<int> as the mask of neighbor atoms should has the same number of element as the number of atoms in current step\n");
        return 1;
    }


    if ( cur_step.wrapped_flag != 1 || cur_step.scaled_flag != 0 )
    {
        fprintf(stderr, "Please calculated the wrapped and unscaled coordinate first!\n");
        return 1;
    }

    /*** prepare bin neighbor list, and assign atoms to bins ***/
    cutoff_ = p_cutoff;
    // number of bins on each direction
    numbinx_ = (int)( (cur_step.Box.xhi - cur_step.Box.xlo) / cutoff_ );
    numbiny_ = (int)( (cur_step.Box.yhi - cur_step.Box.ylo) / cutoff_ );
    numbinz_ = (int)( (cur_step.Box.zhi - cur_step.Box.zlo) / cutoff_ );
    if ( numbinx_ < 2 || numbiny_ < 2 || numbinz_ < 2 )
    {
        printf("Cutoff cannot be larger than 1/2 of box length.\n");
        return 1;
    }
    if ( numbinx_ == 2 || numbiny_ == 2 || numbinz_ == 2 )
    {
        large_cutoff_ = true;
        return 0;
    }
    // length of bin edge on each direction
    binsizex_ = (cur_step.Box.xhi - cur_step.Box.xlo) / numbinx_;
    binsizey_ = (cur_step.Box.yhi - cur_step.Box.ylo) / numbiny_;
    binsizez_ = (cur_step.Box.zhi - cur_step.Box.zlo) / numbinz_;
    // number of bins
    NumOfBins_ = numbinx_ * numbiny_ * numbinz_;

    // reserve memory for vectors
    vector<int> tmp_v_int;
    vector< vector<int> >().swap(bins_atom_index_);
    bins_atom_index_.resize(NumOfBins_, tmp_v_int);
    for (int i = 0; i < NumOfBins_; i++)
        bins_atom_index_[i].reserve( (int)( (binsizex_ + 2) * (binsizey_ + 2) * (binsizez_ + 2) ) );

    tmp_v_int.resize(27,-1);
    vector< vector<int> >().swap(bins_neigh_bins_);
    bins_neigh_bins_.resize(NumOfBins_, tmp_v_int);

    // build bin neighbor lists for bins
    int cent_index, pair_index, pair_i, pair_j, pair_k;
    for (int cent_k = 0; cent_k < numbinz_; cent_k++)
        for (int cent_j = 0; cent_j < numbiny_; cent_j++)
            for (int cent_i = 0; cent_i < numbinx_; cent_i++)
            {
                cent_index = cent_i + cent_j * numbinx_ + cent_k * numbinx_ * numbiny_;
                // iterate over 27 neighbor bins (including ghost bins and the center bin itself)
                for (int delta_k = -1; delta_k <= 1; delta_k++)
                    for (int delta_j = -1; delta_j <= 1; delta_j++)
                        for (int delta_i = -1; delta_i <= 1; delta_i++)
                        {
                            pair_i = cent_i + delta_i;
                            pair_j = cent_j + delta_j;
                            pair_k = cent_k + delta_k;
                            // wrap ghost bins into simulation box
                            if ( pair_i >= numbinx_ )
                                pair_i = pair_i - numbinx_;
                            else if ( pair_i < 0 )
                                pair_i = pair_i + numbinx_;
                            if ( pair_j >= numbiny_ )
                                pair_j = pair_j - numbiny_;
                            else if ( pair_j < 0 )
                                pair_j = pair_j + numbiny_;
                            if ( pair_k >= numbinz_ )
                                pair_k = pair_k - numbinz_;
                            else if ( pair_k < 0 )
                                pair_k = pair_k + numbinz_;
                            // calculate the bin index of neighbor bin
                            pair_index = pair_i + pair_j * numbinx_ + pair_k * numbinx_ * numbiny_;
                            // record the neighbor (including the bin it self)
                            bins_neigh_bins_[cent_index][(delta_i+1) + (delta_j+1) * 3 + (delta_k+1) * 9] = pair_index; 
                        }
            }

    // assign atoms to bins
    if ( cur_step.X == -1 || cur_step.Y == -1 || cur_step.Z == -1 )
    {
        fprintf(stderr, "Cannot find complete (X,Y,Z) properties!\n");
        return 1;
    }
    int X,Y,Z;
    X = cur_step.X;
    Y = cur_step.Y;
    Z = cur_step.Z;

    int i, j, k;
    for (int ai = 0; ai < (int)cur_step.Atoms.size(); ai++)
    {
        if ( neighbor_atoms_flag[ai] == 1 )
        {
            i = (int)( ( cur_step.Atoms[ai][X] - cur_step.Box.xlo ) / binsizex_ );
            j = (int)( ( cur_step.Atoms[ai][Y] - cur_step.Box.ylo ) / binsizey_ );
            k = (int)( ( cur_step.Atoms[ai][Z] - cur_step.Box.zlo ) / binsizez_ );
            bins_atom_index_[i + j * numbinx_ + k * numbinx_ * numbiny_].push_back(ai);
        }
    }

    return 0;
}

int Neigh::Query(const STEP &cur_step, const double q_cutoff, const vector<int> &neighbor_atoms_flag, const int &center_AtomIndex, vector<int> &Vector_destin_index, vector<double> &Vector_destin_r2)
{    
    if ( cur_step.Atoms.size() != neighbor_atoms_flag.size() )
    {
        fprintf(stderr, "Input vector<int> as the mask of neighbor atoms should has the same number of element as the number of atoms in current step\n");
        return 1;
    }

    if ( q_cutoff > cutoff_ )
    {
        fprintf(stderr, "Query cutoff should not be larger than Prepare cutoff.\n");
        return 1;
    }
    if ( cur_step.X == -1 || cur_step.Y == -1 || cur_step.Z == -1 )
    {
        fprintf(stderr, "Cannot find complete (X,Y,Z) properties!\n");
        return 1;
    }

    double boxsize_x, boxsize_y, boxsize_z;
    boxsize_x = cur_step.Box.xhi - cur_step.Box.xlo;
    boxsize_y = cur_step.Box.yhi - cur_step.Box.ylo;
    boxsize_z = cur_step.Box.zhi - cur_step.Box.zlo;
    
    int X,Y,Z;
    X = cur_step.X;
    Y = cur_step.Y;
    Z = cur_step.Z;

    vector<int>().swap(Vector_destin_index);
    Vector_destin_index.reserve( (int)( 4 * pow(q_cutoff,3) ) );

    vector<double>().swap(Vector_destin_r2);
    Vector_destin_r2.reserve( (int)( 4 * pow(q_cutoff,3) ) );

    double dx, dy, dz, d_square, q_cutoff_square;
    q_cutoff_square = q_cutoff * q_cutoff;        // Calculate and compare d^2(distance square) instead of d. No sqrt calculation is required, so it save time.
    
    if ( large_cutoff_ )
    {
        for (int i = 0; i < (int)cur_step.Atoms.size(); i++)
        {
            if ( neighbor_atoms_flag[i] == 1 && !ExcludeNeighByTopology(cur_step, center_AtomIndex, i) )
            {
                dx = cur_step.Atoms[i][X] - cur_step.Atoms[center_AtomIndex][X];
                dy = cur_step.Atoms[i][Y] - cur_step.Atoms[center_AtomIndex][Y];
                dz = cur_step.Atoms[i][Z] - cur_step.Atoms[center_AtomIndex][Z];
                if ( dx > boxsize_x * 0.5 )
                    dx = dx - boxsize_x;
                else if ( dx < - boxsize_x * 0.5 )
                    dx = dx + boxsize_x;
                if ( dy > boxsize_y * 0.5 )
                    dy = dy - boxsize_y;
                else if ( dy < - boxsize_y * 0.5 )
                    dy = dy + boxsize_y;
                if ( dz > boxsize_z * 0.5 )
                    dz = dz - boxsize_z;
                else if ( dz < - boxsize_z * 0.5 )
                    dz = dz + boxsize_z;
                // calculate distance square
                d_square = dx*dx + dy*dy + dz*dz;
            
                if ( d_square <= q_cutoff_square && i != center_AtomIndex )
                {
                    Vector_destin_index.push_back(i);
                    Vector_destin_r2.push_back(d_square);
                }
            }
        }
    }
    else
    {
        // determine which bin is the center atom located at
        int i, j, k, center_BinIndex;
        i = (int)( ( cur_step.Atoms[center_AtomIndex][X] - cur_step.Box.xlo ) / binsizex_ );
        j = (int)( ( cur_step.Atoms[center_AtomIndex][Y] - cur_step.Box.ylo ) / binsizey_ );
        k = (int)( ( cur_step.Atoms[center_AtomIndex][Z] - cur_step.Box.zlo ) / binsizez_ );
        center_BinIndex = i + j * numbinx_ + k * numbinx_ * numbiny_;

        int tmp_AtomIndex;
        for (int i = 0; i < 27; i++)
        {
            if ( i != 13)    // for neighbor bins, we may have to wrap ghost atoms
            {
                for (int j = 0; j < (int)bins_atom_index_[ bins_neigh_bins_[center_BinIndex][i] ].size(); j++)
                {
                    tmp_AtomIndex = bins_atom_index_[ bins_neigh_bins_[center_BinIndex][i] ][j];
                    if ( !ExcludeNeighByTopology(cur_step, center_AtomIndex, tmp_AtomIndex) )
                    {
                        dx = cur_step.Atoms[tmp_AtomIndex][X] - cur_step.Atoms[center_AtomIndex][X];
                        dy = cur_step.Atoms[tmp_AtomIndex][Y] - cur_step.Atoms[center_AtomIndex][Y];
                        dz = cur_step.Atoms[tmp_AtomIndex][Z] - cur_step.Atoms[center_AtomIndex][Z];
                        //  wrap for ghost atom
                        if ( dx > boxsize_x * 0.5 )
                            dx = dx - boxsize_x;
                        else if ( dx < - boxsize_x * 0.5 )
                            dx = dx + boxsize_x;
                        if ( dy > boxsize_y * 0.5 )
                            dy = dy - boxsize_y;
                        else if ( dy < - boxsize_y * 0.5 )
                            dy = dy + boxsize_y;
                        if ( dz > boxsize_z * 0.5 )
                            dz = dz - boxsize_z;
                        else if ( dz < - boxsize_z * 0.5 )
                            dz = dz + boxsize_z;
                        // calculate distance square
                        d_square = dx*dx + dy*dy + dz*dz;

                        if ( d_square <= q_cutoff_square )
                        {
                            Vector_destin_index.push_back(tmp_AtomIndex);
                            Vector_destin_r2.push_back(d_square);
                        }
                    }
                }
            }
            else    // for atoms in the center bin, we don't have to worry about whether they are ghost atoms
            {
                for (int j = 0; j < (int)bins_atom_index_[ bins_neigh_bins_[center_BinIndex][i] ].size(); j++)
                {
                    tmp_AtomIndex = bins_atom_index_[ bins_neigh_bins_[center_BinIndex][i] ][j];
                    if ( !ExcludeNeighByTopology(cur_step, center_AtomIndex, tmp_AtomIndex) )
                    {
                        dx = cur_step.Atoms[tmp_AtomIndex][X] - cur_step.Atoms[center_AtomIndex][X];
                        dy = cur_step.Atoms[tmp_AtomIndex][Y] - cur_step.Atoms[center_AtomIndex][Y];
                        dz = cur_step.Atoms[tmp_AtomIndex][Z] - cur_step.Atoms[center_AtomIndex][Z];
                        // calculate distance square
                        d_square = dx*dx + dy*dy + dz*dz;

                        if ( d_square <= q_cutoff_square && tmp_AtomIndex != center_AtomIndex )
                        {
                            Vector_destin_index.push_back(tmp_AtomIndex);
                            Vector_destin_r2.push_back(d_square);
                        }
                    }
                }
            }
        }
    }

    return 0;
}

int Neigh::BuildNeighborList(const STEP &cur_step, const double cutoff, const vector<int> &center_atoms_flag, const vector<int> &neighbor_atoms_flag, \
vector<int> &Vector_destin_count, vector< vector<int> > &Vector_destin_index, vector< vector<double> > &Vector_destin_r2)
{
    if ( cur_step.Atoms.size() != center_atoms_flag.size() || cur_step.Atoms.size() != neighbor_atoms_flag.size() )
    {
        fprintf(stderr, "Input vector<int> as the mask of center and neighbor atoms should both have the same number of element as the number of atoms in current step\n");
        return 1;
    }

    // prepare bin pair list, and assign atoms to bins
    if( Prepare(cur_step,cutoff,neighbor_atoms_flag) )
        return 1;
    vector<int> neigh_count(cur_step.NumOfAtoms);
    vector< vector<int> > neigh_index(cur_step.NumOfAtoms);
    vector< vector<double> > neigh_r2(cur_step.NumOfAtoms);

    #pragma omp parallel /*num_threads(NUM_THREADS)*/
    {
        vector<int> atom_neigh_index;
        vector<double> atom_neigh_r2;    // square of distance to neighbor atoms
        #pragma omp for
        for (int ai = 0; ai < cur_step.NumOfAtoms; ai++)
        {
            if ( center_atoms_flag[ai] == 1 )
            {
                Query(cur_step,cutoff,neighbor_atoms_flag,ai,atom_neigh_index,atom_neigh_r2);
            }
            else
            {
                vector<int>().swap(atom_neigh_index);
                vector<double>().swap(atom_neigh_r2);
            }
            neigh_count[ai] = (int)atom_neigh_index.size();
            neigh_index[ai] = atom_neigh_index;
            neigh_r2[ai] = atom_neigh_r2;
        }
    }

    vector<int>().swap(Vector_destin_count);
    vector< vector<int> >().swap(Vector_destin_index);
    vector< vector<double> >().swap(Vector_destin_r2);
    Vector_destin_count = neigh_count;
    Vector_destin_index = neigh_index;
    Vector_destin_r2 = neigh_r2;

    return 0;
}

int Neigh::Cluster(const STEP &cur_step, const double cutoff, const vector<int> &cluster_atoms_flag, vector<int> &Vector_destin)
{
    // build neighbor list
    vector<int> neigh_count(cur_step.NumOfAtoms);
    vector< vector<int> > neigh_index(cur_step.NumOfAtoms);
    vector< vector<double> > neigh_r2(cur_step.NumOfAtoms);
    if( BuildNeighborList(cur_step,cutoff,cluster_atoms_flag,cluster_atoms_flag,neigh_count,neigh_index,neigh_r2) )
        return 1;

    // Breath First Search algorithm
    vector<int> clusterid(cur_step.NumOfAtoms, 0);
    vector<int> clustersize(1,0);
    list<int> queue;
    int countcluster = 0;
    int tmpai;
    for (int ai = 0; ai < cur_step.NumOfAtoms; ai++)
    {
        if ( clusterid[ai] == 0 && cluster_atoms_flag[ai] == 1 )
        {
            countcluster++;
            clustersize.push_back(0);
            queue.push_back(ai);
            while ( !queue.empty() )
            {
                tmpai = queue.front();
                queue.pop_front();
                clustersize[countcluster]++;
                clusterid[tmpai] = countcluster;
                for (int i = 0; i < (int)neigh_index[tmpai].size(); i++)
                {
                    if ( clusterid[neigh_index[tmpai][i]] != countcluster )
                    {
                        clusterid[neigh_index[tmpai][i]] = countcluster;
                        queue.push_back(neigh_index[tmpai][i]);
                    }
                }
            }
        }
    }

    // sort by size in decending order
    vector<int> index((int)clustersize.size());
    for (int i = 0; i < (int)index.size(); i++)
        index[i] = i;
    sort(index.begin() + 1, index.end(), [&clustersize](int a, int b){return clustersize[a] > clustersize[b];});
    vector<int> oldid_to_newid((int)clustersize.size());
    for (int i = 1; i < (int)index.size(); i++)
        oldid_to_newid[index[i]] = i;
    for (int ai = 0; ai < cur_step.NumOfAtoms; ai++)
        if ( clusterid[ai] != 0 )
            clusterid[ai] = oldid_to_newid[clusterid[ai]];

    vector<int>().swap(Vector_destin);
    Vector_destin = clusterid;
    //vector<int>().swap(clusterid);

    return 0;
}

int Neigh::RDF(const STEP &cur_step, const double cutoff, const int nbins, const vector<int> &center_atoms_flag, const vector<int> &neighbor_atoms_flag, \
vector<double> &Vector_destin_r, vector<double> &Vector_destin_g)
{
    // build neighbor list
    vector<int> neigh_count(cur_step.NumOfAtoms);
    vector< vector<int> > neigh_index(cur_step.NumOfAtoms);
    vector< vector<double> > neigh_r2(cur_step.NumOfAtoms);
    if( BuildNeighborList(cur_step,cutoff,center_atoms_flag,neighbor_atoms_flag,neigh_count,neigh_index,neigh_r2) )
        return 1;
   
    // set up the vector to record the rdf
    double bin_size = cutoff / nbins;
    double half_bin_size = bin_size * 0.5;
    vector<int> count(nbins, 0);
    vector<double> r(nbins, 0);
    vector<double> g(nbins, 0);
    for (int i = 0; i < nbins; i++)
        r[i] = bin_size * i + half_bin_size;

    // calculate rdf
    double r_tmp;
    int bin_id_tmp;

    int center_atom_count = 0;
    int neighbor_atom_count = 0;
    for (int i = 0; i < (int)center_atoms_flag.size(); i++)
    {
        if ( center_atoms_flag[i] == 1 )
            center_atom_count++;
        if ( neighbor_atoms_flag[i] == 1 )
            neighbor_atom_count++;
    }

    for (int i = 0; i < (int)neigh_r2.size(); i++)
    {
        for (int j = 0; j < (int)neigh_r2[i].size(); j++)
        {
            r_tmp = sqrt(neigh_r2[i][j]);
            bin_id_tmp = (int)(r_tmp / bin_size);
            //printf("cid=%d, nid=%d, r=%.4f, bid=%d\n", (int)cur_step.Atoms[i][cur_step.ID], (int)cur_step.Atoms[neigh_index[i][j]][cur_step.ID], r_tmp, bin_id_tmp);
            if ( bin_id_tmp < nbins )
                count[bin_id_tmp]++;
            else
                printf("r = %f not recorded into rdf.\n", r_tmp);
        }
    }

    double neighbor_atom_density = neighbor_atom_count / (cur_step.Box.xhi - cur_step.Box.xlo) / (cur_step.Box.yhi - cur_step.Box.ylo) / (cur_step.Box.zhi - cur_step.Box.zlo);
    for (int i = 0; i < nbins; i++)
        g[i] = double(count[i]) / neighbor_atom_density / center_atom_count / (4.0 / 3.0 * 3.14159265 * (pow(r[i] + half_bin_size, 3) - pow(r[i] - half_bin_size, 3)));

    vector<double>().swap(Vector_destin_r);
    vector<double>().swap(Vector_destin_g);
    Vector_destin_r = r;
    Vector_destin_g = g;

    return 0;
}

int Neigh::RDFByTypes(const STEP &cur_step, const double cutoff, const int nbins, \
vector< vector<int> > &Vector_destin_pair_index, vector< vector<double> > &Vector_destin_r, vector< vector<double> > &Vector_destin_g)
{
    if (cur_step.MOL == -1)
    {
        fprintf(stderr, "Please assign molecule id first.\n");
        return 1;
    }

    if (cur_step.TYPE == -1)
    {
        fprintf(stderr, "Atom type is required.\n");
        return 1;
    }

    int max_type = 0;
    vector<int> natoms_by_type(1,0);    // natoms_by_type[atom_type] = number of atoms with atom type "atom_type". natoms_by_type[0] = 0 (since no atom type = 0).
    for (int i = 0; i < (int)cur_step.Atoms.size(); i++)
    {
        if (cur_step.Atoms[i][cur_step.TYPE] > max_type)
        {
            max_type = cur_step.Atoms[i][cur_step.TYPE];
            natoms_by_type.resize(max_type+1, 0);
        }
        natoms_by_type[cur_step.Atoms[i][cur_step.TYPE]]++;
    }

    vector< vector<int> > pair_index(max_type + 1, vector<int>(max_type + 1, -1));    // pair_index[type_i][type_j] = pair index (an index of a pair of atom types)
    int num_unique_pairs = 0;
    for (int i = 1; i <= max_type; i++)
        for (int j = i; j <= max_type; j++)
        {
            pair_index[i][j] = num_unique_pairs;
            num_unique_pairs++;
        }
    for (int i = max_type; i >= 1; i--)
        for (int j = i - 1; j >= 1; j--)
            pair_index[i][j] = pair_index[j][i];

    // set up the vector to record the rdf
    double bin_size = cutoff / nbins;
    double half_bin_size = bin_size * 0.5;
    vector< vector<int> > count(num_unique_pairs, vector<int>(nbins, 0));
    vector< vector<double> > r(num_unique_pairs, vector<double>(nbins, 0));
    vector< vector<double> > g(num_unique_pairs, vector<double>(nbins, 0));
    for (int pi = 0; pi < num_unique_pairs; pi++)    // pi is pair index
        for (int i = 0; i < nbins; i++)
            r[pi][i] = bin_size * i + half_bin_size;
    
    // build neighbor list
    vector<int> neigh_count(cur_step.NumOfAtoms);
    vector< vector<int> > neigh_index(cur_step.NumOfAtoms);
    vector< vector<double> > neigh_r2(cur_step.NumOfAtoms);
    vector<int> center_atoms_flag(cur_step.NumOfAtoms, 1);
    vector<int> neighbor_atoms_flag(cur_step.NumOfAtoms, 1);
    if( BuildNeighborList(cur_step, cutoff, center_atoms_flag, neighbor_atoms_flag, neigh_count, neigh_index, neigh_r2) )
        return 1;

    // calculate rdf
    double r_tmp;
    int ni;    // index of a neighbor atom
    int cid, nid;    // id of a center and neighbor atom
    int bin_id_tmp, center_type, neighbor_type;
    for (int i = 0; i < (int)neigh_r2.size(); i++)
    {
        cid = cur_step.Atoms[i][cur_step.ID];
        center_type = cur_step.Atoms[i][cur_step.TYPE];
        for (int j = 0; j < (int)neigh_r2[i].size(); j++)
        {
            ni = neigh_index[i][j];
            nid = cur_step.Atoms[ni][cur_step.ID];
            neighbor_type = cur_step.Atoms[ni][cur_step.TYPE];
            if (center_type <= neighbor_type)
            {
                r_tmp = sqrt(neigh_r2[i][j]);
                bin_id_tmp = (int)(r_tmp / bin_size);
                //printf("cid=%d, nid=%d, r=%.4f, bid=%d\n", (int)cur_step.Atoms[i][cur_step.ID], (int)cur_step.Atoms[neigh_index[i][j]][cur_step.ID], r_tmp, bin_id_tmp);
                if ( bin_id_tmp < nbins )
                    count[pair_index[center_type][neighbor_type]][bin_id_tmp]++;
                else
                    printf("r = %f not recorded into rdf.\n", r_tmp);
            }
        }
    }
    
    for (int i = 1; i <= max_type; i++)
        for (int j = i; j <= max_type; j++)
        {
            int pi = pair_index[i][j];    // pi is pair index
            double neighbor_atom_density = natoms_by_type[j] / (cur_step.Box.xhi - cur_step.Box.xlo) / (cur_step.Box.yhi - cur_step.Box.ylo) / (cur_step.Box.zhi - cur_step.Box.zlo);
            for (int bi = 0; bi < nbins; bi++)    // bi is bin index
                g[pi][bi] = double(count[pi][bi]) / neighbor_atom_density / natoms_by_type[i] / (4.0 / 3.0 * 3.14159265 * (pow(r[pi][bi] + half_bin_size, 3) - pow(r[pi][bi] - half_bin_size, 3)));
        }

    vector< vector<int> >().swap(Vector_destin_pair_index);
    vector< vector<double> >().swap(Vector_destin_r);
    vector< vector<double> >().swap(Vector_destin_g);
    Vector_destin_pair_index = pair_index;
    Vector_destin_r = r;
    Vector_destin_g = g;

    // clear memory
    vector<int>().swap(natoms_by_type);
    vector< vector<int> >().swap(pair_index);
    vector< vector<int> >().swap(count);
    vector< vector<double> >().swap(r);
    vector< vector<double> >().swap(g);
    vector< vector<int> >().swap(neigh_index);
    vector< vector<double> >().swap(neigh_r2);
    vector<int>().swap(center_atoms_flag);
    vector<int>().swap(neighbor_atoms_flag);

    return 0;
}
