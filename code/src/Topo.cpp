/***********************************************************/
/*** Refer to header file for instructions of this class ***/
/***********************************************************/
/*** Version: 12/10/2019 (1.1 Version)  ***/
/*** By: Heyi Liang                     ***/
/******************************************/

#include "Topo.h"

Topo::Topo()
{

}

Topo::~Topo()
{
    vector< vector<int> >().swap(ConnectList_);
}

void Topo::print_version()
{
    printf("Running Topo.cpp, version %s\n", TOPO_VERSION);
}

int Topo::AppendIfNotFound(const int target, vector<int> &array)
{
    for (int i = 0; i < (int)array.size(); i++)
        if ( target == array[i] )
            return 0;
    array.push_back(target);
    return 0;
}

int Topo::BuildConnectivityList(const STEP &cur_step, const bool bond, const bool angle, const bool dihedral, const bool improper)
{
    //vector<int> ConnectNum;
    vector<int> IDToIndexMap;

    // allocate memory for vectors
    //ConnectNum.resize(cur_step.NumOfAtoms, 0);
    int maxID = 0;
    for (int i = 0; i < cur_step.NumOfAtoms; i++)
        if ( cur_step.Atoms[i][cur_step.ID] > maxID )
            maxID = cur_step.Atoms[i][cur_step.ID];
    IDToIndexMap.resize(maxID+1);

    // build the index-to-id-map
    for (int i = 0; i < cur_step.NumOfAtoms; i++)
        IDToIndexMap[ cur_step.Atoms[i][cur_step.ID] ] = i;

    // build connect-list
    vector< vector<int> >().swap(ConnectList_);
    ConnectList_.resize(cur_step.NumOfAtoms,vector<int>());
    int index_a, index_b, index_c, index_d;
    if (bond)
    {
        for (int i = 0; i < cur_step.NumOfBonds; i++)
        {
            index_a = IDToIndexMap[ cur_step.Bonds[i][2] ];
            index_b = IDToIndexMap[ cur_step.Bonds[i][3] ];
            ConnectList_[index_a].push_back(index_b);
            ConnectList_[index_b].push_back(index_a);
        }
    }
    if (angle)
    {
        for (int i = 0; i < cur_step.NumOfAngles; i++)
        {
            index_a = IDToIndexMap[ cur_step.Angles[i][2] ];
            index_b = IDToIndexMap[ cur_step.Angles[i][3] ];
            index_c = IDToIndexMap[ cur_step.Angles[i][4] ];
            AppendIfNotFound(index_b, ConnectList_[index_a]);
            AppendIfNotFound(index_c, ConnectList_[index_a]);
            AppendIfNotFound(index_a, ConnectList_[index_b]);
            AppendIfNotFound(index_c, ConnectList_[index_b]);
            AppendIfNotFound(index_a, ConnectList_[index_c]);
            AppendIfNotFound(index_b, ConnectList_[index_c]);
        }
    }
    if (dihedral)
    {
        for (int i = 0; i < cur_step.NumOfDihedrals; i++)
        {
            index_a = IDToIndexMap[ cur_step.Dihedrals[i][2] ];
            index_b = IDToIndexMap[ cur_step.Dihedrals[i][3] ];
            index_c = IDToIndexMap[ cur_step.Dihedrals[i][4] ];
            index_d = IDToIndexMap[ cur_step.Dihedrals[i][5] ];
            AppendIfNotFound(index_b, ConnectList_[index_a]);
            AppendIfNotFound(index_c, ConnectList_[index_a]);
            AppendIfNotFound(index_d, ConnectList_[index_a]);
            AppendIfNotFound(index_a, ConnectList_[index_b]);
            AppendIfNotFound(index_c, ConnectList_[index_b]);
            AppendIfNotFound(index_d, ConnectList_[index_b]);
            AppendIfNotFound(index_a, ConnectList_[index_c]);
            AppendIfNotFound(index_b, ConnectList_[index_c]);
            AppendIfNotFound(index_d, ConnectList_[index_c]);
            AppendIfNotFound(index_a, ConnectList_[index_d]);
            AppendIfNotFound(index_b, ConnectList_[index_d]);
            AppendIfNotFound(index_c, ConnectList_[index_d]);
        }
    }
    if (improper)
    {
        for (int i = 0; i < cur_step.NumOfImpropers; i++)
        {
            index_a = IDToIndexMap[ cur_step.Impropers[i][2] ];
            index_b = IDToIndexMap[ cur_step.Impropers[i][3] ];
            index_c = IDToIndexMap[ cur_step.Impropers[i][4] ];
            index_d = IDToIndexMap[ cur_step.Impropers[i][5] ];
            AppendIfNotFound(index_b, ConnectList_[index_a]);
            AppendIfNotFound(index_c, ConnectList_[index_a]);
            AppendIfNotFound(index_d, ConnectList_[index_a]);
            AppendIfNotFound(index_a, ConnectList_[index_b]);
            AppendIfNotFound(index_c, ConnectList_[index_b]);
            AppendIfNotFound(index_d, ConnectList_[index_b]);
            AppendIfNotFound(index_a, ConnectList_[index_c]);
            AppendIfNotFound(index_b, ConnectList_[index_c]);
            AppendIfNotFound(index_d, ConnectList_[index_c]);
            AppendIfNotFound(index_a, ConnectList_[index_d]);
            AppendIfNotFound(index_b, ConnectList_[index_d]);
            AppendIfNotFound(index_c, ConnectList_[index_d]);
        }
    }

    // clear memory
    //vector<int>().swap(ConnectNum);
    vector<int>().swap(IDToIndexMap);

    return 0;
}

int Topo::Connectivity(const STEP &cur_step, const bool bond, const bool angle, const bool dihedral, const bool improper, vector< vector<int> > &Vector_destin)
{
    if ( BuildConnectivityList(cur_step, bond, angle, dihedral, improper) )
    {
        fprintf(stderr, "Failed to build connectivity list.\n");
        return 1;
    }

    vector< vector<int> >().swap(Vector_destin);
    Vector_destin = ConnectList_;

    return 0;
}

int Topo::ClusterAnalysis(const STEP &cur_step, vector<int> &Vector_destin)
{
    if ( BuildConnectivityList(cur_step, true, false, false, false) )
    {
        fprintf(stderr, "Failed to build connectivity list.\n");
        return 1;
    }

    vector<int> state(cur_step.NumOfAtoms, 0);
    vector<int> index_next(1,0);
    vector<int> index_cur;
    vector<int> cluster_size(1,0);
    vector<int> cluster_id(cur_step.NumOfAtoms, 0);

    int index;
    int cluster_id_cur = 1;
    int analyzed_num_atom = 0;

    //fprintf(stdout, "Progress:\n");
    //fprintf(stdout, "analyzed/total: %d/%d\n", analyzed_num_atom, cur_step.NumOfAtoms);
    // analyze connectivity
    while ( analyzed_num_atom != cur_step.NumOfAtoms )
    {
        //fprintf(stdout, "analyzed/total: %d/%d\n", analyzed_num_atom, cur_step.NumOfAtoms);
        cluster_size.push_back(0);
        // start a new cluster from one atom
        index = index_next[0];
        state[index] = 1;
        cluster_id[index] = cluster_id_cur;
        cluster_size[cluster_id_cur]++;
        analyzed_num_atom++;
        while ( index_next.size() != 0 )
        {
            index_cur = index_next;
            vector<int>().swap(index_next);
            for (int i = 0; i < (int)index_cur.size(); i++)
            {
                index = index_cur[i];
                for (int j = 0; j < (int)ConnectList_[index].size(); j++)
                {
                    if ( state[ConnectList_[index][j]] == 0 )
                    {
                        state[ConnectList_[index][j]] = 1;
                        cluster_id[ConnectList_[index][j]] = cluster_id_cur;
                        cluster_size[cluster_id_cur]++;
                        analyzed_num_atom++;
                        index_next.push_back(ConnectList_[index][j]);
                    }
                }
            }
        }
        cluster_id_cur++;

        index_next.resize(1);
        for (int i = 0; i < cur_step.NumOfAtoms; i++)
            if ( state[i] != 1 )
            {
                index_next[0] = i;
                break;
            }
    }
    //fprintf(stdout, "analyzed/total: %d/%d\n", analyzed_num_atom, cur_step.NumOfAtoms);

    // sort clusters by size
    //// find the largest size
    int max_size = 0;
    for (int i = 0; i < (int)cluster_size.size(); i++)
        if ( cluster_size[i] > max_size )
            max_size = cluster_size[i];
    //// count the cluster with same size
    vector<int> same_size_total_count(max_size+1, 0);
    for (int i = 1; i < (int)cluster_size.size(); i++)
        same_size_total_count[cluster_size[i]]++;
    vector<int> larger_than_size_total_count(max_size+1, 0);
    for (int i = max_size - 1; i > 0; i--)
        larger_than_size_total_count[i] += larger_than_size_total_count[i+1] + same_size_total_count[i+1];
    //// assign new cluster id
    vector<int> new_cluster_id((int)cluster_size.size(), 0);
    vector<int>().swap(same_size_total_count);
    same_size_total_count.resize(max_size+1, 0);
    for (int i = 1; i < (int)cluster_size.size(); i++)
    {
        new_cluster_id[i] = larger_than_size_total_count[cluster_size[i]] + same_size_total_count[cluster_size[i]] + 1;
        same_size_total_count[cluster_size[i]]++;
    }
    for (int i = 0; i < cur_step.NumOfAtoms; i++)
        cluster_id[i] = new_cluster_id[cluster_id[i]];

    vector<int>().swap(Vector_destin);
    Vector_destin = cluster_id;

    // clear memory
    vector<int>().swap(state);
    vector<int>().swap(index_next);
    vector<int>().swap(index_cur);
    vector<int>().swap(cluster_size);
    vector<int>().swap(same_size_total_count);
    vector<int>().swap(larger_than_size_total_count);
    vector<int>().swap(cluster_id);
    vector<int>().swap(new_cluster_id);

    return 0;
}
