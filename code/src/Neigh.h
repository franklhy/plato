/******************************************************************************************************/
/***This code is used to build atom neighbor list for data stored in container STEP.*******************/
/******************************************************************************************************/
/*** SPECIAL NOTE:                                                                                  ***/
/*** (1) This code is designed for openmp, and it will use all available core during calculation.   ***/
/******************************************************************************************************/
/***                                  HOW TO USE THIS CLASS                                         ***/
/*** 1. Claim an object of Neigh class. No parameter is needed.                                     ***/
/***                                                                                                ***/
/*** 2. Set the topoloy using the following function:                                               ***/
/***    SetTopology(const STEP &cur_step,                                                           ***/
/***                const bool include_inter_molecule, const bool include_intra_molecule,           ***/
/***                const bool exclude_bond, const bool exclude_angle,                              ***/
/***                const bool exclude_dihedral, const bool exclude_improper)                       ***/
/***     The "topology" means whether a pair of atoms are considered in neigh analysis if they:     ***/
/***     (1) belong to the same molecule (according to molecule ID).                                ***/
/***       "include_inter_molecule = true" means atoms on different molecules are considered.       ***/
/***       "include_intra_molecule = true" means atoms on the same molecule are considered.         ***/
/***     (2) are connected by bond, angle, dihedral and/or improper interactions.                   ***/
/***       "exclude_bond/angle/dihedra/improper = true" means atoms connected by the topological    ***/
/***       interactions bond/angle/dihedra/improper are not considered.                             ***/
/***                                                                                                ***/
/***     NOTE: if the topology in cur_step is not changed, there is no need to use this function    ***/
/***       before performing neigh analysis in step 3. However, if the topology in cur_step is      ***/
/***       changed, then user need to reset the topology by SetTopology.                            ***/
/***                                                                                                ***/
/*** 3. Perform neigh analysis. There are many options:                                             ***/
/***   A. Build a neigh list using the following function:                                          ***/
/***     BuildNeighborList(const STEP &cur_step, const double cutoff,                               ***/
/***                       const vector<int> &center_atoms_flag,                                    ***/
/***                       const vector<int> &neighbor_atoms_flag,                                  ***/
/***                       vector< vector<int> > &Vector_destin_index,                              ***/
/***                       vector< vector<double> > &Vector_destin_r2):                             ***/
/***       Build a atom neighbor list for atoms in "cur_step", within a "cutoff" cutoff radius. The ***/
/***       result will be recorded in "Vector_destin". For an atom whose AtomIndex (not AtomID) is  ***/
/***       i in "cur_step", the AtomIndex of its neighbor atoms will be recorded into the array     ***/
/***       Vector_destin_index[i], the square of distance to neighbor atoms will be recored into    ***/   
/***       Vector_destin_r2[i]. The center atoms are selected from the atoms in "cur_step" by the   ***/
/***       mask "center_atoms_flag", whose size is the number of atoms in "cur_step", and an atom   ***/
/***       is selected as center atoms if center_atoms_flag[i] == 1, where i is the AtomIndex. The  ***/
/***       neighbor atoms are selected according to the mask "neighbor_atoms_flag" in similar way.  ***/
/***                                                                                                ***/
/***   B. Build a neighbor list for only one atom/bead in current step:                             ***/
/***     (i) First, call the following function to prepare:                                         ***/
/***     Prepare(const STEP &cur_step, const double cutoff, const vector<int> &neighbor_atoms_flag):***/
/***       Prepare a cell list for center atoms based on the cutoff radius "cutoff" for "cur_step". ***/
/***       The results is recorded to internal of the class. Refer to (2) for the selection of      ***/
/***       neighbor atoms.                                                                          ***/
/***       NOTE: This function only needs to be called once if no new snapshot is recored into      ***/
/***             "cur_step". If "cur_step" is updated for a new snapshot, this function should be   ***/
/***             called to build a new cell list.                                                   ***/
/***                                                                                                ***/
/***     (ii) Then, call the following function to get the atom neighbor list of one atom:          ***/
/***     Query(const STEP &cur_step, const double q_cutoff, const vector<int> &neighbor_atoms_flag, ***/
/***           const int &center_AtomIndex, vector<int> &Vector_destin_index,                       ***/
/***           vector<double> &Vector_destin_r2):                                                   ***/
/***       Build an atom neighbor list of one atom in "cur_step". The AtomIndex of the center atom  ***/
/***       is "center_AtomIndex". The cutoff radius "q_cutoff" should not be larger than the cutoff ***/
/***       radius offered to function "Prepare" previously. The AtomIndex of the neighbor atoms is  ***/
/***       recorded into array Vector_destin_index. The square of distance to neighbor atoms is     ***/
/***       recorded into Vector_destin_r2. Refer to (2) for the selection of neighbor atoms.        ***/
/***       NOTE: The function Query can be called multiple times to build atom neighbor list for    ***/
/***             different atoms. To achieve the best performance, the cutoff radius should be the  ***/
/***             same as the one used in function Prepare.                                          ***/
/***                                                                                                ***/
/***   C. Perform cluster analysis, i.e. find clusters according to neighboring distance criterion, ***/
/***   using the following function:                                                                ***/
/***     Cluster(const STEP &cur_step, const double cutoff, const vector<int> &cluster_atoms_flag,  ***/
/***             vector<int> &Vector_destin)                                                        ***/
/***       Group atoms in "cur_step" into clusters based on the "cutoff" distance. The result is    ***/
/***       written into "Vector_destin" that Vector_destin[AtomIndex] = ClusterID. Clusters are     ***/
/***       sorted by size in descending order, so ClusterID 1 refers to the largest order. Atoms to ***/
/***       be anaylzed can be selected according the mask "cluster_atoms_flag", where cluster_atoms ***/
/***       _flag[AtomIndex] == 1 means seleted, 0 means not seleted. When an atom is not selected   ***/
/***       for cluster analysis, its ClusterID is assigned to 0.                                    ***/
/***                                                                                                ***/
/***   D. Calculate the radial distribution function (rdf) by using one of the following functions: ***/
/***     RDF(const STEP &cur_step, const double cutoff, const int nbins,                            ***/
/***         const vector<int> &center_atoms_flag, const vector<int> &neighbor_atoms_flag,          ***/
/***         vector<double> &Vector_destin_r, vector<double> &Vector_destin_g)                      ***/
/***       Analyze the radial distribution function (rdf) for atoms in "cur_step". The rdf will be  ***/
/***       analyzed within a cutoff distance set by "cutoff" and shown in histograms with "nbins"   ***/
/***       bins. The histogram is recored in "Vector_destin_r" (recording r) and "Vector_destin_g"  ***/
/***       (recording g(r)). Distance r in "Vector_destin_r" is the center of each histogram bin.   ***/
/***       User should also specify the center atoms and neighbor atoms using "center_atoms_flag"   ***/
/***       and "neighbor_atoms_flag". Setting flag[AtomIndex] = 1 means the atom is selected,       ***/
/***       where flag can be "center_atoms_flag" or "neighbor_atoms_flag", and AtomIndex is the     ***/
/***       atom index (not the atom id) of a atom. Similarly, setting flag[AtomIndex] = 0 means the ***/
/***       atom is not selected.                                                                    ***/
/***                                                                                                ***/
/***     RDFbyTypes(const STEP &cur_step, const double cutoff, const int nbins,                     ***/
/***         vector< vector<int> > &Vector_destin_pair_index,                                       ***/
/***         vector< vector<double> > &Vector_destin_r, vector< vector<double> > &Vector_destin_g)  ***/
/***       Calculate the partial rdf of all pairs of atom types. Input parameters 'cur_step',       ***/
/***       'cutoff' and 'nbins' are the same as those in the function RDF(). The 2d vector          ***/
/***       'Vector_destin_pair_index' records the index of pairs of atom types, of which the r and  ***/
/***       g(r) are recorded in the 2d vector 'Vector_destin_r' and 'Vector_destin_g'.              ***/
/***       e.g. partial rdf of atom with types 1 and 2 are:                                         ***/
/***               r    = Vector_destin_r[Vector_destin_pair_index[1][2]]                           ***/
/***               g(r) = Vector_destin_g[Vector_destin_pair_index[1][2]]                           ***/
/******************************************************************************************************/
/*** Major revision:                                                                                ***/
/*** 1. Chanage the functions for getting radial distribution functions.                            ***/
/******************************************************************************************************/
/*** Version: 07/22/2022 (1.6 Version)  ***/
/*** By: Heyi Liang                     ***/
/******************************************/

#ifndef NEIGH_H
#define NEIGH_H

/*
#ifdef _OPENMP
#define NUM_THREADS 12
#else
#define NUM_THREADS 1
#endif
*/
#define NEIGH_VERSION "1.6"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <list>
#include <algorithm>
#include "Container.h"
#include "Topo.h"

class Neigh
{
    private:
        /* variable for preparing cell list */
        double cutoff_;
        bool large_cutoff_;    // when cutoff > 1/3 box length, large_cutoff_ = true, otherwise large_cutoff_ = false
        int numbinx_, numbiny_, numbinz_, NumOfBins_;
        double binsizex_, binsizey_, binsizez_;
        vector< vector<int> > bins_atom_index_;    // atom index of atoms in each bin
        vector< vector<int> > bins_neigh_bins_;    // indices of neighbor bins of each bin

        /* variable for calculating RDF */
        bool include_inter_molecule_;
        bool include_intra_molecule_;
        vector< vector<int> > connect_list_;

        
        bool FindIn(const int target, const vector<int> &array);    // find a target in an array
        bool ExcludeNeighByTopology(const STEP &cur_step, int center_index, int neigh_index);   // determin whether to exclude the neighbor bead due to topology setting

    public:
        Neigh();
        ~Neigh();
        void print_version();

        int SetTopology(const STEP &cur_step, const bool include_inter_molecule, const bool include_intra_molecule, \
            const bool exclude_bond, const bool exclude_angle, const bool exclude_dihedral, const bool exclude_improper);
        int Prepare(const STEP &cur_step, const double p_cutoff, const vector<int> &neighbor_atoms_flag);
        int Query(const STEP &cur_step, const double q_cutoff, const vector<int> &neighbor_atoms_flag, const int &center_AtomIndex, \
            vector<int> &Vector_destin_index, vector<double> &Vector_destin_r2);
        int BuildNeighborList(const STEP &cur_step, const double cutoff, const vector<int> &center_atoms_flag, const vector<int> &neighbor_atoms_flag, \
            vector< vector<int> > &Vector_destin_index, vector< vector<double> > &Vector_destin_r2);
        
        /* functions related to analyzing clustering based on cutoff distance */
        int Cluster(const STEP &cur_step, const double cutoff, const vector<int> &cluster_atoms_flag, vector<int> &Vector_destin);

        /* functions related to calculating radial distribution function (RDF) */
        int RDF(const STEP &cur_step, const double cutoff, const int nbins, const vector<int> &center_atoms_flag, const vector<int> &neighbor_atoms_flag, \
            vector<double> &Vector_destin_r, vector<double> &Vector_destin_g);
        int RDFByTypes(const STEP &cur_step, const double cutoff, const int nbins, \
            vector< vector<int> > &Vector_destin_pair_index, vector< vector<double> > &Vector_destin_r, vector< vector<double> > &Vector_destin_g);
};

#endif // NEIGH_H
