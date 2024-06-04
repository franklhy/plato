/******************************************************************************************************/
/***This code is used to analyze the topology for data stored in container STEP.***********************/
/******************************************************************************************************/
/***                                  HOW TO USE THIS CLASS                                         ***/
/*** 1. Build connectivity list:                                                                    ***/
/***                                                                                                ***/
/*** (1) Claim an object of Topo class. No parameter is needed.                                     ***/
/***                                                                                                ***/
/*** (2) Call the following function to build connectivity list for the current step:               ***/
/***     Connectivity(const STEP &cur_step, const bool bond, const bool angle, const bool dihedral, ***/ 
/***         const bool improper, vector< vector<int> > &Vector_destin)                             ***/
/***       This function will build a connectivity list based on bond, angle, dihedral or/and       ***/
/***       improper interaction recorded in "cur_step". Results will be written to "Vector_destin". ***/
/***       For a center atom with index (not atom id) "i", Vector_destin[i] will record atom index  ***/
/***       (not atom id) of atoms connected to it through bond, angle, dihedral or/and improper.    ***/
/***       Booleans "bond", "angle", "dihedral" and "improper" specify whether to consider these    ***/
/***       interactions.                                                                            ***/
/***                                                                                                ***/
/*** 2. To perform cluster analysis based on bond connectivity:                                     ***/
/***                                                                                                ***/
/*** (1) Claim an object of Topo class. No parameter is needed.                                     ***/
/***                                                                                                ***/
/*** (2) Call the following function to assign cluster id for all atoms in current step:            ***/
/***     ClusterAnalysis(const STEP &cur_step, vector<int> &Vector_destin)                          ***/
/***       This function will analyze the connectivity of atoms in "cur_step" based on the bond     ***/
/***       information stored in "cur_step". The result will be written into "Vector_destin". For   ***/
/***       an atom with atom index (not atom id) "i", the id of cluster which it belongs to is      ***/
/***       Vector_destin[i]. Each cluster is a group of atoms connected together through bonds.     ***/
/***       Clusters are sorted according to the size in descending order, i.e. cluster id = 1 is    ***/
/***       assigned to the largest cluster.                                                         ***/
/******************************************************************************************************/
/*** Major revision:                                                                                ***/
/*** (1) Modify: Change the signature of function Connectivity.                                     ***/
/*** (2) Modify: Change the function ClusterAnalysis, so it will not modify input STEP object.      ***/
/***             Instead, the cluster id of each atom will be recorded into a int vector.           ***/
/******************************************************************************************************/
/*** Version: 12/10/2019 (1.1 Version)  ***/
/*** By: Heyi Liang                     ***/
/******************************************/

#ifndef TOPO_H
#define TOPO_H

#define TOPO_VERSION "1.1"

#include "Container.h"

class Topo
{
    private:
        vector< vector<int> > ConnectList_;
        int AppendIfNotFound(const int target, vector<int> &array);
        int BuildConnectivityList(const STEP &cur_step, const bool bond, const bool angle, const bool dihedral, const bool improper);

    public:
        Topo();
        ~Topo();
        void print_version();

        int Connectivity(const STEP &cur_step, const bool bond, const bool angle, const bool dihedral, const bool improper, vector< vector<int> > &Vector_destin);
        int ClusterAnalysis(const STEP &cur_step, vector<int> &Vector_destin);
};

#endif // TOPO_H
