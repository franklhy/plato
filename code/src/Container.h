/*******************************************************************************************************/
/*** Define some prototypes. These will be used in all cpp files, which deal with data in datafiles  ***/
/*** produced by lammps dump command.                                                                ***/
/*******************************************************************************************************/
/*** Major revision:                                                                                 ***/
/*** (1) Add: Add Masses, which is read from data file.                                              ***/
/*******************************************************************************************************/
/*** Version: 11/12/2018 (1.2 Version)  ***/
/*** By: Heyi Liang                     ***/
/******************************************/

#ifndef TYPES
#define TYPES

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

struct BOX
{
    char xboundary[3], yboundary[3], zboundary[3];  // store boundary type('p', 'f', 's' or 'm', see lammps manual), 3 available position are used for: lower boundary, higher boundary, and '\0'
    double xlo, xhi, ylo, yhi, zlo, zhi;
    double xlo_bound, ylo_bound, zlo_bound, xhi_bound, yhi_bound, zhi_bound;
    double xy, xz, yz;
    bool triclinic;
    double lx, ly, lz, vol;
    double h[6], h_inv[6];    // the matrix of linear transformation between Cartesian and lamda coordinate (used for triclinic box), in Voigt notation
    double d[3];              // distance between parallel boundaries (used for triclinic box)

    BOX();
    ~BOX();
    void set_box(bool bound);
};

struct STEP
{
    /* information provided by dump file */
    int Timestep;
    int NumOfAtoms;
    BOX Box;
    vector<double> Masses;        // record mass of each atom type, Masses[AtomType] = mass. Masses[0]=0 is redundant.
    vector<string> Properties;
    vector< vector<double> > Atoms;        // record atom properties, Atoms[AtomIndex][property's index]. ATTENTION: AtomIndex != AtomID 

    int wrapped_flag;           // 0: unwrapped, 2: wrapped (in LAMMPS dump file), 1: truely wrapped. ATTENTION: wrapped coordinate in dump file does not mean TRUELY wrapped coordinate, refers to lammps manual.
    int scaled_flag;            // 0: unscaled, 1:scaled
    int ID,TYPE,MOL,DIAMETER,DENSITY,X,Y,Z,IX,IY,IZ,VX,VY,VZ,WX,WY,WZ,FX,FY,FZ,Q;        // basic properties mapping (X,Y,Z for coordinates, whether they are wrapped or scaled, depends on wrapped_flag and scaled_flag)

    /* additional information provided by data file */
    int NumOfAtomTypes;
    int NumOfBonds, NumOfBondTypes;
    int NumOfAngles, NumOfAngleTypes;
    int NumOfDihedrals, NumOfDihedralTypes;
    int NumOfImpropers, NumOfImproperTypes;
    vector< vector<int> > Bonds;           // record bond properties, including BondID, Bond Type, Bonded Atom1's ID, Bonded Atom2's ID
    vector< vector<int> > Angles;          // record angle properties, including AngleID, Angle Type, Atom1's ID, Atom2's ID, Atom3's ID 
    vector< vector<int> > Dihedrals;       // record dihedral properties, including DihedralID, Dihedral Type, Atom1's ID, Atom2's ID, Atom3's ID, Atom4's ID
    vector< vector<int> > Impropers;       // record improper properties, including ImproperID, Improper Type, Atom1's ID, Atom2's ID, Atom3's ID, Atom4's ID

    STEP();
    ~STEP();
    void x2lamda();
    void lamda2x();
    int wrap();
    int unwrap();
    void flip();
};        // store all the information in one snapshot

#endif // TYPES
