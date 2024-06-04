/***********************************************************/
/*** Refer to header file for instructions of this class ***/
/***********************************************************/
/*** Version: 11/23/2016 (1.1 Version) ***/
/*** By: Heyi Liang                    ***/
/*****************************************/

#include "WriteData.h"

WriteData::WriteData()
{
}

WriteData::~WriteData()
{
}

void WriteData::print_version()
{
    printf("Running WriteData.cpp, version %s\n", WRITEDATA_VERSION);
}

/** public function that is used to write file **/
int WriteData::Write(string FileName, const STEP& cur_step, const char* atom_style)
{
    FILE *pFile;

    // test whether the file already exist
    pFile = fopen(FileName.c_str(), "r");
    if ( pFile != NULL )
    {
/*        fprintf(stderr, "File already exist, do you want to overwrite it? (y or n)\n");
        int test = getchar();
        if ( test == 'y' )
        {
            fclose(pFile);
            fprintf(stderr, "Overwriting...\n");
        }
        else
        {
            fprintf(stderr, "Quit...\n");
            fclose(pFile);
            exit(1);
        }
*/
        fprintf(stderr,"File already exist will be overwritten.\n");
    }

    // create new file.
    pFile = fopen(FileName.c_str(), "w");
    if ( pFile == NULL )
    {
        fprintf(stderr, "Can not create the file!\n");
        exit(1);
    }


    // write title
    fprintf(pFile, "#LAMMPS data file via WriteData version %s\n\n", WRITEDATA_VERSION);

    // write header
    if ( cur_step.NumOfAtoms > 0 )
        fprintf(pFile,"%d atoms\n", cur_step.NumOfAtoms);
    if ( cur_step.NumOfAtoms > 0 && cur_step.NumOfAtomTypes > 0 )
        fprintf(pFile,"%d atom types\n", cur_step.NumOfAtomTypes);
    if ( cur_step.NumOfAtoms > 0 && cur_step.NumOfBonds > 0 )
        fprintf(pFile,"%d bonds\n", cur_step.NumOfBonds);
    if ( cur_step.NumOfAtoms > 0 && cur_step.NumOfBonds > 0 && cur_step.NumOfBondTypes > 0 )
        fprintf(pFile,"%d bond types\n", cur_step.NumOfBondTypes);
    if ( cur_step.NumOfAtoms > 0 && cur_step.NumOfAngles > 0 )
        fprintf(pFile,"%d angles\n", cur_step.NumOfAngles);
    if ( cur_step.NumOfAtoms > 0 && cur_step.NumOfAngles > 0 && cur_step.NumOfAngleTypes > 0 )
        fprintf(pFile,"%d angle types\n", cur_step.NumOfAngleTypes);
    if ( cur_step.NumOfAtoms > 0 && cur_step.NumOfDihedrals > 0 )
        fprintf(pFile,"%d dihedrals\n", cur_step.NumOfDihedrals);
    if ( cur_step.NumOfAtoms > 0 && cur_step.NumOfDihedrals > 0 && cur_step.NumOfDihedralTypes > 0 )
        fprintf(pFile,"%d dihedral types\n", cur_step.NumOfDihedralTypes);
    if ( cur_step.NumOfAtoms > 0 && cur_step.NumOfImpropers > 0 )
        fprintf(pFile,"%d impropers\n", cur_step.NumOfImpropers);
    if ( cur_step.NumOfAtoms > 0 && cur_step.NumOfImpropers > 0 && cur_step.NumOfImproperTypes > 0 )
        fprintf(pFile,"%d improper types\n", cur_step.NumOfImproperTypes);

    fprintf(pFile, "\n%-1.16e %-1.16e xlo xhi\n", cur_step.Box.xlo, cur_step.Box.xhi);
    fprintf(pFile,   "%-1.16e %-1.16e ylo yhi\n", cur_step.Box.ylo, cur_step.Box.yhi);
    fprintf(pFile,   "%-1.16e %-1.16e zlo zhi\n", cur_step.Box.zlo, cur_step.Box.zhi);
    if (cur_step.Box.triclinic) {
        fprintf(pFile, "%-1.16e %-1.16e %-1.16e xy xz yz\n\n", cur_step.Box.xy, cur_step.Box.xz, cur_step.Box.yz);
    } else {
        fprintf(pFile, "\n");
    }

    if ( strcmp(atom_style,"atomic") != 0 && strcmp(atom_style,"molecular") != 0 && strcmp(atom_style,"sphere") != 0 && strcmp(atom_style,"full") != 0 )
    {
        fprintf(stderr, "Invalid atom style for Write command.\n");
        fclose(pFile);
        return 1;
    }

    // write body
    // write masses section
    if ( (int)cur_step.Masses.size() > 1 )
    {
        fprintf(pFile, "Masses\n\n");
        for (int i = 1; i < (int)cur_step.Masses.size(); i++)
            fprintf(pFile, "%d %g\n", i, cur_step.Masses[i]);
        fprintf(pFile, "\n");
    }

    // write atoms section and velocities section
    if ( cur_step.NumOfAtoms > 0 && (int)cur_step.Atoms.size() == cur_step.NumOfAtoms )
    {
        // atoms section
        fprintf(pFile, "Atoms # %s\n\n", atom_style);
        if ( strcmp(atom_style,"atomic") == 0 )
        {
            if ( cur_step.ID != -1 && cur_step.TYPE != -1 && cur_step.X != -1 && cur_step.Y != -1 && cur_step.Z != -1 )
            {
                if ( cur_step.IX != -1 && cur_step.IY != -1 && cur_step.IZ != -1 )
                    for (int i = 0; i < (int)cur_step.Atoms.size(); i++)
                        fprintf(pFile, "%d %d %-1.16e %-1.16e %-1.16e %d %d %d\n", \
                        (int)cur_step.Atoms[i][cur_step.ID], (int)cur_step.Atoms[i][cur_step.TYPE], \
                        cur_step.Atoms[i][cur_step.X], cur_step.Atoms[i][cur_step.Y], cur_step.Atoms[i][cur_step.Z], \
                        (int)cur_step.Atoms[i][cur_step.IX], (int)cur_step.Atoms[i][cur_step.IY], (int)cur_step.Atoms[i][cur_step.IZ]);
                else
                    for (int i = 0; i < (int)cur_step.Atoms.size(); i++)
                        fprintf(pFile, "%d %d %-1.16e %-1.16e %-1.16e\n", \
                        (int)cur_step.Atoms[i][cur_step.ID], (int)cur_step.Atoms[i][cur_step.TYPE], \
                        cur_step.Atoms[i][cur_step.X], cur_step.Atoms[i][cur_step.Y], cur_step.Atoms[i][cur_step.Z]);
            }
            else
            {
                fprintf(stderr,"No enough properties in your current step.\n");
                fclose(pFile);
                return 1;
            }
        }
        else if ( strcmp(atom_style,"molecular") == 0 )
        {
            if ( cur_step.ID != -1 && cur_step.MOL != -1 && cur_step.TYPE != -1 && cur_step.X != -1 && cur_step.Y != -1 && cur_step.Z != -1 )
            {
                if ( cur_step.IX != -1 && cur_step.IY != -1 && cur_step.IZ != -1 )
                    for (int i = 0; i < (int)cur_step.Atoms.size(); i++)
                        fprintf(pFile, "%d %d %d %-1.16e %-1.16e %-1.16e %d %d %d\n", \
                        (int)cur_step.Atoms[i][cur_step.ID], (int)cur_step.Atoms[i][cur_step.MOL], (int)cur_step.Atoms[i][cur_step.TYPE], \
                        cur_step.Atoms[i][cur_step.X], cur_step.Atoms[i][cur_step.Y], cur_step.Atoms[i][cur_step.Z], \
                        (int)cur_step.Atoms[i][cur_step.IX], (int)cur_step.Atoms[i][cur_step.IY], (int)cur_step.Atoms[i][cur_step.IZ]);
                else               
                    for (int i = 0; i < (int)cur_step.Atoms.size(); i++)
                        fprintf(pFile, "%d %d %d %-1.16e %-1.16e %-1.16e\n", \
                        (int)cur_step.Atoms[i][cur_step.ID], (int)cur_step.Atoms[i][cur_step.MOL], (int)cur_step.Atoms[i][cur_step.TYPE], \
                        cur_step.Atoms[i][cur_step.X], cur_step.Atoms[i][cur_step.Y], cur_step.Atoms[i][cur_step.Z]);
            }
            else
            {   
                fprintf(stderr,"No enough properties in your current step.\n");
                fclose(pFile);
                return 1;
            }
        }
        else if ( strcmp(atom_style,"sphere") == 0 )
        {
            if ( cur_step.ID != -1 && cur_step.TYPE != -1 && cur_step.DIAMETER != -1 && cur_step.DENSITY != -1 && cur_step.X != -1 && cur_step.Y != -1 && cur_step.Z != -1)
            {
                if ( cur_step.IX != -1 && cur_step.IY != -1 && cur_step.IZ != -1 )
                    for (int i = 0; i < (int)cur_step.Atoms.size(); i++)
                        fprintf(pFile, "%d %d %f %f %-1.16e %-1.16e %-1.16e %d %d %d\n", \
                        (int)cur_step.Atoms[i][cur_step.ID], (int)cur_step.Atoms[i][cur_step.TYPE], cur_step.Atoms[i][cur_step.DIAMETER], cur_step.Atoms[i][cur_step.DENSITY],\
                        cur_step.Atoms[i][cur_step.X], cur_step.Atoms[i][cur_step.Y], cur_step.Atoms[i][cur_step.Z], \
                        (int)cur_step.Atoms[i][cur_step.IX], (int)cur_step.Atoms[i][cur_step.IY], (int)cur_step.Atoms[i][cur_step.IZ]);
                else               
                    for (int i = 0; i < (int)cur_step.Atoms.size(); i++)
                        fprintf(pFile, "%d %d %f %f %-1.16e %-1.16e %-1.16e\n", \
                        (int)cur_step.Atoms[i][cur_step.ID], (int)cur_step.Atoms[i][cur_step.TYPE], cur_step.Atoms[i][cur_step.DIAMETER], cur_step.Atoms[i][cur_step.DENSITY],\
                        cur_step.Atoms[i][cur_step.X], cur_step.Atoms[i][cur_step.Y], cur_step.Atoms[i][cur_step.Z]);
            }
            else
            {   
                fprintf(stderr,"No enough properties in your current step.\n");
                fclose(pFile);
                return 1;
            }
        }
        else if ( strcmp(atom_style,"full") == 0 )
        {
            if ( cur_step.ID != -1 && cur_step.MOL != -1 && cur_step.TYPE != -1 && cur_step.X != -1 && cur_step.Y != -1 && cur_step.Z != -1 && cur_step.Q != -1 )
            {
                if ( cur_step.IX != -1 && cur_step.IY != -1 && cur_step.IZ != -1 )
                    for (int i = 0; i < (int)cur_step.Atoms.size(); i++)
                        fprintf(pFile, "%d %d %d %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n", \
                        (int)cur_step.Atoms[i][cur_step.ID], (int)cur_step.Atoms[i][cur_step.MOL], (int)cur_step.Atoms[i][cur_step.TYPE], cur_step.Atoms[i][cur_step.Q],\
                        cur_step.Atoms[i][cur_step.X], cur_step.Atoms[i][cur_step.Y], cur_step.Atoms[i][cur_step.Z], \
                        (int)cur_step.Atoms[i][cur_step.IX], (int)cur_step.Atoms[i][cur_step.IY], (int)cur_step.Atoms[i][cur_step.IZ]);
                else               
                    for (int i = 0; i < (int)cur_step.Atoms.size(); i++)
                        fprintf(pFile, "%d %d %d %-1.16e %-1.16e %-1.16e %-1.16e\n", \
                        (int)cur_step.Atoms[i][cur_step.ID], (int)cur_step.Atoms[i][cur_step.MOL], (int)cur_step.Atoms[i][cur_step.TYPE], cur_step.Atoms[i][cur_step.Q],\
                        cur_step.Atoms[i][cur_step.X], cur_step.Atoms[i][cur_step.Y], cur_step.Atoms[i][cur_step.Z]);
            }
            else
            {
                fprintf(stderr,"No enough properties in your current step.\n");
                fclose(pFile);
                return 1;
            }
        }
        else
        {
            fprintf(stderr, "Invalid atom style!\n");
            fclose(pFile);
            return 1;
        }
        fprintf(pFile,"\n");

        // velocities section
        if ( strcmp(atom_style,"sphere") == 0 && cur_step.VX != -1 && cur_step.VY != -1 && cur_step.VZ != -1 && cur_step.WX != -1 && cur_step.WY != -1 && cur_step.WZ != -1)
        {
            fprintf(pFile, "Velocities\n\n");
            for(int i = 0; i < (int)cur_step.Atoms.size(); i++)
                fprintf(pFile, "%d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n", (int)cur_step.Atoms[i][cur_step.ID],\
                cur_step.Atoms[i][cur_step.VX], cur_step.Atoms[i][cur_step.VY], cur_step.Atoms[i][cur_step.VZ], cur_step.Atoms[i][cur_step.WX], cur_step.Atoms[i][cur_step.WY], cur_step.Atoms[i][cur_step.WZ]);
            fprintf(pFile,"\n");
        }
        else if ( cur_step.VX != -1 && cur_step.VY != -1 && cur_step.VZ != -1 )
        {
            fprintf(pFile, "Velocities\n\n");
            for(int i = 0; i < (int)cur_step.Atoms.size(); i++)
                fprintf(pFile, "%d %-1.16e %-1.16e %-1.16e\n", (int)cur_step.Atoms[i][cur_step.ID], cur_step.Atoms[i][cur_step.VX], cur_step.Atoms[i][cur_step.VY], cur_step.Atoms[i][cur_step.VZ]);
            fprintf(pFile,"\n");
        }
    }
    else
    {
        fprintf(stderr, "No enough atoms!\n");
        fclose(pFile);
        return 1;
    }

    // write bonds, angles, dihedrals and impropers sections
    if ( strcmp(atom_style,"molecular") == 0 || strcmp(atom_style,"full") == 0 )
    {
        // write bonds section
        if ( cur_step.NumOfBonds == 0 );
        else if ( cur_step.NumOfBonds > 0 && (int)cur_step.Bonds.size() == cur_step.NumOfBonds)
        {
            fprintf(pFile, "Bonds\n\n");
            for(int i = 0; i < (int)cur_step.Bonds.size(); i++)
                fprintf(pFile, "%d %d %d %d\n", cur_step.Bonds[i][0], cur_step.Bonds[i][1], cur_step.Bonds[i][2], cur_step.Bonds[i][3]);
            fprintf(pFile,"\n");
        }
        else
        {
            fprintf(stderr, "No enough bonds!\n");
            fclose(pFile);
            return 1;
        }

        // write angles section
        if ( cur_step.NumOfAngles == 0 );
        else if ( cur_step.NumOfAngles > 0 && (int)cur_step.Angles.size() == cur_step.NumOfAngles)
        {
            fprintf(pFile, "Angles\n\n");
            for(int i = 0; i < (int)cur_step.Angles.size(); i++)
                fprintf(pFile, "%d %d %d %d %d\n", cur_step.Angles[i][0], cur_step.Angles[i][1], cur_step.Angles[i][2], cur_step.Angles[i][3], cur_step.Angles[i][4]);
            fprintf(pFile,"\n");
        }
        else
        {
            fprintf(stderr, "No enough angles!\n");
            fclose(pFile);
            return 1;
        }

        // write dihedrals section
        if ( cur_step.NumOfDihedrals == 0 );
        else if ( cur_step.NumOfDihedrals > 0 && (int)cur_step.Dihedrals.size() == cur_step.NumOfDihedrals)
        {
            fprintf(pFile, "Dihedrals\n\n");
            for(int i = 0; i < (int)cur_step.Dihedrals.size(); i++)
                fprintf(pFile, "%d %d %d %d %d %d\n", cur_step.Dihedrals[i][0], cur_step.Dihedrals[i][1], cur_step.Dihedrals[i][2], cur_step.Dihedrals[i][3], cur_step.Dihedrals[i][4], cur_step.Dihedrals[i][5]);
            fprintf(pFile,"\n");
        }
        else
        {
            fprintf(stderr, "No enough dihedrals!\n");
            fclose(pFile);
            return 1;
        }

        // write impropers section
        if ( cur_step.NumOfImpropers == 0 );
        else if ( cur_step.NumOfImpropers > 0 && (int)cur_step.Impropers.size() == cur_step.NumOfImpropers)
        {
            fprintf(pFile, "Impropers\n\n");
            for(int i = 0; i < (int)cur_step.Impropers.size(); i++)
                fprintf(pFile, "%d %d %d %d %d %d\n", cur_step.Impropers[i][0], cur_step.Impropers[i][1], cur_step.Impropers[i][2], cur_step.Impropers[i][3], cur_step.Impropers[i][4], cur_step.Impropers[i][5]);
            fprintf(pFile,"\n");
        }
        else
        {
            fprintf(stderr, "No enough impropers!\n");
            fclose(pFile);
            return 1;
        }
    }

    fclose(pFile);
    return 0;
}
