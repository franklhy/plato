/***********************************************************/
/*** Refer to header file for instructions of this class ***/
/***********************************************************/
/*** Version: 11/12/2016 (1.1 Version)  ***/
/*** By: Heyi Liang                     ***/
/******************************************/

#include "ReadData.h"

#define MAXLINE 1024

enum {ATOMIC,MOLECULAR,SPHERE,FULL};    // atom style 
int NumCol[4] = {5,6,7,7};         // number of column in Atoms section, with corresponding atom style. 3 image flag not included
int NumVel[4] = {3,3,6,3};         // number of column in Velocity section (except the first column, i.e. atom id), with corresponding atom style

ReadData::ReadData()
{
    line = new char[MAXLINE];
}

ReadData::~ReadData()
{
    delete[] line;
}

void ReadData::print_version()
{
    printf("Running ReadData.cpp, version %s\n", READDATA_VERSION);
}

/** functions that are only called when reading the file, they all called by function Read **/
int ReadData::Init()
{
    FilePosition = 0;
    massheadpos = masstailpos = atomheadpos = atomtailpos = velocityheadpos = velocitytailpos = bondheadpos = bondtailpos = angleheadpos = angletailpos = dihedralheadpos = dihedraltailpos = improperheadpos = impropertailpos = 0;
    massflag = atomflag = imageflag = velocityflag = bondflag = angleflag = dihedralflag = improperflag = 0;
    style = -1;

    return 0;
}

int ReadData::ReadHeader(long &FilePos, STEP& step_destin)
{
    // check whether first line start with "LAMMPS data file via write_data"
    // if yes, record the timestep, otherwise set timestep = 0
    fgets(line, MAXLINE, pFile);
    if ( strstr(line,"LAMMPS data file via write_data") )
        step_destin.Timestep = atoi ( strstr(line,"timestep = ") + strlen("timestep = ") );
    else
    {
        step_destin.Timestep = 0;
        //fseek(pFile, 0, SEEK_SET);
    }

    step_destin.Box.triclinic = false;
    // find keywords ("atoms", "atom types", "bonds", "bond types", ... ) and record related informations
    while( fgets(line, MAXLINE, pFile) )
    {
        if ( line[0] != '#' && line[ strspn(line," \t\r") ] != '#' && line[ strspn(line," \t\r") ] != '\n' )    // skip all the comment lines (start with "#")
        {
            if ( strstr(line,"atoms") )
                sscanf(line,"%d",&step_destin.NumOfAtoms);
            else if ( strstr(line,"atom types") )
                sscanf(line,"%d",&step_destin.NumOfAtomTypes);
            else if ( strstr(line,"bonds") )
                sscanf(line,"%d",&step_destin.NumOfBonds);
            else if ( strstr(line,"bond types") )
                sscanf(line,"%d",&step_destin.NumOfBondTypes);
            else if ( strstr(line,"angles") )
                sscanf(line,"%d",&step_destin.NumOfAngles);
            else if ( strstr(line,"angle types") )
                sscanf(line,"%d",&step_destin.NumOfAngleTypes);
            else if ( strstr(line,"dihedrals") )
                sscanf(line,"%d",&step_destin.NumOfDihedrals);
            else if ( strstr(line,"dihedral types") )
                sscanf(line,"%d",&step_destin.NumOfDihedralTypes);
            else if ( strstr(line,"impropers") )
                sscanf(line,"%d",&step_destin.NumOfImpropers);
            else if ( strstr(line,"improper types") )
                sscanf(line,"%d",&step_destin.NumOfImproperTypes);
            else if ( strstr(line,"xlo xhi") )
                sscanf(line,"%lg %lg",&step_destin.Box.xlo,&step_destin.Box.xhi);
            else if ( strstr(line,"ylo yhi") )
                sscanf(line,"%lg %lg",&step_destin.Box.ylo,&step_destin.Box.yhi);
            else if ( strstr(line,"zlo zhi") )
                sscanf(line,"%lg %lg",&step_destin.Box.zlo,&step_destin.Box.zhi);
            else if ( strstr(line,"xy xz yz") ) {
                sscanf(line,"%lg %lg %lg",&step_destin.Box.xy,&step_destin.Box.xz,&step_destin.Box.yz);
                step_destin.Box.triclinic = true;
            }
            /* additional keywords can be added hear */
            else
            {
                fseek(pFile, FilePos, SEEK_SET);
                break;
            }
        }
        FilePos = ftell(pFile);
    }
    if (!step_destin.Box.triclinic) {
        step_destin.Box.xy = 0;
        step_destin.Box.xz = 0;
        step_destin.Box.yz = 0;
    }
    step_destin.Box.set_box(false);

    return 0;
}

int ReadData::FindSections(long &FilePos, STEP& step_destin)
{
    char *keyword;

    fseek(pFile, FilePos, SEEK_SET);
    while ( fgets(line, MAXLINE, pFile) )
    {
        if ( line[0] != '#' && line[ strspn(line," \t\r") ] != '#' && line[ strspn(line," \t\r") ] != '\n' )    // skip all the comment lines (start with "#")
        {
            keyword = strtok(line, "# \t\r\n");
            // whether the first word of the line is one of the keywords: "Atoms", "Velocities, "Bonds", ......
            if ( strcmp(keyword,"Masses") == 0 )
            {
                massheadpos = FilePos;
                massflag = 1;
                for (int i = 0; i < step_destin.NumOfAtomTypes + 1; i++)
                    fgets(line, MAXLINE, pFile);
                masstailpos = ftell(pFile);
            }
            else if ( strcmp(keyword,"Pair") == 0 )
            {
                // do not read pair coeffs
                for (int i = 0; i < step_destin.NumOfAtomTypes + 1; i++)
                    fgets(line, MAXLINE, pFile);
            }
            else if ( strcmp(keyword,"Bond") == 0 )
            {
                // do not read bond coeffs
                for (int i = 0; i < step_destin.NumOfBondTypes + 1; i++)
                    fgets(line, MAXLINE, pFile);
            }
            else if ( strcmp(keyword,"Angle") == 0 )
            {
                // do not read angle coeffs
                for (int i = 0; i < step_destin.NumOfAngleTypes + 1; i++)
                    fgets(line, MAXLINE, pFile);
            }
            else if ( strcmp(keyword,"Dihedral") == 0 )
            {
                // do not read dihedral coeffs
                for (int i = 0; i < step_destin.NumOfDihedralTypes + 1; i++)
                    fgets(line, MAXLINE, pFile);
            }
            else if ( strcmp(keyword,"Improper") == 0 )
            {
                // do not read improper coeffs
                for (int i = 0; i < step_destin.NumOfImproperTypes + 1; i++)
                    fgets(line, MAXLINE, pFile);
            }
            else if ( strcmp(keyword,"Atoms") == 0 )
            {
                atomheadpos = FilePos;
                atomflag = 1;
                if ( step_destin.NumOfAtoms == 0 )
                {
                    fprintf(stderr, "Invalid Atoms section.\n");
                    return 1;
                }
                for (int i = 0; i < step_destin.NumOfAtoms + 1; i++)
                    fgets(line, MAXLINE, pFile);
                atomtailpos = ftell(pFile);
            }
            else if ( strcmp(keyword,"Velocities") == 0 )
            {
                velocityheadpos = FilePos;
                velocityflag = 1;
                if ( step_destin.NumOfAtoms == 0 )
                {
                    fprintf(stderr, "Invalid Velocities section.\n");
                    return 1;
                }
                for (int i = 0; i < step_destin.NumOfAtoms + 1; i++)
                    fgets(line, MAXLINE, pFile);
                velocitytailpos = ftell(pFile);
            }
            else if ( strcmp(keyword,"Bonds") == 0 )
            {
                bondheadpos = FilePos;
                bondflag = 1;
                if ( step_destin.NumOfBonds == 0 )
                {
                    fprintf(stderr, "Invalid Bonds section.\n");
                    return 1;
                }
                for (int i = 0; i < step_destin.NumOfBonds + 1; i++)
                    fgets(line, MAXLINE, pFile);
                bondtailpos = ftell(pFile);
            }
            else if ( strcmp(keyword,"Angles") == 0 )
            {
                angleheadpos = FilePos;
                angleflag = 1;
                if ( step_destin.NumOfAngles == 0 )
                {
                    fprintf(stderr, "Invalid Angles section.\n");
                    return 1;
                }
                for (int i = 0; i < step_destin.NumOfAngles + 1; i++)
                    fgets(line, MAXLINE, pFile);
                angletailpos = ftell(pFile);
            }
            else if ( strcmp(keyword,"Dihedrals") == 0 )
            {
                dihedralheadpos = FilePos;
                dihedralflag = 1;
                if ( step_destin.NumOfDihedrals == 0 )
                {
                    fprintf(stderr, "Invalid Dihedrals section.\n");
                    return 1;
                }
                for (int i = 0; i < step_destin.NumOfDihedrals + 1; i++)
                    fgets(line, MAXLINE, pFile);
                dihedraltailpos = ftell(pFile);
            }
            else if ( strcmp(keyword,"Impropers") == 0 )
            {
                improperheadpos = FilePos;
                improperflag = 1;
                if ( step_destin.NumOfImpropers == 0 )
                {
                    fprintf(stderr, "Invalid Impropers section.\n");
                    return 1;
                }
                for (int i = 0; i < step_destin.NumOfImpropers + 1; i++)
                    fgets(line, MAXLINE, pFile);
                impropertailpos = ftell(pFile);
            }
            else
            {
                PrintTextNearWrongFormat(ftell(pFile));
                return 1;
            }
        }
        FilePos = ftell(pFile);
    }

    return 0;
}

int ReadData::ReadMasses(long &HeadPos, long &TailPos, STEP& step_destin)
{
    int tmpint, count;
    double tmpdouble;

    clearerr(pFile);
    fseek(pFile, HeadPos, SEEK_SET);
    fgets(line, MAXLINE, pFile);

    // one whitespace line is required after the first line of one section
    fgets(line, MAXLINE, pFile);
    if ( line[0] != '\n' )
    {
        PrintTextNearWrongFormat(ftell(pFile));
        //fprintf(stderr,"Wrong format at:\n %s\n", line);
        return 1;
    }

    step_destin.Masses.resize(step_destin.NumOfAtomTypes + 1, 0.0);
    for (int i = 0; i < step_destin.NumOfAtomTypes; i++)
    {
        count = fscanf(pFile, "%d %lf", &tmpint, &tmpdouble);
        if ( count != 2 )
        {
            fprintf(stderr, "Fail to read Masses section. Wrong format.\n");
            PrintTextNearWrongFormat(ftell(pFile));
            return 1;
        }
        fgets(line, MAXLINE, pFile);
        step_destin.Masses[tmpint] = tmpdouble;
    }
    return 0;
}

int ReadData::ReadAtoms(long &HeadPos, long &TailPos, STEP& step_destin, const char* atom_style)
{
    char* pch;
    char word[MAXLINE];
    int filestyleflag = 0;        // = 0 if not atom style is specified in data file. = 1 if atom style is specified in data file
    int userstyleflag = 1;        // = 0 if not atom style is specified by user. = 1 if atom style is specified by user
    int col_count = 0;

    int buffer_size;
//    char *pStart;
    const char* pChar;

    // clear the possible error state of the file stream and set the position
    clearerr(pFile);
    fseek(pFile, HeadPos, SEEK_SET);

    // decide atom style
    // check an appended comment for atom style
    fgets(line, MAXLINE, pFile);
    pch = strchr(line,'#');
    if ( pch != NULL )
    {
        pch++;
        pch = strtok(pch, " \t\r\n");
        if ( pch != NULL && strtok(NULL, " \t\r\n") == NULL )
            if ( strcmp(pch,"atomic") == 0 || strcmp(pch,"molecular") == 0 || strcmp(pch,"sphere") == 0 || strcmp(pch,"full") == 0 )
            {
                filestyleflag = 1;
                strcpy(word, pch);
            }
    }

    if ( strlen(atom_style) == 0 )
        userstyleflag = 0;

    if ( userstyleflag == 1 )
    {
        if ( filestyleflag == 1 && strcmp(word, atom_style) != 0 )
            fprintf(stderr,"WARNNING: LAMMPS atom style required may not match the file format.\n");
        strcpy(word, atom_style);
    }

    if ( userstyleflag == 0 && filestyleflag == 0 )
    {
        fprintf(stderr,"Cannot decide the atom style.\n");
        return 1;
    }
    else
    {
        if ( strcmp(word,"atomic") == 0 )
            style = ATOMIC;
        else if ( strcmp(word,"molecular") == 0 )
            style = MOLECULAR;
        else if ( strcmp(word,"sphere") == 0 )
            style = SPHERE;
        else if ( strcmp(word,"full") == 0 )
            style = FULL;
        else
        {
            if ( userstyleflag == 1 )
                fprintf(stderr,"Invalid atom style.\n");
            else
                fprintf(stderr,"Cannot decide the atom style.\n");
            return 1;
        }
    }

    // one whitespace line is required after the first line of one section
    fgets(line, MAXLINE, pFile);
    if ( line[0] != '\n' )
    {
        PrintTextNearWrongFormat(ftell(pFile));
        //fprintf(stderr,"Wrong format at:\n %s\n", line);
        return 1;
    }

    // count the number of column
    fgets(line, MAXLINE, pFile);
    pch = strtok(line, " \t\r\n");
    col_count++;
    while ( pch != NULL )
    {
        pch = strtok(NULL, " \t\r\n");
        col_count++;
    }
    col_count--;

    if ( col_count == NumCol[style] )
        imageflag = 0;
    else if ( col_count == NumCol[style] + 3 )
        imageflag = 1;
    else
    {
        fprintf(stderr,"Number of column does not match atom style.\n");
        return 1;
    }

    // assign coordinate type
    step_destin.wrapped_flag = 1;
    step_destin.scaled_flag = 0;

    // assign properties and property mapping
    vector<string>().swap(step_destin.Properties);
    if ( style == ATOMIC )
    {
        step_destin.Properties.resize(NumCol[style] + imageflag * 3 + velocityflag * NumVel[style]);
        step_destin.Properties[0] = "id"; step_destin.ID = 0;
        step_destin.Properties[1] = "type"; step_destin.TYPE = 1;
        step_destin.Properties[2] = "x"; step_destin.X = 2;
        step_destin.Properties[3] = "y"; step_destin.Y = 3;
        step_destin.Properties[4] = "z"; step_destin.Z = 4;
        if ( imageflag == 1 )
        {
            step_destin.Properties[5] = "ix"; step_destin.IX = 5;
            step_destin.Properties[6] = "iy"; step_destin.IY = 6;
            step_destin.Properties[7] = "iz"; step_destin.IZ = 7;
        }
        if ( velocityflag == 1 )
        {
            step_destin.Properties[5+imageflag*3] = "vx"; step_destin.VX = 5+imageflag*3;
            step_destin.Properties[6+imageflag*3] = "vy"; step_destin.VY = 6+imageflag*3;
            step_destin.Properties[7+imageflag*3] = "vz"; step_destin.VZ = 7+imageflag*3;
        }
    }
    else if ( style == MOLECULAR )
    {
        step_destin.Properties.resize(NumCol[style] + imageflag * 3 + velocityflag * NumVel[style]);
        step_destin.Properties[0] = "id"; step_destin.ID = 0;
        step_destin.Properties[1] = "mol"; step_destin.MOL = 1;
        step_destin.Properties[2] = "type"; step_destin.TYPE = 2;
        step_destin.Properties[3] = "x"; step_destin.X = 3;
        step_destin.Properties[4] = "y"; step_destin.Y = 4;
        step_destin.Properties[5] = "z"; step_destin.Z = 5;
        if ( imageflag == 1 )
        {
            step_destin.Properties[6] = "ix"; step_destin.IX = 6;
            step_destin.Properties[7] = "iy"; step_destin.IY = 7;
            step_destin.Properties[8] = "iz"; step_destin.IZ = 8;
        }
        if ( velocityflag == 1 )
        {
            step_destin.Properties[6+imageflag*3] = "vx"; step_destin.VX = 6+imageflag*3;
            step_destin.Properties[7+imageflag*3] = "vy"; step_destin.VY = 7+imageflag*3;
            step_destin.Properties[8+imageflag*3] = "vz"; step_destin.VZ = 8+imageflag*3;
        }
    }
    else if ( style == SPHERE )
    {
        step_destin.Properties.resize(NumCol[style] + imageflag * 3 + velocityflag * NumVel[style]);
        step_destin.Properties[0] = "id"; step_destin.ID = 0;
        step_destin.Properties[1] = "type"; step_destin.TYPE = 1;
        step_destin.Properties[2] = "diameter"; step_destin.DIAMETER = 2;
        step_destin.Properties[3] = "density"; step_destin.DENSITY = 3;
        step_destin.Properties[4] = "x"; step_destin.X = 4;
        step_destin.Properties[5] = "y"; step_destin.Y = 5;
        step_destin.Properties[6] = "z"; step_destin.Z = 6;
        if ( imageflag == 1 )
        {
            step_destin.Properties[7] = "ix"; step_destin.IX = 7;
            step_destin.Properties[8] = "iy"; step_destin.IY = 8;
            step_destin.Properties[9] = "iz"; step_destin.IZ = 9;
        }
        if ( velocityflag == 1 )
        {
            step_destin.Properties[7+imageflag*3] = "vx"; step_destin.VX = 7+imageflag*3;
            step_destin.Properties[8+imageflag*3] = "vy"; step_destin.VY = 8+imageflag*3;
            step_destin.Properties[9+imageflag*3] = "vz"; step_destin.VZ = 9+imageflag*3;
            step_destin.Properties[10+imageflag*3] = "wx"; step_destin.WX = 10+imageflag*3;
            step_destin.Properties[11+imageflag*3] = "wy"; step_destin.WY = 11+imageflag*3;
            step_destin.Properties[12+imageflag*3] = "wz"; step_destin.WZ = 12+imageflag*3;
        }
    }
    else if ( style == FULL )
    {
        step_destin.Properties.resize(NumCol[style] + imageflag * 3 + velocityflag * NumVel[style]);
        step_destin.Properties[0] = "id"; step_destin.ID = 0;
        step_destin.Properties[1] = "mol"; step_destin.MOL = 1;
        step_destin.Properties[2] = "type"; step_destin.TYPE = 2;
        step_destin.Properties[3] = "q"; step_destin.Q = 3;
        step_destin.Properties[4] = "x"; step_destin.X = 4;
        step_destin.Properties[5] = "y"; step_destin.Y = 5;
        step_destin.Properties[6] = "z"; step_destin.Z = 6;
        if ( imageflag == 1 )
        {
            step_destin.Properties[7] = "ix"; step_destin.IX = 7;
            step_destin.Properties[8] = "iy"; step_destin.IY = 8;
            step_destin.Properties[9] = "iz"; step_destin.IZ = 9;
        }
        if ( velocityflag == 1 )
        {
            step_destin.Properties[7+imageflag*3] = "vx"; step_destin.VX = 7+imageflag*3;
            step_destin.Properties[8+imageflag*3] = "vy"; step_destin.VY = 8+imageflag*3;
            step_destin.Properties[9+imageflag*3] = "vz"; step_destin.VZ = 9+imageflag*3;
        }
    }

    // read block of data
    clearerr(pFile);
    fseek(pFile, HeadPos, SEEK_SET);
    fgets(line, MAXLINE, pFile);
    fgets(line, MAXLINE, pFile);
    buffer_size = TailPos - ftell(pFile);
    block = new char [buffer_size];
    if ( (int)fread(block,1,buffer_size,pFile) != buffer_size )
    {
        fprintf(stderr, "Reading Error!\n");
        return 1;
    }

    // allocate space for step_destin.Atoms
    // allocation only happens when NumOfAtoms changes
    // Atom information will be written into step_destin.Atoms, step_destin.Atoms[i][j] means the jth property of the ith atom (i = atom index)
    if ( (int)step_destin.Atoms.size() != step_destin.NumOfAtoms || (int)step_destin.Atoms[0].size() != NumCol[style] + imageflag * 3 + velocityflag * 3)
    {
        vector< vector<double> >().swap(step_destin.Atoms);
        step_destin.Atoms.resize(step_destin.NumOfAtoms);
        for (int i = 0; i < step_destin.NumOfAtoms; i++)
            step_destin.Atoms[i].resize(NumCol[style] + imageflag * 3 + velocityflag * NumVel[style]);
    }

//    pStart = block;
    pChar = block;

    for (int i = 0; i < step_destin.NumOfAtoms; i++)
        for (int j = 0; j < NumCol[style] + imageflag * 3; j++)
        {
//            step_destin.Atoms[i][j] = strtod(pStart, &pStart);
            step_destin.Atoms[i][j] = StringToDouble(&pChar);
        }

    delete[] block;

    return 0;
}

int ReadData::ReadVelocities(long &HeadPos, long &TailPos, STEP& step_destin)
{
    int buffer_size;
//    char *pStart;
    const char* pChar;

    // build the atom id mapping
    vector<int> idmapping;
    if ( AtomIDMapping(idmapping, step_destin) )
        return 1;

    // read block of data
    clearerr(pFile);
    fseek(pFile, HeadPos, SEEK_SET);
    fgets(line, MAXLINE, pFile);

    // one whitespace line is required after the first line of one section
    fgets(line, MAXLINE, pFile);
    if ( line[0] != '\n' )
    {
        PrintTextNearWrongFormat(ftell(pFile));
        //fprintf(stderr,"Wrong format at:\n %s\n", line);
        return 1;
    }

    buffer_size = TailPos - ftell(pFile);
    block = new char [buffer_size];
    if ( (int)fread(block,1,buffer_size,pFile) != buffer_size )
    {
        fprintf(stderr, "Reading Error!\n");
        return 1;
    }

//    pStart = block;
    pChar = block;

    for (int i = 0; i < step_destin.NumOfAtoms; i++)
    {
        int tmpAtomID = StringToDouble(&pChar);
        for (int j = 0; j < NumVel[style]; j++)
        {
//            step_destin.Atoms[ idmapping[tmpAtomID] ][j + NumCol[style] + imageflag * 3] = strtod(pStart, &pStart);
            step_destin.Atoms[ idmapping[tmpAtomID] ][j + NumCol[style] + imageflag * 3] = StringToDouble(&pChar);
        }
    }
    delete[] block;

    return 0;
}

int ReadData::ReadBonds(long &HeadPos, long &TailPos, STEP& step_destin)
{
    int buffer_size;
//    char *pStart;
    const char* pChar;

    // read block of data
    clearerr(pFile);
    fseek(pFile, HeadPos, SEEK_SET);
    fgets(line, MAXLINE, pFile);

    // one whitespace line is required after the first line of one section
    fgets(line, MAXLINE, pFile);
    if ( line[0] != '\n' )
    {
        PrintTextNearWrongFormat(ftell(pFile));
        //fprintf(stderr,"Wrong format at:\n %s\n", line);
        return 1;
    }

    buffer_size = TailPos - ftell(pFile);
    block = new char [buffer_size];
    if ( (int)fread(block,1,buffer_size,pFile) != buffer_size )
    {
        fprintf(stderr, "Reading Error!\n");
        return 1;
    }

    // allocate space for step_destin.Bonds
    // allocation only happens when NumOfBonds changes
    // Bond information will be written into step_destin.Bonds, step_destin.Bonds[i][j] means properties the ith bond (i = bond index), j=0 for Bond ID, j=1 for Bond Type. j=2&3 for Atom ID of bonded atom
    if ( (int)step_destin.Bonds.size() != step_destin.NumOfBonds )
    {
        vector< vector<int> >().swap(step_destin.Bonds);
        step_destin.Bonds.resize(step_destin.NumOfBonds);
        for (int i = 0; i < step_destin.NumOfBonds; i++)
            step_destin.Bonds[i].resize(4);
    }

//    pStart = block;
    pChar = block;

    for (int i = 0; i < step_destin.NumOfBonds; i++)
        for (int j = 0; j < 4; j++)
        {
//            step_destin.Bonds[i][j] = strtol(pStart, &pStart);
            step_destin.Bonds[i][j] = StringToInt(&pChar);
        }

    delete[] block;

    return 0;
}

int ReadData::ReadAngles(long &HeadPos, long &TailPos, STEP& step_destin)
{
    int buffer_size;
//    char *pStart;
    const char* pChar;

    // read block of data
    clearerr(pFile);
    fseek(pFile, HeadPos, SEEK_SET);
    fgets(line, MAXLINE, pFile);

    // one whitespace line is required after the first line of one section
    fgets(line, MAXLINE, pFile);
    if ( line[0] != '\n' )
    {
        PrintTextNearWrongFormat(ftell(pFile));
        //fprintf(stderr,"Wrong format at:\n %s\n", line);
        return 1;
    }

    buffer_size = TailPos - ftell(pFile);
    block = new char [buffer_size];
    if ( (int)fread(block,1,buffer_size,pFile) != buffer_size )
    {
        fprintf(stderr, "Reading Error!\n");
        return 1;
    }

    // allocate space for step_destin.Angles
    // allocation only happens when NumOfAngles changes
    // Angle information will be written into step_destin.Angles, step_destin.Angles[i][j] means properties the ith angle (i = angle index), j=0 for AngleID, j=1 for Angle Type. j=2,3&4 for Atom ID
    if ( (int)step_destin.Angles.size() != step_destin.NumOfAngles )
    {
        vector< vector<int> >().swap(step_destin.Angles);
        step_destin.Angles.resize(step_destin.NumOfAngles);
        for (int i = 0; i < step_destin.NumOfAngles; i++)
            step_destin.Angles[i].resize(5);
    }

//    pStart = block;
    pChar = block;

    for (int i = 0; i < step_destin.NumOfAngles; i++)
        for (int j = 0; j < 5; j++)
        {
//            step_destin.Angles[i][j] = strtol(pStart, &pStart);
            step_destin.Angles[i][j] = StringToInt(&pChar);
        }

    delete[] block;

    return 0;
}

int ReadData::ReadDihedrals(long &HeadPos, long &TailPos, STEP& step_destin)
{
    int buffer_size;
//    char *pStart;
    const char* pChar;

    // read block of data
    clearerr(pFile);
    fseek(pFile, HeadPos, SEEK_SET);
    fgets(line, MAXLINE, pFile);

    // one whitespace line is required after the first line of one section
    fgets(line, MAXLINE, pFile);
    if ( line[0] != '\n' )
    {
        PrintTextNearWrongFormat(ftell(pFile));
        //fprintf(stderr,"Wrong format at:\n %s\n", line);
        return 1;
    }

    buffer_size = TailPos - ftell(pFile);
    block = new char [buffer_size];
    if ( (int)fread(block,1,buffer_size,pFile) != buffer_size )
    {
        fprintf(stderr, "Reading Error!\n");
        return 1;
    }

    // allocate space for step_destin.Dihedrals
    // allocation only happens when NumOfDihedrals changes
    // Dihedral information will be written into step_destin.Dihedrals, step_destin.Dihedrals[i][j] means properties the ith dihedral (i = dihedral index), j=0 for DihedralID, j=1 for Dihedral Type. j=2,3,4&5 for Atom ID
    if ( (int)step_destin.Dihedrals.size() != step_destin.NumOfDihedrals )
    {
        vector< vector<int> >().swap(step_destin.Dihedrals);
        step_destin.Dihedrals.resize(step_destin.NumOfDihedrals);
        for (int i = 0; i < step_destin.NumOfDihedrals; i++)
            step_destin.Dihedrals[i].resize(6);
    }

//    pStart = block;
    pChar = block;

    for (int i = 0; i < step_destin.NumOfDihedrals; i++)
        for (int j = 0; j < 6; j++)
        {
//            step_destin.Dihedrals[i][j] = strtol(pStart, &pStart);
            step_destin.Dihedrals[i][j] = StringToInt(&pChar);
        }

    delete[] block;

    return 0;
}

int ReadData::ReadImpropers(long &HeadPos, long &TailPos, STEP& step_destin)
{
    int buffer_size;
//    char *pStart;
    const char* pChar; 

    // read block of data
    clearerr(pFile);
    fseek(pFile, HeadPos, SEEK_SET);
    fgets(line, MAXLINE, pFile);

    // one whitespace line is required after the first line of one section
    fgets(line, MAXLINE, pFile);
    if ( line[0] != '\n' )
    {
        PrintTextNearWrongFormat(ftell(pFile));
        //fprintf(stderr,"Wrong format at:\n %s\n", line);
        return 1;
    }

    buffer_size = TailPos - ftell(pFile);
    block = new char [buffer_size];
    if ( (int)fread(block,1,buffer_size,pFile) != buffer_size )
    {
        fprintf(stderr, "Reading Error!\n");
        return 1;
    }

    // allocate space for step_destin.Impropers
    // allocation only happens when NumOfImpropers changes
    // Improper information will be written into step_destin.Impropers, step_destin.Impropers[i][j] means properties the ith improper (i = improper index), j=0 for ImproperID, j=1 for Improper Type. j=2,3,4&5 for Atom ID
    if ( (int)step_destin.Impropers.size() != step_destin.NumOfImpropers )
    {
        vector< vector<int> >().swap(step_destin.Impropers);
        step_destin.Impropers.resize(step_destin.NumOfImpropers);
        for (int i = 0; i < step_destin.NumOfImpropers; i++)
            step_destin.Impropers[i].resize(6);
    }

//    pStart = block;
    pChar = block;

    for (int i = 0; i < step_destin.NumOfImpropers; i++)
        for (int j = 0; j < 6; j++)
        {
//            step_destin.Impropers[i][j] = strtol(pStart, &pStart);
            step_destin.Impropers[i][j] = StringToInt(&pChar);
        }

    delete[] block;

    return 0;
}


void ReadData::PrintTextNearWrongFormat(long FilePos)
{
    int nline_before = 2;
    int nline_after = 2;
    int count = -1;
    int ch;
    fseek(pFile, FilePos, SEEK_SET);
    do
    {
        ch = fgetc(pFile);
        if ( ch == '\n' )
            count++;
        fseek(pFile, -2, SEEK_CUR);
    } while ( ch != EOF && count <= nline_before );
    fseek(pFile, 2, SEEK_CUR);

    fprintf(stderr,"Wrong format at:\n");
    for (int i = 0 - nline_before; i <= nline_after; i++)
    {
        fgets(line, MAXLINE, pFile);
        if ( i == 0 )
            fprintf(stderr,"> %s",line);
        else
            fprintf(stderr,"  %s",line);
    }
}

double ReadData::StringToDouble(const char **p)
{
    double r = 0.0;
    bool neg = false;

    while ( !(**p >= '0' && **p <= '9') && **p != '-' )
        ++(*p);

    if ( **p == '-' )
    {
        neg = true;
        ++(*p);
    }

    while ( **p >= '0' && **p <= '9' )
    {
        r = (r*10.0) + (**p - '0');
        ++(*p);
    }

    if ( **p == '.' )
    {
        double f = 0.0;
        double n = 1.0;
        ++(*p); 
        while ( **p >= '0' && **p <= '9' )
        {
            n = n / 10.0;
            f = f + (**p - '0') * n;
            ++(*p);
        } 
        r += f;
    }

    if ( **p == 'e' || **p == 'E' )
    {
        bool e_neg = false;
        double e = 0;

        ++(*p);
        if ( **p == '-' )
        {
            e_neg = true;
            ++(*p);
        }
        {
            ++(*p);
        }
        while ( **p >= '0' && **p <= '9' )
        {
            e = (e*10.0) + (**p - '0');
            ++(*p);
        }
        if ( e_neg )
            e = -e;

        r = r * pow(10.0,e);
    }

    if ( neg )
        r = -r;

    return r;
}

int ReadData::StringToInt(const char **p)
{
    int r = 0;
    bool neg = false;

    while ( !(**p >= '0' && **p <= '9') && **p != '-' )
        ++(*p);

    if ( **p == '-' )
    {
        neg = true;
        ++(*p);
    }

    while ( **p >= '0' && **p <= '9' )
    {
        r = (r*10) + (**p - '0');
        ++(*p);
    }

    if ( neg )
        r = -r;

    return r;
}


/** public function that is used to read file **/
int ReadData::Read_File(string FileName, STEP& step_destin, const char* atom_style)
{
    // test whether the file can be opened
    pFile = fopen(FileName.c_str(),"r");
    if ( pFile == NULL )
    {
        fprintf(stderr, "Can not open the file!\n");
        return 1;
    }
    clearerr(pFile);
    fseek(pFile, 0, SEEK_SET);

    if ( Init() )
        return 1;
    
    if ( ReadHeader(FilePosition, step_destin) )
        return 1;
    if ( FindSections(FilePosition, step_destin) )
        return 1;
    if ( massflag == 1 )
        if ( ReadMasses(massheadpos, masstailpos, step_destin) )
            return 1;
    if ( atomflag == 1 )
        if ( ReadAtoms(atomheadpos, atomtailpos, step_destin, atom_style) )
            return 1;
    if ( velocityflag == 1 )
        if ( ReadVelocities(velocityheadpos, velocitytailpos, step_destin) )
            return 1;
    if ( bondflag == 1 )
        if ( ReadBonds(bondheadpos, bondtailpos, step_destin) )
            return 1;
    if ( angleflag == 1 )
        if ( ReadAngles(angleheadpos, angletailpos, step_destin) )
            return 1;
    if ( dihedralflag == 1 )
        if ( ReadDihedrals(dihedralheadpos, dihedraltailpos, step_destin) )
            return 1;
    if ( improperflag == 1 )
        if ( ReadImpropers(improperheadpos, impropertailpos, step_destin) )
            return 1;
    fclose(pFile);
    return 0;
}

/** public functions that are used to get desired information **/
int ReadData::PropertyMapping(int& Index_destin, const string& PropertyName, const STEP& cur_step)
{
    Index_destin = -1; 
    for(int i = 0; i < (int)cur_step.Properties.size(); i++)
    {
        if ( cur_step.Properties[i] == PropertyName )
            Index_destin = i;
    }

    if ( Index_destin == -1 )
    {
        fprintf(stderr, "Required property \"%s\" is not found.\n", PropertyName.c_str());
        return 1;
    }

    return 0;
}

int ReadData::AtomIDMapping(vector<int>& Vector_destin, const STEP& cur_step)
{
    int ID;
    if ( cur_step.ID == -1 )
    {
        fprintf(stderr, "Property \"id\" is required, but not found.\n");
        return 1;
    }
    else
        ID = cur_step.ID;

    vector<int>().swap(Vector_destin);
    int maxID_index = 0;
    for (int i = 0; i < cur_step.NumOfAtoms; i++)
        if ( cur_step.Atoms[i][ID] > cur_step.Atoms[maxID_index][ID] )
            maxID_index = i;
    Vector_destin.resize(cur_step.Atoms[maxID_index][ID]+1,-1);

    for (int i = 0; i < cur_step.NumOfAtoms; i++)
        Vector_destin[cur_step.Atoms[i][ID]] = i;

    return 0;
}

int ReadData::BondIDMapping(vector<int>& Vector_destin, const STEP& cur_step)
{
    if ( bondflag == 0 )
    {
        fprintf(stderr, "No bond is read.\n");
        return 1;
    }

    vector<int>().swap(Vector_destin);
    int maxID = 0;
    for (int i = 0; i < cur_step.NumOfBonds; i++)
        if ( cur_step.Bonds[i][0] > maxID )
            maxID = cur_step.Bonds[i][0];
    Vector_destin.resize(maxID + 1,-1);

    for (int i = 0; i < cur_step.NumOfBonds; i++)
        Vector_destin[cur_step.Bonds[i][0]] = i;

    return 0;
}

/** public tool function **/
int ReadData::CalcWrappedCoord(STEP& cur_step)
{
    return cur_step.wrap();
}

int ReadData::CalcUnwrappedCoord(STEP& cur_step)
{
    return cur_step.unwrap();
}
