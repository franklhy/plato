/***********************************************************/
/*** Refer to header file for instructions of this class ***/
/***********************************************************/
/*** Version: 07/16/2016 (1.0 Version)  ***/
/*** By: Heyi Liang                     ***/
/******************************************/

#include "ReadDump.h"

#define MAXLINE 1024

ReadDump::ReadDump(string FileName)
{
    ID=TYPE=MOL=X=Y=Z=IX=IY=IZ=VX=VY=VZ=WX=WY=WZ=FX=FY=FZ=Q=-1;

    line = new char[MAXLINE];

    // record the file name
    InFileName = FileName;

    // test whether the file can be opened
    pFile = fopen(InFileName.c_str(),"r");
    if ( pFile == NULL )
    {
        fprintf(stderr, "Can not open the file!\n");
        exit(1);
    }

    // Read timesteps. Refer to corresponding functions in header file
    if ( CheckAndScanFile() )
    {
        fprintf(stderr, "Error while running function CheckAndScanFile() in ReadDump.cpp!\n");
        exit(1);
    }

    // Get available atom information
    if ( FindProperties() )
    {
        fprintf(stderr, "Error while running function FindProperties() in ReadDump.cpp!\n");
        exit(1);
    }

    // column mapping
    if ( BasicPropertyMapping() )
    {
        fprintf(stderr, "Fail to do basic mapping!\n");
        exit(1);
    }
    clearerr(pFile);
    fseek(pFile, 0, SEEK_SET);
}

ReadDump::~ReadDump()
{
    delete[] line;
    fclose(pFile);
}

void ReadDump::print_version()
{
    printf("Running ReadDump.cpp, version %s\n", READDUMP_VERSION);
}

/** functions that are only called by constructor **/
int ReadDump::CheckAndScanFile()
{
    long fpos, flength;
    int tmpint;

    // get the size of file
    clearerr(pFile);
    fseek(pFile, 0, SEEK_END);
    flength = ftell(pFile);
    // clear the possible error state of the file stream and set the position to the beginning of the file
    clearerr(pFile);
    fseek(pFile, 0, SEEK_SET);

    // check file format and estimate total number of snapshots (i.e. timesteps)
    fgets(line, MAXLINE, pFile);
    if ( strstr(line, "ITEM: TIMESTEP") != line )        // first line of the file should be "ITEM: TIMESTEP"
    {
        fprintf(stderr, "Wrong file format!\n");
        return 1;
    }
    while ( fgets(line, MAXLINE, pFile) != NULL )
    {
        if ( strstr(line, "ITEM: TIMESTEP") == line )
            break;
    }
    tmpint = flength / ftell(pFile);
    vector<int>().swap(Timesteps);
    Timesteps.reserve(tmpint);
    TimestepsHeadPos.reserve(tmpint);
    // clear the possible error state of the file stream and set the position to the beginning of the file
    clearerr(pFile);
    fseek(pFile, 0, SEEK_SET);
    fpos = 0;

    // read the file until reach the end of file
    while ( fgets(line, MAXLINE, pFile) != NULL )
    {
        if ( strstr(line, "ITEM: TIMESTEP") == line )
        {
            // read next line for timestep, and record the head position for current snapshot
            fgets(line, MAXLINE, pFile);
            Timesteps.push_back(atoi(line));
            TimestepsHeadPos.push_back(fpos);
            printf("Scanning file... %d%%\r",(int)(ftell(pFile) * 100 / flength));

            // read next two line for number of atom, then skip all the data until the next snapshot is reached
            fgets(line, MAXLINE, pFile);
            if ( strstr(line, "ITEM: NUMBER OF ATOMS") == line )
            {
                fgets(line, MAXLINE, pFile);
                tmpint = atoi(line);
                for (int i = 0; i < tmpint + 5; i++)
                    fgets(line, MAXLINE, pFile);
            }
            else
            {
                fprintf(stderr, "Wrong file format!\n");
                return 1;
            }
            
        }
        fpos = ftell(pFile);
    }

    printf("\n");

    clearerr(pFile);
    return 0;
}

int ReadDump::FindProperties()
{
    char* pch;
    string tmps;

    // clear the possible error state of the file stream and set the position to the beginning of the file
    clearerr(pFile);
    fseek(pFile, 0, SEEK_SET);

    // read the file
    while ( fgets(line, MAXLINE, pFile) != NULL )
        // when read "ITEM: ATOMS" from file, get the remaining part of this line
        if ( strstr(line, "ITEM: ATOMS") == line )
            break;
    // tokenize line and skip first two token (which is "ITEM:" and "ATOMS")
    pch = strtok(line, " \t\n");
    pch = strtok(NULL, " \t\n");

    // clear the vector Properties 
    vector<string>().swap(Properties);

    // get available types of atom information	
    pch = strtok(NULL, " \t\n");
    while ( pch != NULL )
    {
        tmps = pch;

        // record the available types of atom information into vector Properties
        Properties.push_back(tmps);

        // if the coordinate information is in the file, get the type of coordinate, and assign it. Since three coordinate have the same type, only need to compare one of them
        if ( tmps == "x" )
        {
            wrapped_flag = 2;
            scaled_flag = 0;
        }
        else if ( tmps == "xs" )
        {
            wrapped_flag = 2;
            scaled_flag = 1;
        }
        else if ( tmps == "xu" )
        {
            wrapped_flag = 0;
            scaled_flag = 0;
        }
        else if ( tmps == "xsu" )
        {
            wrapped_flag = 0;
            scaled_flag = 1;
        }

        pch = strtok(NULL, " \t\n");
    }

    return 0;
}

int ReadDump::BasicPropertyMapping()
{
    for (int i = 0; i < (int)Properties.size(); i++)
    {
        if ( Properties[i] == "id" )
            ID = i;
        else if ( Properties[i] == "mol" )
            MOL = i;
        else if ( Properties[i] == "type" )
            TYPE = i;
        else if ( Properties[i] == "x" || Properties[i] == "xs" || Properties[i] == "xu" || Properties[i] == "xsu" )
            X = i;
        else if ( Properties[i] == "y" || Properties[i] == "ys" || Properties[i] == "yu" || Properties[i] == "ysu" )
            Y = i;
        else if ( Properties[i] == "z" || Properties[i] == "zs" || Properties[i] == "zu" || Properties[i] == "zsu" )
            Z = i;
        else if ( Properties[i] == "ix" )
            IX = i;
        else if ( Properties[i] == "iy" )
            IY = i;
        else if ( Properties[i] == "iz" )
            IZ = i;
        else if ( Properties[i] == "vx" )
            VX = i;
        else if ( Properties[i] == "vy" )
            VY = i;
        else if ( Properties[i] == "vz" )
            VZ = i;
        else if ( Properties[i] == "wx" )
            WX = i;
        else if ( Properties[i] == "wy" )
            WY = i;
        else if ( Properties[i] == "wz" )
            WZ = i;
        else if ( Properties[i] == "fx" )
            FX = i;
        else if ( Properties[i] == "fy" )
            FY = i;
        else if ( Properties[i] == "fz" )
            FZ = i;
        else if ( Properties[i] == "q" )
            Q = i;
    }
    return 0;
}

/** functions that are only called when reading the file, they all called by function Read **/
int ReadDump::ReadOneSnapshot(int Frame, long FilePos, STEP& step_destin)
{
    int tmpint;
    char *pch;
    int buffer_size;
    int NumOfProperties;
    //char *pStart;
    const char* pChar;

    // clear the possible error state of the file stream and set the position to the beginning of the file
    clearerr(pFile);
    fseek(pFile, FilePos, SEEK_SET);

    // assign properties
    vector<string>().swap(step_destin.Properties);
    step_destin.Properties = Properties;

    // assign coordinate type
    step_destin.wrapped_flag = wrapped_flag;
    step_destin.scaled_flag = scaled_flag;

    // assign basic properties mapping
    step_destin.ID = ID;
    step_destin.TYPE = TYPE;
    step_destin.MOL = MOL;
    step_destin.X = X;
    step_destin.Y = Y;
    step_destin.Z = Z;
    step_destin.IX = IX;
    step_destin.IY = IY;
    step_destin.IZ = IZ;
    step_destin.VX = VX;
    step_destin.VY = VY;
    step_destin.VZ = VZ;
    step_destin.WX = WX;
    step_destin.WY = WY;
    step_destin.WZ = WZ;
    step_destin.FX = FX;
    step_destin.FY = FY;
    step_destin.FZ = FZ;
    step_destin.Q = Q;

/* read timestep */
    // make sure it is reading the snapshot with required Timestep
    fgets(line, MAXLINE, pFile);
    if ( strstr(line, "ITEM: TIMESTEP") == line )
    {
        fgets(line, MAXLINE, pFile);
        tmpint = atoi(line);
        if ( tmpint != Timesteps[Frame] )
        {
            fprintf(stderr, "Cannot find snapshot with required Timestep = %d !\n", Timesteps[Frame]);
            return 1;
        }
    }
    else
    {
        fprintf(stderr, "Wrong format.\n");
        return 1;
    }
    step_destin.Timestep = Timesteps[Frame];

/* read number of atoms */
    fgets(line, MAXLINE, pFile);
    if ( strstr(line, "ITEM: NUMBER OF ATOMS") != line )
    {
        fprintf(stderr, "The input file is not in correct format!\n");
        return 1;
    }
    fgets(line, MAXLINE, pFile);
    step_destin.NumOfAtoms = atoi(line);

/* read box */
    fgets(line, MAXLINE, pFile);
    if ( strstr(line, "ITEM: BOX BOUNDS") != line )
    {
        fprintf(stderr, "The input file is not in correct format!\n");
        return 1;
    }

    // tokenize line and discard first three tokens ( i.e. "ITEM:" "BOX" and "BOUNDS")
    pch = strtok(line, " \t\n");
    pch = strtok(NULL, " \t\n");
    pch = strtok(NULL, " \t\n");

    // read the remaining information of this line word by word for boundary information
    pch = strtok(NULL, " \t\n");
    if (pch[0] == 'x' && pch[1] == 'y') {
        step_destin.Box.triclinic = true;
        // skip next two tokens (i.e. "xz" and "yz")
        pch = strtok(NULL, " \t\n");
        pch = strtok(NULL, " \t\n");
    } else {
        step_destin.Box.xboundary[0] = pch[0];
        step_destin.Box.xboundary[1] = pch[1];
        pch = strtok(NULL, " \t\n");
        step_destin.Box.yboundary[0] = pch[0];
        step_destin.Box.yboundary[1] = pch[1];
        pch = strtok(NULL, " \t\n");
        step_destin.Box.zboundary[0] = pch[0];
        step_destin.Box.zboundary[1] = pch[1];
    }

    // read next three lines for box size
    if (step_destin.Box.triclinic) {
        fgets(line, MAXLINE, pFile);
        step_destin.Box.xlo_bound = strtod(line, &pch);
        step_destin.Box.xhi_bound = strtod(pch, &pch);
        step_destin.Box.xy = strtod(pch, NULL);
        fgets(line, MAXLINE, pFile);
        step_destin.Box.ylo_bound = strtod(line, &pch);
        step_destin.Box.yhi_bound = strtod(pch, &pch);
        step_destin.Box.xz = strtod(pch, NULL);
        fgets(line, MAXLINE, pFile);
        step_destin.Box.zlo_bound = strtod(line, &pch);
        step_destin.Box.zhi_bound = strtod(pch, &pch);
        step_destin.Box.yz = strtod(pch, NULL);
        step_destin.Box.triclinic = true;
    } else {
        fgets(line, MAXLINE, pFile);
        step_destin.Box.xlo_bound = strtod(line, &pch);
        step_destin.Box.xhi_bound = strtod(pch, NULL);
        fgets(line, MAXLINE, pFile);
        step_destin.Box.ylo_bound = strtod(line, &pch);
        step_destin.Box.yhi_bound = strtod(pch, NULL);
        fgets(line, MAXLINE, pFile);
        step_destin.Box.zlo_bound = strtod(line, &pch);
        step_destin.Box.zhi_bound = strtod(pch, NULL);
        step_destin.Box.xy = 0;
        step_destin.Box.xz = 0;
        step_destin.Box.yz = 0;
        step_destin.Box.triclinic = false;
    }
    step_destin.Box.set_box(true);

/* read atoms */
    // read the next lines to insure that it is atoms information
    fgets(line, MAXLINE, pFile);
    if ( strstr(line, "ITEM: ATOMS") != line )
    {
        fprintf(stderr, "The input file is not in correct format!\n");
        return 1;
    }

    // read the whole snapshot into "block"
    long tmpFilePos = ftell(pFile);
    if ( Frame == (int)Timesteps.size() - 1 )
    {
        fseek(pFile,0,SEEK_END);
        buffer_size = ftell(pFile) - tmpFilePos;
        clearerr(pFile);
        fseek(pFile, tmpFilePos, SEEK_SET);
    }
    else
    {
        buffer_size = TimestepsHeadPos[Frame+1] - tmpFilePos;
    }
    block = new char [buffer_size];
    if ( (int)fread(block,1,buffer_size,pFile) != buffer_size )
    {
        fprintf(stderr, "Reading Error!\n");
        return 1;
    }

    // allocate space for step_destin.Atoms
    // allocation only happens when NumOfAtoms changes
    // Atom information will be written into step_destin.Atoms, step_destin.Atoms[i][j] means the jth property of the ith atom (i = atom index)
    NumOfProperties = (int)Properties.size();
    if ( (int)step_destin.Atoms.size() != step_destin.NumOfAtoms || (int)step_destin.Atoms[0].size() != NumOfProperties )
    {
        vector< vector<double> >().swap(step_destin.Atoms);
        step_destin.Atoms.resize(step_destin.NumOfAtoms);
        for (int i = 0; i < step_destin.NumOfAtoms; i++)
        step_destin.Atoms[i].resize(NumOfProperties);
    }

//    pStart = block;
    pChar = block;

    for (int i = 0; i < step_destin.NumOfAtoms; i++)
        for (int j = 0; j < NumOfProperties; j++)
        {
//            step_destin.Atoms[i][j] = strtod(pStart, &pStart);
            step_destin.Atoms[i][j] = StringToDouble(&pChar);
        }

    delete[] block;

    return 0;
}

double ReadDump::StringToDouble(const char **p)
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
        else if ( **p == '+' )
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

/** public function that is used to read file **/
int ReadDump::Read_Timestep(int Timestep, STEP& step_destin)
{
    long FilePos;
    int FindState, Frame;

    FindState = 0;
    for (int i = 0; i < (int)Timesteps.size(); i++)
    {
        if ( Timesteps[i] == Timestep )
        {
            FilePos = TimestepsHeadPos[i];
            FindState = 1;
            Frame = i;
            break;
        }
    }
    if ( FindState == 0 )
    {
        fprintf(stderr, "Cannot find Timestep = %d\n", Timestep);
        return 1;
    }

    // read one snapshot
    if ( ReadOneSnapshot(Frame,FilePos,step_destin) )
    {
        fprintf(stderr, "Error when reading Timestep = %d\n", Timestep);
        return 1;
    }

    return 0;
}

int ReadDump::Read_Frame(int Frame, STEP& step_destin)
{
    long FilePos;

    if ( Frame < 0 || Frame > (int)Timesteps.size() - 1 )
    {
        fprintf(stderr, "There is no Frame No. %d\n", Frame);
        return 1;
    }

    // read one snapshot
    FilePos = TimestepsHeadPos[Frame];
    if ( ReadOneSnapshot(Frame,FilePos,step_destin) )
    {
        fprintf(stderr, "Error when reading Frame No. %d\n", Frame);
        return 1;
    }

    return 0;

}

/** public functions that are used to get desired information **/
int ReadDump::GetTimesteps(vector<int>& Vector_destin)
{
    Vector_destin = Timesteps;
    return 0;
}

int ReadDump::GetProperties(vector<string>& Vector_destin)
{
    Vector_destin = Properties;
    return 0;
}

int ReadDump::PropertyMapping(int& Index_destin, const string& PropertyName, const STEP& cur_step)
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

int ReadDump::IDMapping(vector<int>& Vector_destin, const STEP& cur_step)
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

/** public tool function **/
int ReadDump::CalcWrappedCoord(STEP& cur_step)
{
    return cur_step.wrap();
}

int ReadDump::CalcUnwrappedCoord(STEP& cur_step)
{
    return cur_step.unwrap();
}
