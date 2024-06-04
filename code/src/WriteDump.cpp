/***********************************************************/
/*** Refer to header file for instructions of this class ***/
/***********************************************************/
/*** Version: 11/15/2019 (1.1 Version) ***/
/*** By: Heyi Liang                    ***/
/*****************************************/

#include "WriteDump.h"

WriteDump::WriteDump(string FileName, const vector<string>& Properties, const string& mode)
{
    // store the parameters
    OutFileName = FileName;
    vector<string>().swap(OutProperties);
    OutProperties = Properties;

    if ( mode.compare("w") != 0 && mode.compare("a") != 0 )
    {
        fprintf(stderr, "Mode should be either \"w\"(create/overwrite) or \"a\"(append).\n");
        exit(1);
    }

    FILE* pFile = fopen(OutFileName.c_str(), mode.c_str());
    if ( pFile == NULL )
    {
        fprintf(stderr, "Can not create/open the file!\n");
        exit(1);
    }
    fclose(pFile);

    // Check whether OutProperties is empty
    if ( (int)OutProperties.size() == 0 )
    {
        fprintf(stderr, "Please specify the properties to be written to file.\n");
        exit(1);
    }

}

WriteDump::~WriteDump()
{

}

void WriteDump::print_version()
{
    printf("Running WriteDump.cpp, version %s\n", WRITEDUMP_VERSION);
}

/** public function that is used to write file **/
int WriteDump::Write(const STEP& cur_step)
{
    FILE* pFile = fopen(OutFileName.c_str(), "a");
    if ( pFile == NULL )
    {
        fprintf(stderr, "Can not create the file!\n");
        exit(1);
    }
    //int MaxBufsiz = 409600;
    //int BufWarn = 400000;
    //char* sbuf = new char [MaxBufsiz];
    //int offset = 0;
    int NumOfProp = (int)OutProperties.size();
    Index_OutProperties.resize(NumOfProp, -1);
    for (int i = 0; i < NumOfProp; i++)
    {
        for (int j = 0; j < (int)cur_step.Properties.size(); j++)
              if ( OutProperties[i] == cur_step.Properties[j] )
                  Index_OutProperties[i] = j;
        if ( Index_OutProperties[i] == -1 )
        {
            fprintf(stderr, "Property \"%s\" not available.\n", OutProperties[i].c_str());
            return 1;
        }
    }

    fprintf(pFile, "ITEM: TIMESTEP\n%d\n", cur_step.Timestep);
    fprintf(pFile, "ITEM: NUMBER OF ATOMS\n%d\n", cur_step.NumOfAtoms);
    //offset += sprintf(&sbuf[offset], "ITEM: TIMESTEP\n%d\n", cur_step.Timestep);
    //offset += sprintf(&sbuf[offset], "ITEM: NUMBER OF ATOMS\n%d\n", cur_step.NumOfAtoms);

    if (cur_step.Box.triclinic) {
        fprintf(pFile, "ITEM: BOX BOUNDS xy xz yz %s %s %s\n", cur_step.Box.xboundary, cur_step.Box.yboundary, cur_step.Box.zboundary);
        fprintf(pFile, "%-1.16e %-1.16e %-1.16e\n", cur_step.Box.xlo_bound, cur_step.Box.xhi_bound, cur_step.Box.xy);
        fprintf(pFile, "%-1.16e %-1.16e %-1.16e\n", cur_step.Box.ylo_bound, cur_step.Box.yhi_bound, cur_step.Box.xz);
        fprintf(pFile, "%-1.16e %-1.16e %-1.16e\n", cur_step.Box.zlo_bound, cur_step.Box.zhi_bound, cur_step.Box.yz);
    } else {
        fprintf(pFile, "ITEM: BOX BOUNDS %s %s %s\n", cur_step.Box.xboundary, cur_step.Box.yboundary, cur_step.Box.zboundary);
        fprintf(pFile, "%-1.16e %-1.16e\n", cur_step.Box.xlo_bound, cur_step.Box.xhi_bound);
        fprintf(pFile, "%-1.16e %-1.16e\n", cur_step.Box.ylo_bound, cur_step.Box.yhi_bound);
        fprintf(pFile, "%-1.16e %-1.16e\n", cur_step.Box.zlo_bound, cur_step.Box.zhi_bound);
    }
    //offset += sprintf(&sbuf[offset], "ITEM: BOX BOUNDS %s %s %s\n", cur_step.Box.xboundary, cur_step.Box.yboundary, cur_step.Box.zboundary);
    //offset += sprintf(&sbuf[offset], "%g %g\n", cur_step.Box.xlo, cur_step.Box.xhi);
    //offset += sprintf(&sbuf[offset], "%g %g\n", cur_step.Box.ylo, cur_step.Box.yhi);
    //offset += sprintf(&sbuf[offset], "%g %g\n", cur_step.Box.zlo, cur_step.Box.zhi);

    fprintf(pFile, "ITEM: ATOMS ");
    for (int i = 0; i < NumOfProp; i++)
        fprintf(pFile, "%s ", OutProperties[i].c_str());
    fprintf(pFile, "\n");
    //offset += sprintf(&sbuf[offset], "ITEM: ATOMS ");
    //for (int i = 0; i < NumOfProp; i++)
    //    offset += sprintf(&sbuf[offset], "%s ", OutProperties[i].c_str());
    //offset += sprintf(&sbuf[offset], "\n");

    //fwrite(sbuf, offset, 1, pFile);

    //offset = 0;
    for (int i = 0; i < cur_step.NumOfAtoms; i++)
    {
        for (int j = 0; j < NumOfProp; j++)
        {
            if ( OutProperties[j] == "id" || OutProperties[j] == "mol" || OutProperties[j] == "type" )
                fprintf(pFile, "%d ", (int)cur_step.Atoms[i][Index_OutProperties[j]]);
                //offset += sprintf(&sbuf[offset], "%d ", (int)cur_step.Atoms[i][Index_OutProperties[j]]);
            else
                fprintf(pFile, "%g ", cur_step.Atoms[i][Index_OutProperties[j]]);
                //offset += sprintf(&sbuf[offset], "%g ", cur_step.Atoms[i][Index_OutProperties[j]]);
        }
        fprintf(pFile, "\n");
        //offset += sprintf(&sbuf[offset], "\n");
        //if ( offset > BufWarn )
        //{
        //    fwrite(sbuf, offset ,1, pFile);
        //    offset = 0;
        //}
    }

    fclose(pFile);
    return 0;
}

