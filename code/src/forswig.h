#ifndef FORSWIG
#define FORSWIG

#include <vector>
#include <string>
#include <map>
#include <cmath>

using namespace std;

void ConvertToNumpy(const vector<int>& cur_vec, int** Numpy_destin, int* dim1)
{
    *dim1 = (int)cur_vec.size();

    int* tmp;
    tmp = (int *)malloc((*dim1)*sizeof(int));
    for (int i = 0; i < (*dim1); i++)
        tmp[i] = cur_vec[i];

    *Numpy_destin = tmp;
}

void ConvertToNumpy(const vector<double>& cur_vec, double** Numpy_destin, int* dim1)
{
    *dim1 = (int)cur_vec.size();

    double* tmp;
    tmp = (double *)malloc((*dim1)*sizeof(double));
    for (int i = 0; i < (*dim1); i++)
        tmp[i] = cur_vec[i];

    *Numpy_destin = tmp;
}

void ConvertToNumpy(const vector< vector<int> >& cur_vec, int** Numpy_destin, int* dim1, int* dim2)
{
    int* tmp;
    *dim1 = (int)cur_vec.size();

    if ( (*dim1) == 0 )
    {
        *dim2 = 0;
        tmp = (int *)malloc(0);

        *Numpy_destin = tmp;
    }
    else
    {
        int matrix_flag = 1;
        int tmpdim2 = (int)cur_vec[0].size();
        for (int i = 1; i < (*dim1); i++)
        {
            if ( tmpdim2 == (int)cur_vec[i].size())
                tmpdim2 = (int)cur_vec[i].size();
            else
            {    matrix_flag = 0;
                 break;
            }
        }

        if ( matrix_flag == 0 )
        {
            fprintf(stderr, "Not a matrix, cannot convert.\n");
            *dim1 = 0;
            *dim2 = 0;
            tmp = (int *)malloc(0);

            *Numpy_destin = tmp;
        }
        else
        {
            *dim2 = tmpdim2;
            tmp = (int *)malloc((*dim1)*(*dim2)*sizeof(int));
            for (int i = 0; i < (*dim1); i++)
                for (int j = 0; j < (*dim2); j++)
                    tmp[i*(*dim2)+j] = cur_vec[i][j];

            *Numpy_destin = tmp;
        }
    }
}

void ConvertToNumpy(const vector< vector<double> >& cur_vec, double** Numpy_destin, int* dim1, int* dim2)
{
    double* tmp;
    *dim1 = (int)cur_vec.size();

    if ( (*dim1) == 0 )
    {
        *dim2 = 0;
        tmp = (double *)malloc(0);

        *Numpy_destin = tmp;
    }
    else
    {
        int matrix_flag = 1;
        int tmpdim2 = (int)cur_vec[0].size();
        for (int i = 1; i < (*dim1); i++)
        {
            if ( tmpdim2 == (int)cur_vec[i].size())
                tmpdim2 = (int)cur_vec[i].size();
            else
            {    matrix_flag = 0;
                 break;
            }
        }

        if ( matrix_flag == 0 )
        {
            fprintf(stderr, "Not a matrix, cannot convert.\n");
            *dim1 = 0;
            *dim2 = 0;
            tmp = (double *)malloc(0);

            *Numpy_destin = tmp;
        }
        else
        {
            *dim2 = tmpdim2;
            tmp = (double *)malloc((*dim1)*(*dim2)*sizeof(double));
            for (int i = 0; i < (*dim1); i++)
                for (int j = 0; j < (*dim2); j++)
                    tmp[i*(*dim2)+j] = cur_vec[i][j];

            *Numpy_destin = tmp;
        }
    }
}

void ConvertToC(vector<double>& Vector_destin, double* cur_Numpy, int dim1)
{
    vector<double>().swap(Vector_destin);
    Vector_destin.resize(dim1);
    for (int i = 0; i < dim1; i++)
        Vector_destin[i] = cur_Numpy[i];
}

void ConvertToC(vector<int>& Vector_destin, int* cur_Numpy, int dim1)
{
    vector<int>().swap(Vector_destin);
    Vector_destin.resize(dim1);
    for (int i = 0; i < dim1; i++)
        Vector_destin[i] = cur_Numpy[i];
}

void ConvertToC(vector< vector<double> >& VectorVector_destin, double* cur_Numpy, int dim1, int dim2)
{
    vector< vector<double> >().swap(VectorVector_destin);
    VectorVector_destin.resize(dim1);
    for (int i = 0; i < dim1; i++)
        VectorVector_destin[i].resize(dim2);
    for (int i = 0; i < dim1; i++)
        for (int j = 0; j < dim2; j++)
            VectorVector_destin[i][j] = cur_Numpy[i*dim2 + j];
}

void ConvertToC(vector< vector<int> >& VectorVector_destin, int* cur_Numpy, int dim1, int dim2)
{
    vector< vector<int> >().swap(VectorVector_destin);
    VectorVector_destin.resize(dim1);
    for (int i = 0; i < dim1; i++)
        VectorVector_destin[i].resize(dim2);
    for (int i = 0; i < dim1; i++)
        for (int j = 0; j < dim2; j++)
            VectorVector_destin[i][j] = cur_Numpy[i*dim2 + j];
}

void Sqrt_2dDouble(vector< vector<double> >& cur_vec)
{
    for (int i = 0; i < (int)cur_vec.size(); i++)
        for (int j = 0; j < (int)cur_vec[i].size(); j++)
            cur_vec[i][j] = sqrt(cur_vec[i][j]);
}
#endif
