%module plato
%{
#define SWIG_FILE_WITH_INIT
#include "Container.h"
#include "ReadDump.h"
#include "WriteDump.h"
#include "ReadData.h"
#include "WriteData.h"
#include "Neigh.h"
#include "Topo.h"
#include "forswig.h"
%}

// typemapping of primitive data type
%apply int& OUTPUT { int& Index_destin }    // this will be applied to function PropertyMapping in class ReadData and ReadDump


// Include the STL, and make some extension
%include "std_vector.i"
%include "std_string.i"

namespace std{

%template(_VecInt) vector<int>;
%template(_VecVecInt) vector< vector<int> >;
%template(_VecDouble) vector<double>;
%template(_VecVecDouble) vector< vector<double> >;
%template(_VecString) vector<string>;

}


// Include numpy library and make some typemapping
%include "numpy.i"

%init %{
    import_array();
%}

%apply (int** ARGOUTVIEWM_ARRAY1, int* DIM1) {(int** Numpy_destin, int* dim1)}
%apply (int** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(int** Numpy_destin, int* dim1, int* dim2)}
%apply (double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** Numpy_destin, int* dim1)}
%apply (double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(double** Numpy_destin, int* dim1, int* dim2)}
%apply (int* IN_ARRAY1, int DIM1) {(int* cur_Numpy, int dim1)}
%apply (double* IN_ARRAY1, int DIM1) {(double* cur_Numpy, int dim1)}
%apply (int* IN_ARRAY2, int DIM1, int DIM2) {(int* cur_Numpy, int dim1, int dim2)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* cur_Numpy, int dim1, int dim2)}

// hide cpp code and write python code as access to them
// (1) from forswig.h
%rename(_ConvertToNumpy_VecInt) ConvertToNumpy(const vector<int>&, int**, int*);
%rename(_ConvertToNumpy_VecVecInt) ConvertToNumpy(const vector< vector<int> >&, int**, int*, int*);
%rename(_ConvertToNumpy_VecDouble) ConvertToNumpy(const vector<double>&, double**, int*);
%rename(_ConvertToNumpy_VecVecDouble) ConvertToNumpy(const vector< vector<double> >&, double**, int*, int*);
%rename(_ConvertToC_Int) ConvertToC(vector<int>&, int*, int);
%rename(_ConvertToC_Double) ConvertToC(vector<double>&, double*, int);
%rename(_ConvertToC_2dInt) ConvertToC(vector< vector<int> >&, int*, int, int);
%rename(_ConvertToC_2dDouble) ConvertToC(vector< vector<double> >&, double*, int, int);
%rename(_Sqrt_2dDouble) Sqrt_2dDouble(vector< vector<double> >&);

// (2) from Container.h
%copyctor STEP;        // generate an copy constructor by swig
%rename(_STEP) STEP;
%rename(_BOX) BOX;

// (3) from ReadDump.h
%rename(_ReadDump) ReadDump;
/*
%rename(_Read_Timestep) Read_Timestep;
%rename(_Read_Frame) Read_Frame;
%rename(_GetTimesteps) GetTimesteps;
%rename(_GetProperties) GetProperties;
%rename(_PropertyMapping) PropertyMapping;
%rename(_IDMapping) IDMapping;
%rename(_CalcWrappedCoord) CalcWrappedCoord;
%rename(_CalcUnwrappedCoord) CalcUnwrappedCoord;
*/

// (4) from WriteDump.h
%rename(_WriteDump) WriteDump;
/*
%rename(_Write) Write;
*/

// (5) from ReadData.h
%rename(_ReadData) ReadData;
/*
%rename(_Read_File) Read_File;
%rename(_PropertyMapping) PropertyMapping;
%rename(_IDMapping) IDMapping;
%rename(_CalcWrappedCoord) CalcWrappedCoord;
%rename(_CalcUnwrappedCoord) CalcUnwrappedCoord;
*/

// (6) from WriteData.h
%rename(_WriteData) WriteData;
/*
%rename(_Write) Write;
*/

// (7) from Neigh.h
%rename(_Neigh) Neigh;
/*
%rename(_Prepare) Prepare;
%rename(_Query) Query;
%rename(_BuildNeighborList) BuildNeighborList;
%reneame(_Cluster) Cluster;
%reneame(_RadialDistributionFunction) RadialDistributionFunction;
*/

// (8) from Topo.h
%rename(_Topo) Topo;
/*
%rename (_ClusterAnalysis) ClusterAnalysis;
%rename (_Connectivity) Connectivity;
*/

// Include the python codes
%pythoncode "dump.py"
%pythoncode "data.py"
%pythoncode "neigh.py"
%pythoncode "topo.py"
%pythoncode "hbond.py"

/*%feature("autodoc", "2");*/

// Include the header file with above prototypes
%include "Container.h"
%include "ReadDump.h"
%include "WriteDump.h"
%include "ReadData.h"
%include "WriteData.h"
%include "Neigh.h"
%include "Topo.h"
%include "forswig.h"


%pythonbegin %{
import numpy as np
np.set_printoptions(suppress=True)

import gc
gc.enable()
%}
