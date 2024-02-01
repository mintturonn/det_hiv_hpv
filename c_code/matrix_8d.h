#ifndef MATRIX8D_H
#define MATRIX8D_H

#include <stdio.h>

typedef double***** DMatrix5D;
typedef double**** DMatrix4D;
typedef double*** DMatrix3D;
typedef double** DMatrix;
typedef double* DArray;

// 1D

#define CREATE_DARRAY_START(_DS, _NAME) \
  new double[_DS + 1]; \
  for (int d1 = 1; d1 <= _DS; d1++) {

#define CREATE_DARRAY_END }

#define CREATE_DARRAY(_DS, _NAME, _VAL) \
  CREATE_DARRAY_START(_DS, _NAME) \
   _NAME[d1] = _VAL; \
  CREATE_DARRAY_END

#define FREE_DARRAY(_NAME) delete [] _NAME;

#define CLEAR_DARRAY(_NAME, _DS) for (int d1 = 1; d1 <= _DS; d1++) _NAME[_DS] = 0.0;

#define COPY_DARRAY(_FROM, _TO, _DS) \
  CREATE_DARRAY_START(_DS, _TO) \
    _TO[d1] = _FROM[d1]; \
  CREATE_DARRAY_END

void writeIntArray(FILE* f, int* a, int length);
void writeArray(FILE * f, DArray a, int length);

// 2D

#define CREATE_DMATRIX_START(_DS, _NAME) \
  new DArray[_DS[1] + 1]; \
  for (int d1 = 1; d1 <= _DS[1]; d1++) { \
      _NAME[d1] = new double[_DS[2] + 1]; \
      for (int d2 = 1; d2 <= _DS[2]; d2++) {

#define CREATE_DMATRIX_END }}

#define CREATE_DMATRIX(_DS, _NAME, _VAL) \
  CREATE_DMATRIX_START(_DS, _NAME) \
                _NAME[d1][d2] = _VAL; \
  CREATE_DMATRIX_END

void writeMatrix(FILE* f, DMatrix m, int* dims);

// 3D

#define CREATE_DMATRIX3D_START(_DS, _NAME) \
  new DMatrix[_DS[1] + 1]; \
  for (int d1 = 1; d1 <= _DS[1]; d1++) { \
    _NAME[d1] = new DArray[_DS[2] + 1]; \
    for (int d2 = 1; d2 <= _DS[2]; d2++) { \
      _NAME[d1][d2] = new double[_DS[3] + 1]; \
      for (int d3 = 1; d3 <= _DS[3]; d3++) {

#define CREATE_DMATRIX3D_END }}}

#define CREATE_DMATRIX3D(_DS, _NAME, _VAL) \
  CREATE_DMATRIX3D_START(_DS, _NAME) \
                _NAME[d1][d2][d3] = _VAL; \
  CREATE_DMATRIX3D_END

#define LOOP_DMATRIX3D_START(_DS) \
  for (int d1 = 1; d1 <= _DS[1]; d1++) { \
    for (int d2 = 1; d2 <= _DS[2]; d2++) { \
      for (int d3 = 1; d3 <= _DS[3]; d3++) {

#define LOOP_DMATRIX3D_END }}}

#define APPLY_3D(_VAR, _OP, _VAL, _F1, _T1, _F2, _T2, _F3, _T3) \
  for (int d1 =_F1; d1 <=_T1; d1++) \
    for (int d2 = _F2; d2 <= _T2; d2++) \
      for (int d3 = _F3; d3 <= _T3; d3++) \
        _VAR[d1][d2][d3] _OP _VAL;

#define COPY_3D(_FROM, _TO, _DS) \
  LOOP_DMATRIX3D_START(_DS) \
    _TO[d1][d2][d3] = _FROM[d1][d2][d3]; \
  LOOP_DMATRIX3D_END

void writeMatrix3D(FILE* f, DMatrix3D m, int* dims);

// 4D

#define CREATE_DMATRIX4D_START(_DS, _NAME) \
  new DMatrix3D[_DS[1] + 1]; \
  for (int d1 = 1; d1 <= _DS[1]; d1++) { \
    _NAME[d1] = new DMatrix[_DS[2] + 1]; \
    for (int d2 = 1; d2 <= _DS[2]; d2++) { \
      _NAME[d1][d2] = new DArray[_DS[3] + 1]; \
      for (int d3 = 1; d3 <= _DS[3]; d3++) { \
        _NAME[d1][d2][d3] = new double[_DS[4] + 1]; \
        for (int d4 = 1; d4 <= _DS[4]; d4++) {

#define CREATE_DMATRIX4D_END }}}}

#define CREATE_DMATRIX4D(_DS, _NAME, _VAL) \
  CREATE_DMATRIX4D_START(_DS, _NAME) \
                _NAME[d1][d2][d3][d4] = _VAL; \
  CREATE_DMATRIX4D_END

#define APPLY_4D(_VAR, _OP, _VAL, _F1, _T1, _F2, _T2, _F3, _T3, _F4, _T4) \
  for (int d1 =_F1; d1 <=_T1; d1++) \
    for (int d2 = _F2; d2 <= _T2; d2++) \
      for (int d3 = _F3; d3 <= _T3; d3++) \
        for (int d4 = _F4; d4 <= _T4; d4++) \
          _VAR[d1][d2][d3][d4] _OP _VAL;

void writeMatrix4D(FILE * f, DMatrix4D m, int* dims);

// 5D

#define CREATE_DMATRIX5D_START(_DS, _NAME) \
  new DMatrix4D[_DS[1] + 1]; \
  for (int d1 = 1; d1 <= _DS[1]; d1++) { \
    _NAME[d1] = new DMatrix3D[_DS[2] + 1]; \
    for (int d2 = 1; d2 <= _DS[2]; d2++) { \
      _NAME[d1][d2] = new DMatrix[_DS[3] + 1]; \
      for (int d3 = 1; d3 <= _DS[3]; d3++) { \
        _NAME[d1][d2][d3] = new DArray[_DS[4] + 1]; \
        for (int d4 = 1; d4 <= _DS[4]; d4++) { \
          _NAME[d1][d2][d3][d4] = new double[_DS[5] + 1]; \
          for (int d5 = 1; d5 <= _DS[5]; d5++) {

#define CREATE_DMATRIX5D_END }}}}}

#define CREATE_DMATRIX5D(_DS, _NAME, _VAL) \
  CREATE_DMATRIX5D_START(_DS, _NAME) \
                _NAME[d1][d2][d3][d4][d5] = _VAL; \
  CREATE_DMATRIX5D_END

void writeMatrix5D(FILE* f, DMatrix5D m, int* dims);

#endif
