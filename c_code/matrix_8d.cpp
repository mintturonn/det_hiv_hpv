#include "matrix_8d.h"

void writeMatrix(FILE* f, DMatrix m, int* dims) {
  fwrite(&dims[1], 4, 1, f);
  fwrite(&dims[2], 4, 1, f);
  for (int i = 1; i <= dims[1]; i++)
    for (int j = 1; j <= dims[2]; j++)
      fwrite(&m[i][j], 8, 1, f);
}

void writeIntArray(FILE* f, int* a, int length) {
  fwrite(&length, 4, 1, f);
  for (int i=1; i<=length; i++) fwrite(&a[i], 4, 1, f);
}

void writeArray(FILE* f, DArray a, int length) {
  fwrite(&length, 4, 1, f);
  for (int i = 1; i <= length; i++) fwrite(&a[i], 8, 1, f);
}

void writeMatrix3D(FILE* f, DMatrix3D m, int* dims) {
  fwrite(&dims[1], 4, 1, f);
  fwrite(&dims[2], 4, 1, f);
  fwrite(&dims[3], 4, 1, f);
  for (int i = 1; i <= dims[1]; i++)
    for (int j = 1; j <= dims[2]; j++)
      for (int k = 1; k <= dims[3]; k++)
        fwrite(&m[i][j][k], 8, 1, f);
}

void writeMatrix4D(FILE* f, DMatrix4D m, int* dims) {
  fwrite(&dims[1], 4, 1, f);
  fwrite(&dims[2], 4, 1, f);
  fwrite(&dims[3], 4, 1, f);
  fwrite(&dims[4], 4, 1, f);
  for (int i = 1; i <= dims[1]; i++)
    for (int j = 1; j <= dims[2]; j++)
      for (int k = 1; k <= dims[3]; k++)
        for (int l = 1; l <= dims[4]; l++)
        fwrite(&m[i][j][k][l], 8, 1, f);

}

void writeMatrix5D(FILE* f, DMatrix5D m, int* dims) {
  fwrite(&dims[1], 4, 1, f);
  fwrite(&dims[2], 4, 1, f);
  fwrite(&dims[3], 4, 1, f);
  fwrite(&dims[4], 4, 1, f);
  fwrite(&dims[5], 4, 1, f);
  for (int i = 1; i <= dims[1]; i++)
    for (int j = 1; j <= dims[2]; j++)
      for (int k = 1; k <= dims[3]; k++)
        for (int l = 1; l <= dims[4]; l++)
          for (int n = 1; n <= dims[5]; n++)
            fwrite(&m[i][j][k][l][n], 8, 1, f);
}
