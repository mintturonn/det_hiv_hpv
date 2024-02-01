#include "lmatrix_8d.h"

void write8D(FILE* f, dMAT_8D m) {
  LOOP_ALL fwrite(&GET_ALL(m), 8, 1, f);
}
