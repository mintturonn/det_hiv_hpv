#include <stdio.h>
#include <stdlib.h>
#include "popfile.h"
#include "matrix_8d.h"

void read_popfile(char* file, DArray basepop) {
  FILE* f;
  if (f = fopen(file, "rb")) {
    LOOP_A
      fread(&basepop[_da], 1, 8, f);
    fclose(f);
  } else {
    printf("Pop File not found, %s\n", file);
    exit(-1);
  }
}
