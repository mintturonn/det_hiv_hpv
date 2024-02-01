// HPV_HIV
// Conversion From Matlab Code
// Wes Hinsley, December 2019

#include "omp.h"
#include "params.h"
#include "popfile.h"

#include "hpv_hiv_model_vacc.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
  if (argc != 6) {
    printf("Usage: \n");
    printf("        HPV_HIV path/to/params.bin path/to/popfile.bin path/to/scenario.bin path/to/outputfile.bin no_threads\n");
    printf("  eg.   HPV_HIV params.bin popfile.bin scenario.bin output.bin 8\n");
    printf("        or to use as many threads a the machine has, set no_threads to 0");
    exit(-1);
  }

  params* p = new params(argv[1], argv[3]);
  p->threads = atoi(argv[5]);
  if (p->threads == 0) p->threads = omp_get_max_threads();
  omp_set_num_threads(p->threads);
  hiv_hpv_model_vacc(p, argv[2], argv[4]);
  //fflush(stdout);
}
