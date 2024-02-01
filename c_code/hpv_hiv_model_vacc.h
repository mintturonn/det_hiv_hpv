#ifndef HPV_HIV_MODEL_VACC_H
#define HPV_HIV_MODEL_VACC_H

#include "diffeq_vacc.h"
#include "params.h"
#include "popfile.h"
#include "results.h"
#include "matrix_8d.h"
#include "lmatrix_8d.h"
#include "math.h"
#include "state.h"


void hiv_hpv_model_vacc(params* p, char* pop_file, char* output_file);
#endif
