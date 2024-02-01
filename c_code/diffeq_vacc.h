#ifndef DIFFEQ_VACC_H
#define DIFFEQ_VACC_H

#include "math.h"
#include "matlab.h"
#include "lmatrix_8d.h"
#include "matrix_8d.h"
#include "params.h"
#include "results.h"
#include "state.h"
#include "tvarying.h"

results* diffeq_vacc(int tend, int tstep, dMAT_8D pop0, params* p, state* S, int* dims);

#endif
