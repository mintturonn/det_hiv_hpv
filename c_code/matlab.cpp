#include "matlab.h"

// Sum across h1

void sum_h1(dMAT_8D src, dMAT_8D dest) {
  LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1 {
    GET(dest, 1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1) = GET(src, 1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);
    for (int _dh1 = 2; _dh1 <= _SIZE_H1; _dh1++)
      GET(dest, 1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1) += GET(src, _dh1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);
  }
}

// Sum across h2, assuming h1 already summed.

void sum_h2(dMAT_8D src, dMAT_8D dest) {
  LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1 {
    GET(dest, 1, 1, _dh3, _di, _ds, _dr, _da, _dv1) = GET(src, 1, 1, _dh3, _di, _ds, _dr, _da, _dv1);
    for (int _dh2 = 2; _dh2 <= _SIZE_H2; _dh2++)
      GET(dest, 1, 1, _dh3, _di, _ds, _dr, _da, _dv1) += GET(src, 1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);
  }
}

// Sum across h3, assuming h1 and h2 are already summed.

void sum_h3(dMAT_8D src, dMAT_8D dest) {
  LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1 {
    GET(dest, 1, 1, 1, _di, _ds, _dr, _da, _dv1) = GET(src, 1, 1, 1, _di, _ds, _dr, _da, _dv1);
    for (int _dh3 = 2; _dh3 <= _SIZE_H3; _dh3++)
      GET(dest, 1, 1, 1, _di, _ds, _dr, _da, _dv1) += GET(src, 1, 1, _dh3, _di, _ds, _dr, _da, _dv1);
  }
}

// Sum across i, assuming h1, h2 and h3 are already summed.

void sum_i(dMAT_8D src, dMAT_8D dest) {
  LOOP_S LOOP_R LOOP_A LOOP_V1 {
    GET(dest, 1, 1, 1, 1, _ds, _dr, _da, _dv1) = GET(src, 1, 1, 1, 1, _ds, _dr, _da, _dv1);
      for (int _di = 2; _di <= _SIZE_I; _di++)
        GET(dest, 1, 1, 1, 1, _ds, _dr, _da, _dv1) += GET(src, 1, 1, 1, _di, _ds, _dr, _da, _dv1);
  }
}




void sum_h1_h2(dMAT_8D src, dMAT_8D dest) {
  LOOP_LINEAR GET_LIN_ALL(dest) = GET_LIN_ALL(src);
  for (int _dh1 = 2; _dh1 <= _SIZE_H1; _dh1++) LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1
    GET(dest, 1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1) +=
    GET(dest, _dh1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);

  for (int _dh2 = 2; _dh2 <= _SIZE_H2; _dh2++) LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1
    GET(dest, 1, 1, _dh3, _di, _ds, _dr, _da, _dv1) +=
    GET(dest, 1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);
}

void sum_h1_h2_h3_i(dMAT_8D src, dMAT_8D dest) {
  LOOP_LINEAR GET_LIN_ALL(dest) = GET_LIN_ALL(src);

  for (int _dh1 = 2; _dh1 <= _SIZE_H1; _dh1++) LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1
    GET(dest, 1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1) +=
    GET(dest, _dh1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);

  for (int _dh2 = 2; _dh2 <= _SIZE_H2; _dh2++) LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1
    GET(dest, 1, 1, _dh3, _di, _ds, _dr, _da, _dv1) +=
    GET(dest, 1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);

  for (int _dh3 = 2; _dh3 <= _SIZE_H3; _dh3++) LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1
    GET(dest, 1, 1, 1, _di, _ds, _dr, _da, _dv1) +=
    GET(dest, 1, 1, _dh3, _di, _ds, _dr, _da, _dv1);

  for (int _di = 2; _di <= _SIZE_I; _di++) LOOP_S LOOP_R LOOP_A LOOP_V1
    GET(dest, 1, 1, 1, 1, _ds, _dr, _da, _dv1) +=
    GET(dest, 1, 1, 1, _di, _ds, _dr, _da, _dv1);
}



void sum_h1_h2_h3(dMAT_8D src, dMAT_8D dest) {

  LOOP_LINEAR GET_LIN_ALL(dest) = GET_LIN_ALL(src);

  for (int _dh1 = 2; _dh1 <= _SIZE_H1; _dh1++) LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1
    GET(dest, 1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1) +=
    GET(dest, _dh1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);

  for (int _dh2 = 2; _dh2 <= _SIZE_H2; _dh2++) LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1
    GET(dest, 1, 1, _dh3, _di, _ds, _dr, _da, _dv1) +=
    GET(dest, 1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);

  for (int _dh3 = 2; _dh3 <= _SIZE_H3; _dh3++) LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1
    GET(dest, 1, 1, 1, _di, _ds, _dr, _da, _dv1) +=
    GET(dest, 1, 1, _dh3, _di, _ds, _dr, _da, _dv1);
}

void sum_prod_h1_h2_h3(dMAT_8D src1, dMAT_8D src2, dMAT_8D dest) {

  LOOP_LINEAR GET_LIN_ALL(dest) = GET_LIN_ALL(src1) * GET_LIN_ALL(src2);

  for (int _dh1 = 2; _dh1 <= _SIZE_H1; _dh1++) LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1
    GET(dest, 1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1) +=
    GET(dest, _dh1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);

  for (int _dh2 = 2; _dh2 <= _SIZE_H2; _dh2++) LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1
    GET(dest, 1, 1, _dh3, _di, _ds, _dr, _da, _dv1) +=
    GET(dest, 1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);

  for (int _dh3 = 2; _dh3 <= _SIZE_H3; _dh3++) LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1
    GET(dest, 1, 1, 1, _di, _ds, _dr, _da, _dv1) +=
    GET(dest, 1, 1, _dh3, _di, _ds, _dr, _da, _dv1);
}


DArray linspace(double start, double end, int steps) {
  DArray result = CREATE_DARRAY_START(steps, result) CREATE_DARRAY_END
  result[1] = start;
  result[steps] = end;
  double dstep = (end - start) / (steps - 1.0);
  for (int i = 2; i < steps; i++) {
    result[i] = start + (dstep * (double) (i-1));
  }
  return result;
}

double mlpow(double v1, double v2) {
  if ((fabs(v1)<1E-16) && (fabs(v2)<1E-16)) return 1.0;
  else return pow(v1, v2);
}

double sum_second_dim(dMAT_2D_RARA6 m, int i1) {
  double tot = 0;
  int index = GET_2D_RARA6_INDEX(i1, 1);
  for (int i=1; i <= (_SIZE_R * _SIZE_A6); i++) tot += m[index++];
  return tot;
}

double sum_ageinf(dMAT_8D tmp, int _s, int _r) {

  // Compute something like this:
  //  sum(sum(sum(sum(sum(sum(agein1(:, : , : , : , 1, 1, : , : )))))))
  //   where the supplied matrix already has the first 4 indices reduced.

  for (int _da = 2; _da <= _SIZE_A; _da++)
    LOOP_V1
      GET(tmp, 1, 1, 1, 1, _s, _r, 1, _dv1) += GET(tmp, 1, 1, 1, 1, _s, _r, _da, _dv1);

  for (int _dv1 = 2; _dv1 <= _SIZE_V1; _dv1++)
    GET(tmp, 1, 1, 1, 1, _s, _r, 1, 1) += GET(tmp, 1, 1, 1, 1, _s, _r, 1, _dv1);

  return GET(tmp,1,1,1,1,_s,_r,1,1);
}

// Sum everything except s, with specific i and a ranges, assuming h1, h2 h3 are already summed, returning array of size s.

void sum_123_4678_ia(dMAT_8D spare, dMAT_8D m, int di1, int di2, int da1, int da2, DArray result) {
  LOOP_LINEAR  GET_LIN_ALL(spare) = GET_LIN_ALL(m);

  for (int _di = di1 + 1; _di <= di2; _di++) LOOP_S LOOP_R
    for (int _da = da1; _da <= da2; _da++) LOOP_V1
      GET(spare, 1, 1, 1, di1, _ds, _dr, _da, _dv1) +=
      GET(spare, 1, 1, 1, _di, _ds, _dr, _da, _dv1);

  LOOP_S for (int _dr = 2; _dr <= _SIZE_R; _dr++)
    for (int _da = da1; _da <= da2; _da++) LOOP_V1
      GET(spare, 1, 1, 1, di1, _ds, 1, _da, _dv1) +=
      GET(spare, 1, 1, 1, di1, _ds, _dr, _da, _dv1);

  LOOP_S for (int _da = da1+1; _da <= da2; _da++) LOOP_V1
    GET(spare, 1, 1, 1, di1, _ds, 1, da1, _dv1) +=
    GET(spare, 1, 1, 1, di1, _ds, 1, _da, _dv1);

  LOOP_S for (int _dv1 = 2; _dv1 <= _SIZE_V1; _dv1++)
    GET(spare, 1, 1, 1, di1, _ds, 1, da1, 1) +=
    GET(spare, 1, 1, 1, di1, _ds, 1, da1, _dv1);

  LOOP_S result[_ds] = GET(spare, 1, 1, 1, di1, _ds, 1, da1, 1);
}

// Sum on everything except a, for given i and s ranges, assuming h1, h2, h3 are already summed, and returning array for a.

void sum_123_4568_is(dMAT_8D spare, dMAT_8D m, int di1, int di2, int ds1, int ds2, DArray result) {
  LOOP_LINEAR  GET_LIN_ALL(spare) = GET_LIN_ALL(m);

  for (int _di = di1 + 1; _di <= di2; _di++)
    for (int _ds = ds1; _ds <= ds2; _ds++) LOOP_R LOOP_A LOOP_V1
      GET(spare, 1, 1, 1, di1, _ds, _dr, _da, _dv1) +=
      GET(spare, 1, 1, 1, _di, _ds, _dr, _da, _dv1);

  for (int _ds = ds1 + 1; _ds <= ds2; _ds++) LOOP_R LOOP_A LOOP_V1
    GET(spare, 1, 1, 1, di1, ds1, _dr, _da, _dv1) +=
    GET(spare, 1, 1, 1, di1, _ds, _dr, _da, _dv1);

  for (int _dr = 2; _dr <= _SIZE_R; _dr++) LOOP_A LOOP_V1
    GET(spare, 1, 1, 1, di1, ds1, 1, _da, _dv1) +=
    GET(spare, 1, 1, 1, di1, ds1, _dr, _da, _dv1);

  LOOP_A for (int _dv1 = 2; _dv1 <= _SIZE_V1; _dv1++)
    GET(spare, 1, 1, 1, di1, ds1, 1, _da, 1) +=
    GET(spare, 1, 1, 1, di1, ds1, 1, _da, _dv1);

  LOOP_A result[_da] = GET(spare, 1, 1, 1, di1, ds1, 1, _da, 1);
}

// Sum everything except a, for given i and v1 ranges, assuming h1, h2, h3 are already summed.

void sum_123_4568_iv1(dMAT_8D spare, dMAT_8D m, int di1, int di2, int dv11, int dv12, DArray result) {
  LOOP_LINEAR  GET_LIN_ALL(spare) = GET_LIN_ALL(m);

  for (int _di = di1 + 1; _di <= di2; _di++) LOOP_S LOOP_R LOOP_A
    for (int _dv1 = dv11; _dv1 <= dv12; _dv1++)
      GET(spare, 1, 1, 1, di1, _ds, _dr, _da, _dv1) +=
      GET(spare, 1, 1, 1, _di, _ds, _dr, _da, _dv1);

  for (int _ds = 2; _ds <= _SIZE_S; _ds++) LOOP_R LOOP_A
    for (int _dv1 = dv11; _dv1 <= dv12; _dv1++)
      GET(spare, 1, 1, 1, di1, 1, _dr, _da, _dv1) +=
      GET(spare, 1, 1, 1, di1, _ds, _dr, _da, _dv1);

  for (int _dr = 2; _dr <= _SIZE_R; _dr++) LOOP_A
    for (int _dv1 = dv11; _dv1 <= dv12; _dv1++)
      GET(spare, 1, 1, 1, di1, 1, 1, _da, _dv1) +=
      GET(spare, 1, 1, 1, di1, 1, _dr, _da, _dv1);

  for (int _da = 2; _da <= _SIZE_A; _da++)
    for (int _dv1 = dv11; _dv1 <= dv12; _dv1++)
      GET(spare, 1, 1, 1, di1, 1, 1, 1, _dv1) +=
      GET(spare, 1, 1, 1, di1, 1, 1, _da, _dv1);

  for (int _dv1 = dv11 + 1; _dv1 <= dv12; _dv1++)
    GET(spare, 1, 1, 1, di1, 1, 1, 1, dv11) +=
    GET(spare, 1, 1, 1, di1, 1, 1, 1, _dv1);

  LOOP_A result[_da] = GET(spare, 1, 1, 1, di1, 1, 1, _da, dv11);
}

// Sum everything except a, with specified i, s and v1, assuming h1, h2, h3 already summed, returning an array across a
void sum_123_4568_isv1(dMAT_8D spare, dMAT_8D m, int di1, int di2, int ds1, int ds2, int dv11, int dv12, DArray result) {
  LOOP_LINEAR GET_LIN_ALL(spare) = GET_LIN_ALL(m);

  for (int _di = di1 + 1; _di <= di2; _di++)
    for (int _ds = ds1; _ds <= ds2; _ds++) LOOP_R LOOP_A
      for (int _dv1 = dv11; _dv1 <= dv12; _dv1++)
        GET(spare, 1, 1, 1, di1, _ds, _dr, _da, _dv1) +=
        GET(spare, 1, 1, 1, _di, _ds, _dr, _da, _dv1);

  for (int _ds = ds1 + 1; _ds <= ds2; _ds++) LOOP_R LOOP_A
    for (int _dv1 = dv11; _dv1 <= dv12; _dv1++)
      GET(spare, 1, 1, 1, di1, ds1, _dr, _da, _dv1) +=
      GET(spare, 1, 1, 1, di1, _ds, _dr, _da, _dv1);

  for (int _dr = 2; _dr <= _SIZE_R; _dr++) LOOP_A
    for (int _dv1 = dv11; _dv1 <= dv12; _dv1++)
      GET(spare, 1, 1, 1, di1, ds1, 1, _da, _dv1) +=
      GET(spare, 1, 1, 1, di1, ds1, _dr, _da, _dv1);

  LOOP_A for (int _dv1 = dv11 + 1; _dv1 <= dv12; _dv1++)
      GET(spare, 1, 1, 1, di1, ds1, 1, _da, dv11) +=
      GET(spare, 1, 1, 1, di1, ds1, 1, _da, _dv1);

  LOOP_A result[_da] = GET(spare, 1, 1, 1, di1, ds1, 1, _da, dv11);
}

// Sum everything for particular i, s, r, a, where h1, h2 and h3 are already summed.

double sum_all_isra_ch123(dMAT_8D spare, dMAT_8D m, int di1, int di2, int ds1, int ds2, int dr1, int dr2, int da1, int da2) {
  LOOP_LINEAR GET_LIN_ALL(spare) = GET_LIN_ALL(m);
  for (int _di = di1 + 1; _di <= di2; _di++)
    for (int _ds = ds1; _ds <= ds2; _ds++)
      for (int _dr = dr1; _dr <= dr2; _dr++)
        for (int _da = da1; _da <= da2; _da++) LOOP_V1
          GET(spare, 1, 1, 1, di1, _ds, _dr, _da, _dv1) +=
          GET(spare, 1, 1, 1, _di, _ds, _dr, _da, _dv1);

  for (int _ds = ds1 + 1; _ds <= ds2; _ds++)
    for (int _dr = dr1; _dr <= dr2; _dr++)
      for (int _da = da1; _da <= da2; _da++) LOOP_V1
        GET(spare, 1, 1, 1, di1, ds1, _dr, _da, _dv1) +=
        GET(spare, 1, 1, 1, di1, _ds, _dr, _da, _dv1);

  for (int _dr = dr1 + 1; _dr <= dr2; _dr++)
    for (int _da = da1; _da <= da2; _da++) LOOP_V1
      GET(spare, 1, 1, 1, di1, ds1, dr1, _da, _dv1) +=
      GET(spare, 1, 1, 1, di1, ds1, _dr, _da, _dv1);

  for (int _da = da1 + 1; _da <= da2; _da++) LOOP_V1
    GET(spare, 1, 1, 1, di1, ds1, dr1, da1, _dv1) +=
    GET(spare, 1, 1, 1, di1, ds1, dr1, _da, _dv1);

  for (int _dv1 = 2; _dv1 <= _SIZE_V1; _dv1++)
    GET(spare, 1, 1, 1, di1, ds1, dr1, da1, 1) +=
    GET(spare, 1, 1, 1, di1, ds1, dr1, da1, _dv1);

  return GET(spare, 1, 1, 1, di1, ds1, dr1, da1, 1);
}

void sum_hivsus_pop(dMAT_8D spare, dMAT_8D m, DMatrix3D result) {
  // Special case-
  //      hivsus_pop(gg, :, : , : , : ) = reshape(sum(sum(sum(sum(popn(:, : , : , 1, : , : , : , : ), 1), 2), 3), 6), 1, s, a, v1);

  LOOP_LINEAR GET_LIN_ALL(spare) = GET_LIN_ALL(m);

  for (int _dh1 = 2; _dh1 <= _SIZE_H1; _dh1++) LOOP_H2 LOOP_H3 LOOP_S LOOP_R LOOP_A LOOP_V1
       GET(spare, 1,    _dh2, _dh3, 1, _ds, _dr, _da, _dv1) +=
       GET(spare, _dh1, _dh2, _dh3, 1, _ds, _dr, _da, _dv1);

  for (int _dh2 = 2; _dh2 <= _SIZE_H2; _dh2++) LOOP_H3 LOOP_S LOOP_R LOOP_A LOOP_V1
      GET(spare, 1, 1,    _dh3, 1, _ds, _dr, _da, _dv1) +=
      GET(spare, 1, _dh2, _dh3, 1, _ds, _dr, _da, _dv1);

  for (int _dh3 = 2; _dh3 <= _SIZE_H3; _dh3++) LOOP_S LOOP_R LOOP_A LOOP_V1
      GET(spare, 1, 1, 1,    1, _ds, _dr, _da, _dv1) +=
      GET(spare, 1, 1, _dh3, 1, _ds, _dr, _da, _dv1);

  LOOP_S for (int _dr = 2; _dr <= _SIZE_R; _dr++) LOOP_A LOOP_V1
    GET(spare, 1, 1, 1, 1, _ds, 1, _da, _dv1) +=
    GET(spare, 1, 1, 1, 1, _ds, _dr, _da, _dv1);

  LOOP_S LOOP_A LOOP_V1 result[_ds][_da][_dv1] = GET(spare, 1, 1, 1, 1, _ds, 1, _da, _dv1);

}

//      ccinc_nvtft(n, :) = ccinc_nvtft(n - 1, :) + reshape(sum(sum(sum(sum(sum(sum(ccinc_pop(:, : , 4, : , : , : , : , : ), 1), 2), 4), 5), 6), 8), 1, a);
void sum_12_4568_h3is(dMAT_8D spare, dMAT_8D m, int dh31, int dh32, int di1, int di2, int ds1, int ds2, DArray result) {
  LOOP_LINEAR GET_LIN_ALL(spare) = GET_LIN_ALL(m);

  for (int _dh3 = dh31 + 1; _dh3 <= dh32; _dh3++)
    for (int _di = di1; _di <= di2; _di++)
      for (int _ds = ds1; _ds <= ds2; _ds++) LOOP_R LOOP_A LOOP_V1
        GET(spare, 1, 1, dh31, _di, _ds, _dr, _da, _dv1) +=
        GET(spare, 1, 1, _dh3, _di, _ds, _dr, _da, _dv1);

  for (int _di = di1 + 1; _di <= di2; _di++)
    for (int _ds = ds1; _ds <= ds2; _ds++) LOOP_R LOOP_A LOOP_V1
      GET(spare, 1, 1, dh31, di1, _ds, _dr, _da, _dv1) +=
      GET(spare, 1, 1, dh31, _di, _ds, _dr, _da, _dv1);

  for (int _ds = ds1 + 1; _ds <= ds2; _ds++) LOOP_R LOOP_A LOOP_V1
    GET(spare, 1, 1, dh31, di1, ds1, _dr, _da, _dv1) +=
    GET(spare, 1, 1, dh31, di1, _ds, _dr, _da, _dv1);

  for (int _dr = 2; _dr <= _SIZE_R; _dr++) LOOP_A LOOP_V1
    GET(spare, 1, 1, dh31, di1, ds1, 1, _da, _dv1) +=
    GET(spare, 1, 1, dh31, di1, ds1, _dr, _da, _dv1);

  LOOP_A for (int _dv1 = 2; _dv1 <= _SIZE_V1; _dv1++)
    GET(spare, 1, 1, dh31, di1, ds1, 1, _da, 1) +=
    GET(spare, 1, 1, dh31, di1, ds1, 1, _da, _dv1);

  LOOP_A result[_da] = GET(spare, 1, 1, dh31, di1, ds1, 1, _da, 1);

}

void sum_1234_68_s(dMAT_8D spare, dMAT_8D m, int ds, DArray result) {
  LOOP_LINEAR GET_LIN_ALL(spare) = GET_LIN_ALL(m);

  for (int _dr = 2; _dr <= _SIZE_R;_dr++) LOOP_A LOOP_V1
    GET(spare, 1, 1, 1, 1, ds, 1, _da, _dv1) +=
    GET(spare, 1, 1, 1, 1, ds, _dr, _da, _dv1);

    LOOP_A for (int _dv1 = 2; _dv1 <= _SIZE_V1; _dv1++)
      GET(spare, 1, 1, 1, 1, ds, 1, _da, 1) +=
      GET(spare, 1, 1, 1, 1, ds, 1, _da, _dv1);

  LOOP_A result[_da] = GET(spare, 1, 1, 1, 1, ds, 1, _da, 1);
}

// Sum across v1, for given s, r, a, returning array of s,r,a, assuming h1, h2, h3, i, already summed.

void sum_1234_8_sra_a6(DArray spare_v1, dMAT_8D m, int ds1, int ds2, int dr1, int dr2, int da1, int da2, dMAT_3D_SRA6 result, int target_a) {
  for (int _ds = ds1; _ds <= ds2; _ds++)
    for (int _dr = dr1; _dr <= dr2; _dr++) {
      LOOP_V1 {
        spare_v1[_dv1] = GET(m, 1, 1, 1, 1, _ds, _dr, da1, _dv1);
        for (int _da = da1 + 1; _da <= da2; _da++)
          spare_v1[_dv1] += GET(m, 1, 1, 1, 1, _ds, _dr, _da, _dv1);
      }
      GET_3D_SRA6(result, _ds, _dr, target_a) = spare_v1[1];
      for (int _dv1 = 2; _dv1 <= _SIZE_V1; _dv1++) GET_3D_SRA6(result, _ds, _dr, target_a) += spare_v1[_dv1];
    }
}

// Sum across v1, for given i, s, r, a, returning array of s, r, a, assuming h1, h2, h3 already summed.

void sum123_8_isr_a(DArray spare_v1, dMAT_8D m, int di1, int di2, int ds1, int ds2, int dr1, int dr2, int da1, int da2, dMAT_4D_ISRA6 result, int target_a) {
  for (int _di = di1; _di <= di2; _di++)
    for (int _ds = ds1; _ds <= ds2; _ds++)
      for (int _dr = dr1; _dr <= dr2; _dr++) {
        LOOP_V1 {
          spare_v1[_dv1] = GET(m, 1, 1, 1, _di, _ds, _dr, da1, _dv1);
          for (int _da = da1 + 1; _da <= da2; _da++)
            spare_v1[_dv1] += GET(m, 1, 1, 1, _di, _ds, _dr, _da, _dv1);
        }

        GET_4D_ISRA6(result, _di, _ds, _dr, target_a) = spare_v1[1];
        for (int _dv1 = 2; _dv1 <= _SIZE_V1; _dv1++) GET_4D_ISRA6(result, _di, _ds, _dr, target_a) += spare_v1[_dv1];
      }
}

void sum12_48_h3sr_a(dMAT_2D_AV spare_av, dMAT_8D m, int dh31, int dh32, int ds1, int ds2, int dr1, int dr2, int da1, int da2, dMAT_4D_H3SRA6 result, int target_a) {
  for (int _dh3 = dh31; _dh3 <= dh32; _dh3++)
    for (int _ds = ds1; _ds <= ds2; _ds++)
      for (int _dr = dr1; _dr <= dr2; _dr++) {
        LOOP_A LOOP_V1 {
          GET_2D_AV(spare_av, _da, _dv1) = GET(m, 1, 1, _dh3, 1, _ds, _dr, _da, _dv1);
          for (int _di = 2; _di <= _SIZE_I; _di++)
            GET_2D_AV(spare_av, _da, _dv1) += GET(m, 1, 1, _dh3, _di, _ds, _dr, _da, _dv1);
        }
        LOOP_V1 {
          for (int _da = da1 + 1; _da <= da2; _da++)
            GET_2D_AV(spare_av, da1, _dv1) += GET_2D_AV(spare_av, _da, _dv1);
        }
        GET_4D_H3SRA6(result, _dh3, _ds, _dr, target_a) = GET_2D_AV(spare_av, da1, 1);
        for (int _dv1 = 2; _dv1 <= _SIZE_V1; _dv1++)
          GET_4D_H3SRA6(result, _dh3, _ds, _dr, target_a) += GET_2D_AV(spare_av, da1, _dv1);
      }
}

void sum1_348_h2sr_a(dMAT_3D_IAV spare_iav, dMAT_8D m, int dh21, int dh22, int ds1, int ds2, int dr1, int dr2, int da1, int da2, dMAT_4D_H2SRA6 result, int target_a) {
  for (int _dh2 = dh21; _dh2 <= dh22; _dh2++)
    for (int _ds = ds1; _ds <= ds2; _ds++)
      for (int _dr = dr1; _dr <= dr2; _dr++) {
        LOOP_I LOOP_A LOOP_V1 {
          GET_3D_IAV(spare_iav, _di, _da, _dv1) = GET(m, 1, _dh2, 1, _di, _ds, _dr, _da, _dv1);
          for (int _dh3 = 2; _dh3 <= _SIZE_H3; _dh3++)
            GET_3D_IAV(spare_iav, _di, _da, _dv1) += GET(m, 1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);
        }

        LOOP_A LOOP_V1
          for (int _di = 2; _di <= _SIZE_I; _di++)
            GET_3D_IAV(spare_iav, 1, _da, _dv1) += GET_3D_IAV(spare_iav, _di, _da, _dv1);

        for (int _da = da1 + 1; _da <= da2; _da++)
          LOOP_V1
            GET_3D_IAV(spare_iav, 1, da1, _dv1) += GET_3D_IAV(spare_iav, 1, _da, _dv1);

        GET_4D_H2SRA6(result, _dh2, _ds, _dr, target_a) = GET_3D_IAV(spare_iav, 1, da1, 1);

        for (int _dv1 = 2; _dv1 <= _SIZE_V1; _dv1++)
          GET_4D_H2SRA6(result, _dh2, _ds, _dr, target_a) += GET_3D_IAV(spare_iav, 1, da1, _dv1);
      }
}

void sum_2348_h1sr_a(dMAT_4D_H3IAV spare_h3iav, dMAT_8D m, int dh11, int dh12, int ds1, int ds2, int dr1, int dr2, int da1, int da2, dMAT_4D_H1SRA6 result, int target_a) {
  for (int _dh1 = dh11; _dh1 <= dh12; _dh1++)
    for (int _ds = ds1; _ds <= ds2; _ds++)
      for (int _dr = dr1; _dr <= dr2; _dr++) {
        LOOP_H3 LOOP_I LOOP_A LOOP_V1 {
          GET_4D_H3IAV(spare_h3iav, _dh3, _di, _da, _dv1) = GET(m, _dh1, 1, _dh3, _di, _ds, _dr, _da, _dv1);
          for (int _dh2 = 2; _dh2 <= _SIZE_H2; _dh2++)
            GET_4D_H3IAV(spare_h3iav, _dh3, _di, _da, _dv1) += GET(m, _dh1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1);
        }

        for (int _dh3 = 2; _dh3 <= _SIZE_H3; _dh3++)
          LOOP_I LOOP_A LOOP_V1
            GET_4D_H3IAV(spare_h3iav, 1, _di, _da, _dv1) += GET_4D_H3IAV(spare_h3iav, _dh3, _di, _da, _dv1);

        for (int _di = 2; _di <= _SIZE_I; _di++)
          LOOP_A LOOP_V1
            GET_4D_H3IAV(spare_h3iav, 1, 1, _da, _dv1) += GET_4D_H3IAV(spare_h3iav, 1, _di, _da, _dv1);

        for (int _da = da1+1; _da <= da2; _da++)
          LOOP_V1
            GET_4D_H3IAV(spare_h3iav, 1, 1, da1, _dv1) += GET_4D_H3IAV(spare_h3iav, 1, 1, _da, _dv1);

        GET_4D_H1SRA6(result, _dh1, _ds, _dr, target_a) = GET_4D_H3IAV(spare_h3iav, 1, 1, da1, 1);
        for (int _dv1 = 2; _dv1 <= _SIZE_V1; _dv1++)
          GET_4D_H1SRA6(result, _dh1, _ds, _dr, target_a) += GET_4D_H3IAV(spare_h3iav, 1, 1, da1, _dv1);
      }
}

// Sum everything, specifying ds and dv, return a. Assuming h1, h2, h3 i summed already.

void sum_1234_568_sv1(dMAT_8D spare, dMAT_8D m, int ds, int dv1, DArray result) {
  LOOP_LINEAR GET_LIN_ALL(spare) = GET_LIN_ALL(m);
  for (int _dr = 2; _dr <= _SIZE_R; _dr++) LOOP_A
    GET(spare, 1, 1, 1, 1, ds, 1, _da, dv1) += GET(spare, 1, 1, 1, 1, ds, _dr, _da, dv1);

  LOOP_A  result[_da] = GET(spare,1, 1, 1, 1, ds, 1, _da, dv1);
}

// Sum everything, return s, assumg h1,h2,h3 i  summed already.
void sum_all_1234_s(dMAT_8D spare, dMAT_8D m, DArray result) {
  LOOP_LINEAR GET_LIN_ALL(spare) = GET_LIN_ALL(m);
  LOOP_S {
    for (int _dr = 2; _dr <= _SIZE_R; _dr++) LOOP_A LOOP_V1
      GET(spare, 1, 1, 1, 1, _ds, 1, _da, _dv1) += GET(spare, 1, 1, 1, 1, _ds, _dr, _da, _dv1);
    for (int _da = 2; _da <= _SIZE_A; _da++) LOOP_V1
      GET(spare, 1, 1, 1, 1, _ds, 1, 1, _dv1) += GET(spare, 1, 1, 1, 1, _ds, 1, _da, _dv1);
    for (int _dv1 = 2; _dv1 <= _SIZE_V1; _dv1++)
      GET(spare, 1, 1, 1, 1, _ds, 1, 1, 1) += GET(spare, 1, 1, 1, 1, _ds, 1, 1, _dv1);
    result[_ds] = GET(spare, 1, 1, 1, 1, _ds, 1, 1, 1);
  }
}

void sum_1234_568_s(dMAT_8D spare, dMAT_8D m, int ds, DArray result) {
  LOOP_LINEAR GET_LIN_ALL(spare) = GET_LIN_ALL(m);
  for (int _dr = 2; _dr <= _SIZE_R; _dr++) LOOP_A LOOP_V1
    GET(spare, 1, 1, 1, 1, ds, 1, _da, _dv1) += GET(spare, 1, 1, 1, 1, ds, _dr, _da, _dv1);

  LOOP_A for (int _dv1 = 2; _dv1 <= _SIZE_V1; _dv1++)
    GET(spare, 1, 1, 1, 1, ds, 1, _da, 1) += GET(spare, 1, 1, 1, 1, ds, 1, _da, _dv1);

  LOOP_A  result[_da] = GET(spare, 1, 1, 1, 1, ds, 1, _da, 1);
}

#define LOOP_DI for (int _di = _di1; _di <= _di2; _di++)

void sum_prod_2346_578(dMAT_4D_H2H3IR spare, dMAT_8D foi, dMAT_8D pop, dMAT_4D_NSAV1 res, int n, int _di1, int _di2, double dt) {
  int _dh1 = 1;
  LOOP_S LOOP_A LOOP_V1 {
      LOOP_H2 LOOP_H3 LOOP_DI LOOP_R GET_4D_H2H3IR(spare, _dh2, _dh3, _di, _dr) = GET_ALL(foi) * GET_ALL(pop);
      for (int _dh2 = 2; _dh2 <= _SIZE_H2; _dh2++) LOOP_H3 LOOP_DI LOOP_R GET_4D_H2H3IR(spare, 1, _dh3, _di, _dr) += GET_4D_H2H3IR(spare, _dh2, _dh3, _di, _dr);
      for (int _dh3 = 2; _dh3 <= _SIZE_H3; _dh3++) LOOP_DI LOOP_R         GET_4D_H2H3IR(spare, 1, 1, _di, _dr)    += GET_4D_H2H3IR(spare, 1, _dh3, _di, _dr);
      for (int _di = _di1 + 1; _di <= _di2; _di++) LOOP_R                 GET_4D_H2H3IR(spare, 1, 1, _di1, _dr)      += GET_4D_H2H3IR(spare, 1, 1, _di, _dr);
      double sum = 0;
      LOOP_R sum += GET_4D_H2H3IR(spare, 1, 1, _di1, _dr);
      GET_4D_NSAV1(res, n, _ds, _da, _dv1) = GET_4D_NSAV1(res, n - 1 , _ds, _da, _dv1) + (sum * dt);
  }
}

void sum_prod_1346_578(dMAT_4D_H1H3IR spare, dMAT_8D foi, dMAT_8D pop, dMAT_4D_NSAV1 res, int n, int _di1, int _di2, double dt) {
  int _dh2 = 1;
  LOOP_S LOOP_A LOOP_V1 {
      LOOP_H1 LOOP_H3 LOOP_DI LOOP_R GET_4D_H2H3IR(spare, _dh1, _dh3, _di, _dr) = GET_ALL(foi) * GET_ALL(pop);
      for (int _dh1 = 2; _dh1 <= _SIZE_H1; _dh1++) LOOP_H3 LOOP_DI LOOP_R GET_4D_H1H3IR(spare, 1, _dh3, _di, _dr) += GET_4D_H1H3IR(spare, _dh1, _dh3, _di, _dr);
      for (int _dh3 = 2; _dh3 <= _SIZE_H3; _dh3++) LOOP_DI LOOP_R         GET_4D_H1H3IR(spare, 1, 1, _di, _dr)    += GET_4D_H1H3IR(spare, 1, _dh3, _di, _dr);
      for (int _di = _di1 + 1; _di <= _di2; _di++) LOOP_R                 GET_4D_H1H3IR(spare, 1, 1, _di1, _dr)   += GET_4D_H1H3IR(spare, 1, 1, _di, _dr);
      double sum = 0;
      LOOP_R sum += GET_4D_H1H3IR(spare, 1, 1, _di1, _dr);
      GET_4D_NSAV1(res, n , _ds, _da, _dv1) = GET_4D_NSAV1(res, n - 1, _ds, _da, _dv1) + (sum * dt);
  }
}

void sum_prod_1246_578(dMAT_4D_H1H2IR spare, dMAT_8D foi, dMAT_8D pop, dMAT_4D_NSAV1 res, int n, int _di1, int _di2, double dt) {
  int _dh3 = 1;
  LOOP_S LOOP_A LOOP_V1 {
      LOOP_H1 LOOP_H2 LOOP_DI LOOP_R GET_4D_H1H2IR(spare, _dh1, _dh2, _di, _dr) = GET_ALL(foi) * GET_ALL(pop);
      for (int _dh1 = 2; _dh1 <= _SIZE_H1; _dh1++) LOOP_H2 LOOP_DI LOOP_R GET_4D_H1H2IR(spare, 1, _dh2, _di, _dr) += GET_4D_H1H2IR(spare, _dh1, _dh2, _di, _dr);
      for (int _dh2 = 2; _dh2 <= _SIZE_H2; _dh2++) LOOP_DI LOOP_R         GET_4D_H1H2IR(spare, 1, 1,    _di, _dr) += GET_4D_H1H2IR(spare, 1,    _dh2, _di, _dr);
      for (int _di = _di1 + 1; _di <= _di2; _di++) LOOP_R                 GET_4D_H1H2IR(spare, 1, 1,   _di1, _dr) += GET_4D_H1H2IR(spare, 1,       1, _di, _dr);
      double sum = 0;
      LOOP_R sum += GET_4D_H1H2IR(spare, 1, 1, _di1, _dr);
      GET_4D_NSAV1(res, n, _ds, _da, _dv1) = GET_4D_NSAV1(res, n - 1, _ds, _da, _dv1) + (sum * dt);
  }
}
