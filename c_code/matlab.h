#ifndef MATLAB_H
#define MATLAB_H
  #include "matrix_8d.h"
  #include "math.h"
  #include "lmatrix_8d.h"
  #include "params.h"

  void sum_h1(dMAT_8D src, dMAT_8D dest);
  void sum_h2(dMAT_8D src, dMAT_8D dest);
  void sum_h3(dMAT_8D src, dMAT_8D dest);
  void sum_i(dMAT_8D src, dMAT_8D dest);

  void sum_h1_h2(dMAT_8D src, dMAT_8D dest);
  void sum_h1_h2_h3_i(dMAT_8D src, dMAT_8D dest);
  void sum_h1_h2_h3(dMAT_8D src, dMAT_8D dest);
  void sum_prod_h1_h2_h3(dMAT_8D src1, dMAT_8D src2, dMAT_8D dest);

  DArray linspace(double start, double end, int steps);
  double sum_second_dim(dMAT_2D_RARA6  m, int i1);
  double sum_ageinf(dMAT_8D tmp, int _s, int _r);

  double sum_all_isra_ch123(dMAT_8D spare, dMAT_8D m, int di1, int di2, int ds1, int ds2, int dr1, int dr2, int da1, int da2);
  void sum_hivsus_pop(dMAT_8D spare, dMAT_8D m, DMatrix3D result);

  void sum_123_4678_ia(dMAT_8D spare, dMAT_8D m, int di1, int di2, int da1, int da2, DArray result);
  void sum_123_4568_is(dMAT_8D spare, dMAT_8D m, int di1, int di2, int ds1, int ds2, DArray result);
  void sum_123_4568_iv1(dMAT_8D spare, dMAT_8D m, int di1, int di2, int dv11, int dv12, DArray result);
  void sum_123_4568_isv1(dMAT_8D spare, dMAT_8D m, int di1, int di2, int ds1, int ds2, int dv11, int dv12, DArray result);
  void sum_12_4568_h3is(dMAT_8D spare, dMAT_8D m, int dh31, int dh32, int di1, int di2, int ds1, int ds2, DArray result);
  void sum_1234_68_s(dMAT_8D spare, dMAT_8D m, int ds, DArray result);
  void sum_1234_8_sra_a6(DArray spare_v1, dMAT_8D m, int ds1, int ds2, int dr1, int dr2, int da1, int da2, dMAT_3D_SRA6 result, int target_a);
  void sum123_8_isr_a(DArray spare_v1, dMAT_8D m, int di1, int di2, int ds1, int ds2, int dr1, int dr2, int da1, int da2, dMAT_4D_ISRA6 result, int target_a);
  void sum12_48_h3sr_a(dMAT_2D_AV spare_av, dMAT_8D m, int dh31, int dh32, int ds1, int ds2, int dr1, int dr2, int da1, int da2, dMAT_4D_H3SRA6 result, int target_a);
  void sum1_348_h2sr_a(dMAT_3D_IAV spare_iav, dMAT_8D m, int dh21, int dh22, int ds1, int ds2, int dr1, int dr2, int da1, int da2, dMAT_4D_H2SRA6 result, int target_a);
  void sum_2348_h1sr_a(dMAT_4D_H3IAV spare_h3iav, dMAT_8D m, int dh11, int dh12, int ds1, int ds2, int dr1, int dr2, int da1, int da2, dMAT_4D_H1SRA6 result, int target_a);
  void sum_1234_568_sv1(dMAT_8D spare, dMAT_8D m, int ds, int dv1, DArray result);
  void sum_all_1234_s(dMAT_8D spare, dMAT_8D m, DArray result);
  void sum_1234_568_s(dMAT_8D spare, dMAT_8D m, int ds, DArray result);
  void sum_prod_2346_578(dMAT_4D_H2H3IR spare, dMAT_8D foi, dMAT_8D pop, dMAT_4D_NSAV1 res, int n, int _di1, int _di2, double dt);
  void sum_prod_1346_578(dMAT_4D_H1H3IR spare, dMAT_8D foi, dMAT_8D pop, dMAT_4D_NSAV1 res, int n, int _di1, int _di2, double dt);
  void sum_prod_1246_578(dMAT_4D_H1H2IR spare, dMAT_8D foi, dMAT_8D pop, dMAT_4D_NSAV1 res, int n, int _di1, int _di2, double dt);

  // Expand a matrix of size 18 x 18 (ie, R*A6 squared), into one 39 x 18 (ie, R*A x R*A6) - for use with tempfoih1_full etc.
  #define EXPAND_MATRIX_A6_A6(MOUT, MIN) {\
    int p = 0; \
    int q = 0; \
    for (int _x = 1; _x <= 3 * _SIZE_R * _SIZE_A6; _x++) MOUT[p++] = MIN[q++]; \
    for (int _y = 1; _y <= 2; _y++) { \
      for (int _x = 1; _x <= 3 * _SIZE_R * _SIZE_A6; _x++) { \
        MOUT[p] = MIN[q]; \
        MOUT[p + (3 * _SIZE_R * _SIZE_A6)] = MIN[q]; \
        p++; \
        q++; \
      } \
      p += (3 * _SIZE_R * _SIZE_A6); \
    } \
    for (int _x = 1; _x <= 3 * _SIZE_R * _SIZE_A6; _x++) { \
      MOUT[p] = MIN[q]; \
      MOUT[p + (3 * _SIZE_R * _SIZE_A6)] = MIN[q]; \
      MOUT[p + (6 * _SIZE_R * _SIZE_A6)] = MIN[q]; \
      p++; \
      q++; \
    } \
    p += (6 * _SIZE_R * _SIZE_A6); \
    for (int _x = 1; _x <= 3 * _SIZE_R * _SIZE_A6; _x++) { \
      MOUT[p] = MIN[q]; \
      MOUT[p + (3 * _SIZE_R * _SIZE_A6)] = MIN[q]; \
      p++; \
      q++; \
    } \
    p += (3 * _SIZE_R * _SIZE_A6); \
    for (int _x = 1; _x <= 3 * _SIZE_R * _SIZE_A6; _x++) { \
      MOUT[p] = MIN[q]; \
      MOUT[p + (3 * _SIZE_R * _SIZE_A6)] = MIN[q]; \
      MOUT[p + (6 * _SIZE_R * _SIZE_A6)] = MIN[q]; \
      p++; \
      q++; \
    } \
  }

  // Expand an array of size 18 x 1 (ie, R*A6), into an array 39 x 1 (ie, R*A6) - for use with tempfoih1_comm_full etc.
#define EXPAND_MATRIX_A6_1(MOUT, MIN) \
  MOUT[1]  = MIN[1];   MOUT[2] = MIN[2];   MOUT[3] = MIN[3]; \
  MOUT[4]  = MIN[4];   MOUT[5] = MIN[5];   MOUT[6] = MIN[6];   MOUT[7] = MIN[4];   MOUT[8] = MIN[5];   MOUT[9] = MIN[6]; \
  MOUT[10] = MIN[7];  MOUT[11] = MIN[8];  MOUT[12] = MIN[9];  MOUT[13] = MIN[7];  MOUT[14] = MIN[8];  MOUT[15] = MIN[9]; \
  MOUT[16] = MIN[10]; MOUT[17] = MIN[11]; MOUT[18] = MIN[12]; MOUT[19] = MIN[10]; MOUT[20] = MIN[11]; MOUT[21] = MIN[12]; MOUT[22] = MIN[10]; MOUT[23] = MIN[11]; MOUT[24] = MIN[12]; \
  MOUT[25] = MIN[13]; MOUT[26] = MIN[14]; MOUT[27] = MIN[15]; MOUT[28] = MIN[13]; MOUT[29] = MIN[14]; MOUT[30] = MIN[15]; \
  MOUT[31] = MIN[16]; MOUT[32] = MIN[17]; MOUT[33] = MIN[18]; MOUT[34] = MIN[16]; MOUT[35] = MIN[17]; MOUT[36] = MIN[18]; MOUT[37] = MIN[16]; MOUT[38] = MIN[17]; MOUT[39] = MIN[18];

#endif
