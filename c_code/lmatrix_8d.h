#ifndef LMATRIX8D_H
#define LMATRIX8D_H

#include <stdio.h>

#define _NITER 62050
#define _SIZE_H1 5
#define _SIZE_H2 5
#define _SIZE_H3 5
#define _SIZE_I 4
#define _SIZE_S 2
#define _SIZE_R 3
#define _SIZE_A 13
#define _SIZE_V1 2

#define _SIZE_A6 6

#define LOOP_H1 for (int _dh1 = 1; _dh1 <= _SIZE_H1; _dh1++)
#define LOOP_H2 for (int _dh2 = 1; _dh2 <= _SIZE_H2; _dh2++)
#define LOOP_H3 for (int _dh3 = 1; _dh3 <= _SIZE_H3; _dh3++)
#define LOOP_I  for (int _di = 1; _di <= _SIZE_I; _di++)
#define LOOP_S  for (int _ds = 1; _ds <= _SIZE_S; _ds++)
#define LOOP_R  for (int _dr = 1; _dr <= _SIZE_R; _dr++)
#define LOOP_A  for (int _da = 1; _da <= _SIZE_A; _da++)
#define LOOP_V1 for (int _dv1 = 1; _dv1 <= _SIZE_V1; _dv1++)
#define LOOP_ALL LOOP_H1 LOOP_H2 LOOP_H3 LOOP_I LOOP_S LOOP_R LOOP_A LOOP_V1

#define LOOP_A6  for (int _da = 1; _da <= _SIZE_A6; _da++)

// Pre-calculate some positions. Start with V1 (2), multiply by A, then R, etc...

#define _PP7 _SIZE_V1
#define _PP6 (_PP7 * _SIZE_A)
#define _PP5 (_PP6 * _SIZE_R)
#define _PP4 (_PP5 * _SIZE_S)
#define _PP3 (_PP4 * _SIZE_I)
#define _PP2 (_PP3 * _SIZE_H3)
#define _PP1 (_PP2 * _SIZE_H2)
#define _PRODS (_PP1 * _SIZE_H1)

#define PAR_LOOP_LINEAR for (int _thread=0; _thread < p->threads; _thread++) \
for (int _dall = _thread; _dall < _PRODS; _dall += p->threads)

#define LOOP_LINEAR for (int _dall = 0; _dall < _PRODS; _dall++)

#define dMAT_8D double*
#define dNEW_COPY_8D(_NEW, _OLD) new double[_PRODS]; for (int __z = 0; __z<_PRODS; __z++) _NEW[__z] = _OLD[__z];
#define dNEW_8D(_V) new double[_PRODS]; for (int __z = 0; __z<_PRODS; __z++) _V[__z] = 0.0;

#define GET_INDEX(D1, D2, D3, D4, D5, D6, D7, D8) (D8-1) + ((D7-1)* _PP7) + ((D6-1)*_PP6) + ((D5-1)*_PP5) + ((D4-1)*_PP4) + ((D3-1)*_PP3) + ((D2-1)*_PP2) + ((D1-1)*_PP1)
#define GET(M, D1, D2, D3, D4, D5, D6, D7, D8) M[GET_INDEX(D1,D2,D3,D4,D5,D6,D7,D8)]
#define SET(M, D1, D2, D3, D4, D5, D6, D7, D8, V) GET(M,D1,D2,D3,D4,D5,D6,D7,D8) = V;
#define GET_ALL(V) GET(V, _dh1, _dh2, _dh3, _di, _ds, _dr, _da, _dv1)

#define cMAT_8D char*
#define cNEW_COPY_8D(_NEW, _OLD) new char[_PRODS]; for (int __z = 0; __z<_PRODS; __z++) _NEW[__z] = _OLD[__z];
#define cNEW_8D(_V) new char[_PRODS]; for (int __z = 0; __z<_PRODS; __z++) _V[__z] = 0;

#define IS(VAR,VAL) (VAL == GET_ALL(VAR))
#define LT(VAR,VAL) (GET_ALL(VAR) < VAL)
#define GT(VAR,VAL) (GET_ALL(VAR) > VAL)

#define IS_LIN(VAR,VAL) (VAL == VAR[_dall])
#define LT_LIN(VAR,VAL) (VAR[_dall] < VAL)
#define GT_LIN(VAR,VAL) (VAR[_dall] > VAL)

#define GET_LIN_ALL(V) V[_dall]

void write8D(FILE* f, dMAT_8D m);

// 3-D matrices with dimensions S, R and A (with 13 age groups).

#define _SRA_P2 _SIZE_A
#define _SRA_P1 (_SRA_P2 * _SIZE_R)
#define _SRA_PRODS (_SRA_P1 * _SIZE_S)

#define dMAT_3D_SRA double*
#define dNEW_COPY_3D_SRA(_NEW, _OLD) new double[_SRA_PRODS]; for (int __z = 0; __z<_SRA_PRODS; __z++) _NEW[__z] = _OLD[__z];
#define dNEW_3D_SRA(_V) new double[_SRA_PRODS]; for (int __z = 0; __z<_SRA_PRODS; __z++) _V[__z] = 0.0;

#define GET_3D_SRA_INDEX(D1, D2, D3) (D3-1) + ((D2-1)*_SRA_P2) + ((D1-1)*_SRA_P1)
#define GET_3D_SRA(M, D1, D2, D3) M[GET_3D_SRA_INDEX(D1,D2,D3)]

// 3-D matrices with dimensions S, R and A - with 6 age groups

#define _SRA6_P2 _SIZE_A6
#define _SRA6_P1 (_SRA6_P2 * _SIZE_R)
#define _SRA6_PRODS (_SRA6_P1 * _SIZE_S)

#define dMAT_3D_SRA6 double*
#define dNEW_COPY_3D_SRA6(_NEW, _OLD) new double[_SRA6_PRODS]; for (int __z = 0; __z<_SRA6_PRODS; __z++) _NEW[__z] = _OLD[__z];
#define dNEW_3D_SRA6(_V) new double[_SRA6_PRODS]; for (int __z = 0; __z<_SRA6_PRODS; __z++) _V[__z] = 0.0;

#define GET_3D_SRA6_INDEX(D1, D2, D3) (D3-1) + ((D2-1)*_SRA6_P2) + ((D1-1)*_SRA6_P1)
#define GET_3D_SRA6(M, D1, D2, D3) M[GET_3D_SRA6_INDEX(D1,D2,D3)]
#define GET_ALL_3D_SRA6(V) GET_3D_SRA6(V, _ds, _dr, _da)
#define LOOP_LINEAR_3D_SRA6 for (int _dall = 0; _dall < _SRA6_PRODS; _dall++)

void write3D_SRA(FILE f, dMAT_3D_SRA m);

// 4-D matrices, with dimensions I, S, R, A6

#define _ISRA6_P3 _SIZE_A6
#define _ISRA6_P2 (_ISRA6_P3 * _SIZE_R)
#define _ISRA6_P1 (_ISRA6_P2 * _SIZE_S)
#define _ISRA6_PRODS (_ISRA6_P1 * _SIZE_I)

#define dMAT_4D_ISRA6 double*
#define dNEW_4D_ISRA6(_V) new double[_ISRA6_PRODS]; for (int __z = 0; __z<_ISRA6_PRODS; __z++) _V[__z] = 0.0;

#define GET_4D_ISRA6_INDEX(D1, D2, D3, D4) (D4-1) + ((D3-1)*_ISRA6_P3) + ((D2-1)*_ISRA6_P2) + ((D1-1)*_ISRA6_P1)
#define GET_4D_ISRA6(M, D1, D2, D3, D4) M[GET_4D_ISRA6_INDEX(D1,D2,D3,D4)]

// 4-D matrices, with dimensions H3, S, R, A6

#define _H3SRA6_P3 _SIZE_A6
#define _H3SRA6_P2 (_H3SRA6_P3 * _SIZE_R)
#define _H3SRA6_P1 (_H3SRA6_P2 * _SIZE_S)
#define _H3SRA6_PRODS (_H3SRA6_P1 * _SIZE_H3)

#define dMAT_4D_H3SRA6 double*
#define dNEW_4D_H3SRA6(_V) new double[_H3SRA6_PRODS]; for (int __z = 0; __z<_H3SRA6_PRODS; __z++) _V[__z] = 0.0;

#define GET_4D_H3SRA6_INDEX(D1, D2, D3, D4) (D4-1) + ((D3-1)*_H3SRA6_P3) + ((D2-1)*_H3SRA6_P2) + ((D1-1)*_H3SRA6_P1)
#define GET_4D_H3SRA6(M, D1, D2, D3, D4) M[GET_4D_H3SRA6_INDEX(D1,D2,D3,D4)]

// 4-D matrices, with dimensions H2, S, R, A6

#define _H2SRA6_P3 _SIZE_A6
#define _H2SRA6_P2 (_H2SRA6_P3 * _SIZE_R)
#define _H2SRA6_P1 (_H2SRA6_P2 * _SIZE_S)
#define _H2SRA6_PRODS (_H2SRA6_P1 * _SIZE_H2)

#define dMAT_4D_H2SRA6 double*
#define dNEW_4D_H2SRA6(_V) new double[_H2SRA6_PRODS]; for (int __z = 0; __z<_H2SRA6_PRODS; __z++) _V[__z] = 0.0;

#define GET_4D_H2SRA6_INDEX(D1, D2, D3, D4) (D4-1) + ((D3-1)*_H2SRA6_P3) + ((D2-1)*_H2SRA6_P2) + ((D1-1)*_H2SRA6_P1)
#define GET_4D_H2SRA6(M, D1, D2, D3, D4) M[GET_4D_H2SRA6_INDEX(D1,D2,D3,D4)]

// 4-D matrices, with dimensions H1, S, R, A6

#define _H1SRA6_P3 _SIZE_A6
#define _H1SRA6_P2 (_H1SRA6_P3 * _SIZE_R)
#define _H1SRA6_P1 (_H1SRA6_P2 * _SIZE_S)
#define _H1SRA6_PRODS (_H1SRA6_P1 * _SIZE_H1)

#define dMAT_4D_H1SRA6 double*
#define dNEW_4D_H1SRA6(_V) new double[_H1SRA6_PRODS]; for (int __z = 0; __z<_H1SRA6_PRODS; __z++) _V[__z] = 0.0;

#define GET_4D_H1SRA6_INDEX(D1, D2, D3, D4) (D4-1) + ((D3-1)*_H1SRA6_P3) + ((D2-1)*_H1SRA6_P2) + ((D1-1)*_H1SRA6_P1)
#define GET_4D_H1SRA6(M, D1, D2, D3, D4) M[GET_4D_H1SRA6_INDEX(D1,D2,D3,D4)]

// 2-D matrices, square with dimensions (r*a6), (r*a6)

#define _RA6RA6_P1 (_SIZE_R * _SIZE_A6)
#define _RA6RA6_PRODS (_RA6RA6_P1 * _RA6RA6_P1)

#define dMAT_2D_RA6RA6 double*
#define dNEW_2D_RA6RA6(_V) new double[_RA6RA6_PRODS]; for (int __z = 0; __z<_RA6RA6_PRODS; __z++) _V[__z] = 0.0;

#define GET_2D_RA6RA6_INDEX(D1, D2) (D2-1) + ((D1-1)*_RA6RA6_P1)
#define GET_2D_RA6RA6(M, D1, D2) M[GET_2D_RA6RA6_INDEX(D1,D2)]

// 2-D matrices, square with dimensions (r*a), (r*a6)

#define _RARA6_P1 (_SIZE_R * _SIZE_A6)
#define _RARA6_PRODS (_SIZE_R * _SIZE_A * _RARA6_P1)

#define dMAT_2D_RARA6 double*
#define dNEW_2D_RARA6(_V) new double[_RARA6_PRODS]; for (int __z = 0; __z<_RARA6_PRODS; __z++) _V[__z] = 0.0;

#define GET_2D_RARA6_INDEX(D1, D2) (D2-1) + ((D1-1)*_RARA6_P1)
#define GET_2D_RARA6(M, D1, D2) M[GET_2D_RARA6_INDEX(D1,D2)]


// 2-D matrices, niterations * age groups (n=62050, a = 5)

#define _NA_P1 (_SIZE_A)
#define _NA_PRODS (_NA_P1 * _NITER)

#define dMAT_2D_NA double*
#define dNEW_2D_NA(_V) new double[_NA_PRODS]; for (int __z = 0; __z<_NA_PRODS; __z++) _V[__z] = 0.0;

#define GET_2D_NA_INDEX(D1, D2) (D2-1) + ((D1-1)*_NA_P1)
#define GET_2D_NA(M, D1, D2) M[GET_2D_NA_INDEX(D1,D2)]
#define LOOP_LINEAR_2D_NA for (int _dall = 0; _dall < _NA_PRODS; _dall++)

// 2-D matrices, A * V

#define _AV_P1 (_SIZE_V1)
#define _AV_PRODS (_SIZE_V1 * _SIZE_A)

#define dMAT_2D_AV double*
#define dNEW_2D_AV(_V) new double[_AV_PRODS]; for (int __z = 0; __z<_AV_PRODS; __z++) _V[__z] = 0.0;

#define GET_2D_AV_INDEX(D1, D2) (D2-1) + ((D1-1)*_AV_P1)
#define GET_2D_AV(M, D1, D2) M[GET_2D_AV_INDEX(D1,D2)]
#define LOOP_LINEAR_AV_NA for (int _dall = 0; _dall < _AV_PRODS; _dall++)

// 2-D matrices, S x (R * A6)

#define _S_RA6_P1 (_SIZE_R * _SIZE_A6)
#define _S_RA6_PRODS (_SIZE_S * _S_RA6_P1)

#define dMAT_2D_S_RA6 double*
#define dNEW_2D_S_RA6(_V) new double[_S_RA6_PRODS]; for (int __z = 0; __z<_S_RA6_PRODS; __z++) _V[__z] = 0.0;

#define GET_2D_S_RA6_INDEX(D1, D2) (D2-1) + ((D1-1)*_S_RA6_P1)
#define GET_2D_S_RA6(M, D1, D2) M[GET_2D_S_RA6_INDEX(D1,D2)]

// 4-D matrix, niterations, s, a, v1

#define _NSAV1_P3 _SIZE_V1
#define _NSAV1_P2 (_NSAV1_P3 * _SIZE_A)
#define _NSAV1_P1 (_NSAV1_P2 * _SIZE_S)
#define _NSAV1_PRODS (_NSAV1_P1 * _NITER)

#define dMAT_4D_NSAV1 double*
#define dNEW_4D_NSAV1(_V) new double[_NSAV1_PRODS]; for (int __z = 0; __z<_NSAV1_PRODS; __z++) _V[__z] = 0.0;

#define GET_4D_NSAV1_INDEX(D1, D2, D3, D4) (D4-1) + ((D3-1)*_NSAV1_P3) + ((D2-1)*_NSAV1_P2) + ((D1-1)*_NSAV1_P1)
#define GET_4D_NSAV1(M, D1, D2, D3, D4) M[GET_4D_NSAV1_INDEX(D1,D2,D3,D4)]
#define LOOP_LINEAR_4D_NSAV1 for (int _dall = 0; _dall < _NSAV1_PRODS; _dall++)

// 4-D matrix, niterations, h3, i, a, v1

#define _H3IAV_P3 _SIZE_V1
#define _H3IAV_P2 (_H3IAV_P3 * _SIZE_A)
#define _H3IAV_P1 (_H3IAV_P2 * _SIZE_I)
#define _H3IAV_PRODS (_H3IAV_P1 * _SIZE_H3)

#define dMAT_4D_H3IAV double*
#define dNEW_4D_H3IAV(_V) new double[_H3IAV_PRODS]; for (int __z = 0; __z<_H3IAV_PRODS; __z++) _V[__z] = 0.0;

#define GET_4D_H3IAV_INDEX(D1, D2, D3, D4) (D4-1) + ((D3-1)*_H3IAV_P3) + ((D2-1)*_H3IAV_P2) + ((D1-1)*_H3IAV_P1)
#define GET_4D_H3IAV(M, D1, D2, D3, D4) M[GET_4D_H3IAV_INDEX(D1,D2,D3,D4)]
#define LOOP_LINEAR_4D_H3IAV for (int _dall = 0; _dall < _H3IAV_PRODS; _dall++)

// 4-D matrix, niterations, h2, h3, i, r

#define _H2H3IR_P3 _SIZE_R
#define _H2H3IR_P2 (_H2H3IR_P3 * _SIZE_I)
#define _H2H3IR_P1 (_H2H3IR_P2 * _SIZE_H3)
#define _H2H3IR_PRODS (_H2H3IR_P1 * _SIZE_H2)

#define dMAT_4D_H2H3IR double*
#define dNEW_4D_H2H3IR(_V) new double[_H2H3IR_PRODS]; for (int __z = 0; __z<_H2H3IR_PRODS; __z++) _V[__z] = 0.0;

#define GET_4D_H2H3IR_INDEX(D1, D2, D3, D4) (D4-1) + ((D3-1)*_H2H3IR_P3) + ((D2-1)*_H2H3IR_P2) + ((D1-1)*_H2H3IR_P1)
#define GET_4D_H2H3IR(M, D1, D2, D3, D4) M[GET_4D_H2H3IR_INDEX(D1,D2,D3,D4)]

// 4-D matrix, niterations, h1, h3, i, r

#define _H1H3IR_P3 _SIZE_R
#define _H1H3IR_P2 (_H1H3IR_P3 * _SIZE_I)
#define _H1H3IR_P1 (_H1H3IR_P2 * _SIZE_H3)
#define _H1H3IR_PRODS (_H1H3IR_P1 * _SIZE_H1)

#define dMAT_4D_H1H3IR double*
#define dNEW_4D_H1H3IR(_V) new double[_H1H3IR_PRODS]; for (int __z = 0; __z<_H1H3IR_PRODS; __z++) _V[__z] = 0.0;

#define GET_4D_H1H3IR_INDEX(D1, D2, D3, D4) (D4-1) + ((D3-1)*_H1H3IR_P3) + ((D2-1)*_H1H3IR_P2) + ((D1-1)*_H1H3IR_P1)
#define GET_4D_H1H3IR(M, D1, D2, D3, D4) M[GET_4D_H1H3IR_INDEX(D1,D2,D3,D4)]

// 4-D matrix, niterations, h1, h2, i, r

#define _H1H2IR_P3 _SIZE_R
#define _H1H2IR_P2 (_H1H2IR_P3 * _SIZE_I)
#define _H1H2IR_P1 (_H1H2IR_P2 * _SIZE_H2)
#define _H1H2IR_PRODS (_H1H2IR_P1 * _SIZE_H1)

#define dMAT_4D_H1H2IR double*
#define dNEW_4D_H1H2IR(_V) new double[_H1H2IR_PRODS]; for (int __z = 0; __z<_H1H2IR_PRODS; __z++) _V[__z] = 0.0;

#define GET_4D_H1H2IR_INDEX(D1, D2, D3, D4) (D4-1) + ((D3-1)*_H1H2IR_P3) + ((D2-1)*_H1H2IR_P2) + ((D1-1)*_H1H2IR_P1)
#define GET_4D_H1H2IR(M, D1, D2, D3, D4) M[GET_4D_H1H2IR_INDEX(D1,D2,D3,D4)]

// 3-D Matrix for screencov.

#define _N32_P2 2
#define _N32_P1 6
#define _N32_PRODS (_N32_P1 * _NITER)

#define dMAT_3D_N32 double*
#define dNEW_3D_N32(_V) new double[_N32_PRODS]; for (int __z = 0; __z<_N32_PRODS; __z++) _V[__z] = 0.0;

#define GET_3D_N32_INDEX(D1, D2, D3) (D3-1) + ((D2-1)*_N32_P2) + ((D1-1)*_N32_P1)
#define GET_3D_N32(M, D1, D2, D3) M[GET_3D_N32_INDEX(D1,D2,D3)]

// 3-D Matrix for IAV calculation.

#define _IAV_P2 _SIZE_V1
#define _IAV_P1 (_SIZE_V1 * _SIZE_A)
#define _IAV_PRODS (_IAV_P1 * _SIZE_I)

#define dMAT_3D_IAV double*
#define dNEW_3D_IAV(_V) new double[_IAV_PRODS]; for (int __z = 0; __z<_IAV_PRODS; __z++) _V[__z] = 0.0;

#define GET_3D_IAV_INDEX(D1, D2, D3) (D3-1) + ((D2-1)*_IAV_P2) + ((D1-1)*_IAV_P1)
#define GET_3D_IAV(M, D1, D2, D3) M[GET_3D_IAV_INDEX(D1,D2,D3)]


#endif
