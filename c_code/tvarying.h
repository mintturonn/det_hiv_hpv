#ifndef TVARYING_H
#define TVARYING_H

#include "math.h"
#include "matlab.h"
#include "matrix_8d.h"

DArray tvarying_plat(int tlengthtot, int tvlength, double tmin, double tmax, double grate, int grplat);
DArray tvarying_artcon(int tlengthtot, int tvlength, int tvlength2, double tmin, double tmax);
DArray tvarying_theta(int tlengthtot, int tyr1, int tyr2, int tyr3, int tyr4, int tyr5,
                      double ttheta1, double ttheta2, double ttheta3, double ttheta4, double ttheta5);
DArray tvarying_theta_add(int tlengthtot, int tyr1, int tyr2, int tyr3, int tyr4, int tyr5, int tyr6, int tyr7, int tyr8,
                          double ttheta1, double ttheta2, double ttheta3, double ttheta4, double ttheta5, double ttheta6,
                          double ttheta7, double ttheta8);


#endif
