#include "tvarying.h"

/****************************************************************************************/
/* tvarying_plat.m */

DArray tvarying_plat(int tlengthtot, int tvlength, double tmin, double tmax, double grate, int grplat) {

  //  % tlengthtot = total parameter vector
  //  % tvlength = timevarying component of the parameter vector
  //  % tmax = max of parameter
  //  % tmin = min / start of the parameter
  //  % grate = growth rate

  //  % %change starts 1 tstep after tchange as tchange(1) is still defined as the starting conditions

  //  tfinal = zeros(tlengthtot, 1);

    // Going to skip tfinal here, since we fill it later.

  //  % tfinal = linspace(tmax, tmin, tlengthtot);

  //  ttime = linspace(0, 1, tlengthtot - tvlength);
  //  tchange = tvlength;

  DArray ttime = linspace(0, 1, tlengthtot - tvlength);
  int tchange = tvlength;

  //  if tmin == 0
  //    tmin2 = 10 ^ -5;
  //  else
  //    tmin2 = tmin;
  //  end

  double tmin2 = (tmin == 0) ? (1E-5) : tmin;

  //    tlogistic(:, 1) = (tmax * tmin2.*exp(grate.*ttime)). / (tmax + tmin2.*(exp(grate.*ttime) - 1));

  DArray tlogistic = CREATE_DARRAY(tlengthtot, tlogistic, (tmax * tmin2 * exp(grate * ttime[d1])) / (tmax + tmin2 * (exp(grate * ttime[d1]) - 1.0)))

    //  tfinal(1:tchange) = tmin;
    //  tfinal((tchange + 1) : end) = tlogistic(:, 1);

      /* Not quite sure what the above means - what does ":" mean in tlogistic? Does it match index with tfinal(...) - or does it start from 1? */

    DArray tfinal = CREATE_DARRAY(tlengthtot, tfinal, tmin)
    for (int i = tchange + 1; i <= tlengthtot; i++) tfinal[i] = tlogistic[i];

  //  % to restrict the end, to plateu for the last
  //    plateuphase = linspace(tfinal(grplat) / 1.2, 0.5, tlengthtot - grplat);

  DArray plateuphase = linspace(tfinal[grplat] / 1.2, 0.5, tlengthtot - grplat);

  //  tfinal(grplat + 1:end) = plateuphase;%;

  for (int i = grplat + 1; i <= tlengthtot; i++) tfinal[i] = plateuphase[i - grplat];

  //  tvars = tfinal;

  return tfinal;
}


/***************************************************************/
/* tvarying_artcon.m */

//function tvars = tvarying_artcon(tlengthtot, tvlength, tvlength2, tmin, tmax)

DArray tvarying_artcon(int tlengthtot, int tvlength, int tvlength2, double tmin, double tmax) {

  //  % tlengthtot = total parameter vector
  //  % tvlength = timevarying component of the parameter vector
  //  % tmax = max of parameter
  //  % tmin = min / start of the parameter
  //  % grate = growth rate

  //  % %change starts 1 tstep after tchange as tchange(1) is still defined as the starting conditions

  //tfinal = zeros(tlengthtot, 1);
  //ttime = tvlength2 - tvlength;
  //tchange = linspace(tmin, tmax, ttime);
  //tfinal(1:tvlength) = tmin;

  DArray tfinal = CREATE_DARRAY(tlengthtot, tfinal, tmin)
  DArray tchange = linspace(tmin, tmax, tvlength2 - tvlength);

  //tfinal(tvlength + 1:tvlength2) = tchange;
  //tfinal(tvlength2 + 1:end) = tmax;

  for (int i = tvlength + 1; i <= tvlength2; i++) tfinal[i] = tchange[i - tvlength];
  for (int i = tvlength2 + 1; i <= tlengthtot; i++) tfinal[i] = tmax;

  //tvars = tfinal;
  //end

  FREE_DARRAY(tchange)
  return tfinal;
}

/***************************************************************/
/* tvarying_theta.m */


// function tvars = tvarying_theta(tlengthtot, tyr1, tyr2, tyr3, tyr4, tyr5, ttheta1, ttheta2, ttheta3, ttheta4, ttheta5)

DArray tvarying_theta(int tlengthtot, int tyr1, int tyr2, int tyr3, int tyr4, int tyr5,
                      double ttheta1, double ttheta2, double ttheta3, double ttheta4, double ttheta5) {

//   tfinal = zeros(tlengthtot, 1);

  DArray tfinal = CREATE_DARRAY(tlengthtot, tfinal, 0.0);

//   tfinal(tyr1:tyr2) = ttheta1;% here there might be a bit of a jump

  for (int i = tyr1; i <= tyr2; i++) tfinal[i] = ttheta1;

//   tfinal(tyr2 + 1:tyr3) = linspace(ttheta2, ttheta3, length(tyr2 + 1:tyr3));

  DArray lin2_3 = linspace(ttheta2, ttheta3, tyr3 - tyr2);
  for (int i = 1; i <= (tyr3 - tyr2); i++) tfinal[tyr2 + i] = lin2_3[i];

//   tfinal(tyr3 + 1:tyr4) = linspace(ttheta3, ttheta4, length(tyr3 + 1:tyr4));

  DArray lin3_4 = linspace(ttheta3, ttheta4, tyr4 - tyr3);
  for (int i = 1; i <= (tyr4 - tyr3); i++) tfinal[tyr3 + i] = lin3_4[i];

//   tfinal(tyr4 + 1:tyr5) = linspace(ttheta4, ttheta5, length(tyr4 + 1:tyr5));

  DArray lin4_5 = linspace(ttheta4, ttheta5, tyr5 - tyr4);
  for (int i = 1; i <= (tyr5 - tyr4); i++) tfinal[tyr4 + i] = lin4_5[i];

//   tfinal(tyr5 + 1:end) = ttheta5;

  for (int i= tyr5 + 1; i <= tlengthtot; i++) tfinal[i] = ttheta5;

//   tvars = tfinal;
// end

  FREE_DARRAY(lin2_3)
  FREE_DARRAY(lin3_4)
  FREE_DARRAY(lin4_5)
  return tfinal;
}


//function tvars = tvarying_theta_add(tlengthtot, tyr1, tyr2, tyr3, tyr4, tyr5, tyr6, tyr7, tyr8, ttheta1, ttheta2, ttheta3, ttheta4, ttheta5, ttheta6, ttheta7, ttheta8)

DArray tvarying_theta_add(int tlengthtot, int tyr1, int tyr2, int tyr3, int tyr4, int tyr5, int tyr6, int tyr7, int tyr8,
                          double ttheta1, double ttheta2, double ttheta3, double ttheta4, double ttheta5, double ttheta6, double ttheta7, double ttheta8) {

// tfinal = zeros(tlengthtot, 1);
  DArray tfinal = CREATE_DARRAY(tlengthtot, tfinal, 0.0);


  // tfinal(tyr1:tyr2) = ttheta1;% here there might be a bit of a jump
  for (int i=tyr1; i <= tyr2; i++) tfinal[i] = ttheta1;

  // tfinal(tyr2 + 1:tyr3) = linspace(ttheta2, ttheta3, length(tyr2 + 1:tyr3));
  DArray lin2_3 = linspace(ttheta2, ttheta3, tyr3 - tyr2);
  for (int i=1; i <= (tyr3 - tyr2); i++) tfinal[tyr2 + i] = lin2_3[i];

  // tfinal(tyr3 + 1:tyr4) = linspace(ttheta3, ttheta4, length(tyr3 + 1:tyr4));
  DArray lin3_4 = linspace(ttheta3, ttheta4, tyr4 - tyr3);
  for (int i = 1; i <= (tyr4 - tyr3); i++) tfinal[tyr3 + i] = lin3_4[i];

  // tfinal(tyr4 + 1:tyr5) = linspace(ttheta4, ttheta5, length(tyr4 + 1:tyr5));
  DArray lin4_5 = linspace(ttheta4, ttheta5, tyr5 - tyr4);
  for (int i = 1; i <= (tyr5 - tyr4); i++) tfinal[tyr4 + i] = lin4_5[i];

  // tfinal(tyr5 + 1:tyr6) = linspace(ttheta5, ttheta6, length(tyr5 + 1:tyr6));
  DArray lin5_6 = linspace(ttheta5, ttheta6, tyr6 - tyr5);
  for (int i = 1; i <= (tyr6 - tyr5); i++) tfinal[tyr5 + i] = lin5_6[i];

  // tfinal(tyr6 + 1:tyr7) = linspace(ttheta6, ttheta7, length(tyr6 + 1:tyr7));
  DArray lin6_7 = linspace(ttheta6, ttheta7, tyr7 - tyr6);
  for (int i = 1; i <= (tyr7 - tyr6); i++) tfinal[tyr6 + i] = lin6_7[i];

  // tfinal(tyr7 + 1:tyr8) = linspace(ttheta7, ttheta8, length(tyr7 + 1:tyr8));
  DArray lin7_8 = linspace(ttheta7, ttheta8, tyr8 - tyr7);
  for (int i = 1; i <= (tyr8 - tyr7); i++) tfinal[tyr7 + i] = lin7_8[i];

  // tfinal(tyr8 + 1:end) = ttheta8;
  for (int i = tyr8 + 1; i <= tlengthtot; i++) tfinal[i] = ttheta8;


  // tvars = tfinal;
  // end

  delete[] lin2_3;
  delete[] lin3_4;
  delete[] lin4_5;
  delete[] lin5_6;
  delete[] lin6_7;
  delete[] lin7_8;

  return tfinal;
}
