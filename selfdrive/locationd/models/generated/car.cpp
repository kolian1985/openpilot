#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7193476806311499376) {
   out_7193476806311499376[0] = delta_x[0] + nom_x[0];
   out_7193476806311499376[1] = delta_x[1] + nom_x[1];
   out_7193476806311499376[2] = delta_x[2] + nom_x[2];
   out_7193476806311499376[3] = delta_x[3] + nom_x[3];
   out_7193476806311499376[4] = delta_x[4] + nom_x[4];
   out_7193476806311499376[5] = delta_x[5] + nom_x[5];
   out_7193476806311499376[6] = delta_x[6] + nom_x[6];
   out_7193476806311499376[7] = delta_x[7] + nom_x[7];
   out_7193476806311499376[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8103533003814411399) {
   out_8103533003814411399[0] = -nom_x[0] + true_x[0];
   out_8103533003814411399[1] = -nom_x[1] + true_x[1];
   out_8103533003814411399[2] = -nom_x[2] + true_x[2];
   out_8103533003814411399[3] = -nom_x[3] + true_x[3];
   out_8103533003814411399[4] = -nom_x[4] + true_x[4];
   out_8103533003814411399[5] = -nom_x[5] + true_x[5];
   out_8103533003814411399[6] = -nom_x[6] + true_x[6];
   out_8103533003814411399[7] = -nom_x[7] + true_x[7];
   out_8103533003814411399[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_4590973562022669193) {
   out_4590973562022669193[0] = 1.0;
   out_4590973562022669193[1] = 0;
   out_4590973562022669193[2] = 0;
   out_4590973562022669193[3] = 0;
   out_4590973562022669193[4] = 0;
   out_4590973562022669193[5] = 0;
   out_4590973562022669193[6] = 0;
   out_4590973562022669193[7] = 0;
   out_4590973562022669193[8] = 0;
   out_4590973562022669193[9] = 0;
   out_4590973562022669193[10] = 1.0;
   out_4590973562022669193[11] = 0;
   out_4590973562022669193[12] = 0;
   out_4590973562022669193[13] = 0;
   out_4590973562022669193[14] = 0;
   out_4590973562022669193[15] = 0;
   out_4590973562022669193[16] = 0;
   out_4590973562022669193[17] = 0;
   out_4590973562022669193[18] = 0;
   out_4590973562022669193[19] = 0;
   out_4590973562022669193[20] = 1.0;
   out_4590973562022669193[21] = 0;
   out_4590973562022669193[22] = 0;
   out_4590973562022669193[23] = 0;
   out_4590973562022669193[24] = 0;
   out_4590973562022669193[25] = 0;
   out_4590973562022669193[26] = 0;
   out_4590973562022669193[27] = 0;
   out_4590973562022669193[28] = 0;
   out_4590973562022669193[29] = 0;
   out_4590973562022669193[30] = 1.0;
   out_4590973562022669193[31] = 0;
   out_4590973562022669193[32] = 0;
   out_4590973562022669193[33] = 0;
   out_4590973562022669193[34] = 0;
   out_4590973562022669193[35] = 0;
   out_4590973562022669193[36] = 0;
   out_4590973562022669193[37] = 0;
   out_4590973562022669193[38] = 0;
   out_4590973562022669193[39] = 0;
   out_4590973562022669193[40] = 1.0;
   out_4590973562022669193[41] = 0;
   out_4590973562022669193[42] = 0;
   out_4590973562022669193[43] = 0;
   out_4590973562022669193[44] = 0;
   out_4590973562022669193[45] = 0;
   out_4590973562022669193[46] = 0;
   out_4590973562022669193[47] = 0;
   out_4590973562022669193[48] = 0;
   out_4590973562022669193[49] = 0;
   out_4590973562022669193[50] = 1.0;
   out_4590973562022669193[51] = 0;
   out_4590973562022669193[52] = 0;
   out_4590973562022669193[53] = 0;
   out_4590973562022669193[54] = 0;
   out_4590973562022669193[55] = 0;
   out_4590973562022669193[56] = 0;
   out_4590973562022669193[57] = 0;
   out_4590973562022669193[58] = 0;
   out_4590973562022669193[59] = 0;
   out_4590973562022669193[60] = 1.0;
   out_4590973562022669193[61] = 0;
   out_4590973562022669193[62] = 0;
   out_4590973562022669193[63] = 0;
   out_4590973562022669193[64] = 0;
   out_4590973562022669193[65] = 0;
   out_4590973562022669193[66] = 0;
   out_4590973562022669193[67] = 0;
   out_4590973562022669193[68] = 0;
   out_4590973562022669193[69] = 0;
   out_4590973562022669193[70] = 1.0;
   out_4590973562022669193[71] = 0;
   out_4590973562022669193[72] = 0;
   out_4590973562022669193[73] = 0;
   out_4590973562022669193[74] = 0;
   out_4590973562022669193[75] = 0;
   out_4590973562022669193[76] = 0;
   out_4590973562022669193[77] = 0;
   out_4590973562022669193[78] = 0;
   out_4590973562022669193[79] = 0;
   out_4590973562022669193[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_738104619053175027) {
   out_738104619053175027[0] = state[0];
   out_738104619053175027[1] = state[1];
   out_738104619053175027[2] = state[2];
   out_738104619053175027[3] = state[3];
   out_738104619053175027[4] = state[4];
   out_738104619053175027[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_738104619053175027[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_738104619053175027[7] = state[7];
   out_738104619053175027[8] = state[8];
}
void F_fun(double *state, double dt, double *out_7442358059839394283) {
   out_7442358059839394283[0] = 1;
   out_7442358059839394283[1] = 0;
   out_7442358059839394283[2] = 0;
   out_7442358059839394283[3] = 0;
   out_7442358059839394283[4] = 0;
   out_7442358059839394283[5] = 0;
   out_7442358059839394283[6] = 0;
   out_7442358059839394283[7] = 0;
   out_7442358059839394283[8] = 0;
   out_7442358059839394283[9] = 0;
   out_7442358059839394283[10] = 1;
   out_7442358059839394283[11] = 0;
   out_7442358059839394283[12] = 0;
   out_7442358059839394283[13] = 0;
   out_7442358059839394283[14] = 0;
   out_7442358059839394283[15] = 0;
   out_7442358059839394283[16] = 0;
   out_7442358059839394283[17] = 0;
   out_7442358059839394283[18] = 0;
   out_7442358059839394283[19] = 0;
   out_7442358059839394283[20] = 1;
   out_7442358059839394283[21] = 0;
   out_7442358059839394283[22] = 0;
   out_7442358059839394283[23] = 0;
   out_7442358059839394283[24] = 0;
   out_7442358059839394283[25] = 0;
   out_7442358059839394283[26] = 0;
   out_7442358059839394283[27] = 0;
   out_7442358059839394283[28] = 0;
   out_7442358059839394283[29] = 0;
   out_7442358059839394283[30] = 1;
   out_7442358059839394283[31] = 0;
   out_7442358059839394283[32] = 0;
   out_7442358059839394283[33] = 0;
   out_7442358059839394283[34] = 0;
   out_7442358059839394283[35] = 0;
   out_7442358059839394283[36] = 0;
   out_7442358059839394283[37] = 0;
   out_7442358059839394283[38] = 0;
   out_7442358059839394283[39] = 0;
   out_7442358059839394283[40] = 1;
   out_7442358059839394283[41] = 0;
   out_7442358059839394283[42] = 0;
   out_7442358059839394283[43] = 0;
   out_7442358059839394283[44] = 0;
   out_7442358059839394283[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_7442358059839394283[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_7442358059839394283[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7442358059839394283[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7442358059839394283[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_7442358059839394283[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_7442358059839394283[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_7442358059839394283[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_7442358059839394283[53] = -9.8000000000000007*dt;
   out_7442358059839394283[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_7442358059839394283[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_7442358059839394283[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7442358059839394283[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7442358059839394283[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_7442358059839394283[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_7442358059839394283[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_7442358059839394283[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7442358059839394283[62] = 0;
   out_7442358059839394283[63] = 0;
   out_7442358059839394283[64] = 0;
   out_7442358059839394283[65] = 0;
   out_7442358059839394283[66] = 0;
   out_7442358059839394283[67] = 0;
   out_7442358059839394283[68] = 0;
   out_7442358059839394283[69] = 0;
   out_7442358059839394283[70] = 1;
   out_7442358059839394283[71] = 0;
   out_7442358059839394283[72] = 0;
   out_7442358059839394283[73] = 0;
   out_7442358059839394283[74] = 0;
   out_7442358059839394283[75] = 0;
   out_7442358059839394283[76] = 0;
   out_7442358059839394283[77] = 0;
   out_7442358059839394283[78] = 0;
   out_7442358059839394283[79] = 0;
   out_7442358059839394283[80] = 1;
}
void h_25(double *state, double *unused, double *out_5800239379178490244) {
   out_5800239379178490244[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7167661522270920967) {
   out_7167661522270920967[0] = 0;
   out_7167661522270920967[1] = 0;
   out_7167661522270920967[2] = 0;
   out_7167661522270920967[3] = 0;
   out_7167661522270920967[4] = 0;
   out_7167661522270920967[5] = 0;
   out_7167661522270920967[6] = 1;
   out_7167661522270920967[7] = 0;
   out_7167661522270920967[8] = 0;
}
void h_24(double *state, double *unused, double *out_3363033020538636619) {
   out_3363033020538636619[0] = state[4];
   out_3363033020538636619[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4990447098663770994) {
   out_4990447098663770994[0] = 0;
   out_4990447098663770994[1] = 0;
   out_4990447098663770994[2] = 0;
   out_4990447098663770994[3] = 0;
   out_4990447098663770994[4] = 1;
   out_4990447098663770994[5] = 0;
   out_4990447098663770994[6] = 0;
   out_4990447098663770994[7] = 0;
   out_4990447098663770994[8] = 0;
   out_4990447098663770994[9] = 0;
   out_4990447098663770994[10] = 0;
   out_4990447098663770994[11] = 0;
   out_4990447098663770994[12] = 0;
   out_4990447098663770994[13] = 0;
   out_4990447098663770994[14] = 1;
   out_4990447098663770994[15] = 0;
   out_4990447098663770994[16] = 0;
   out_4990447098663770994[17] = 0;
}
void h_30(double *state, double *unused, double *out_366874660792644059) {
   out_366874660792644059[0] = state[4];
}
void H_30(double *state, double *unused, double *out_4649328563763672340) {
   out_4649328563763672340[0] = 0;
   out_4649328563763672340[1] = 0;
   out_4649328563763672340[2] = 0;
   out_4649328563763672340[3] = 0;
   out_4649328563763672340[4] = 1;
   out_4649328563763672340[5] = 0;
   out_4649328563763672340[6] = 0;
   out_4649328563763672340[7] = 0;
   out_4649328563763672340[8] = 0;
}
void h_26(double *state, double *unused, double *out_576635192109018376) {
   out_576635192109018376[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7537579232564574425) {
   out_7537579232564574425[0] = 0;
   out_7537579232564574425[1] = 0;
   out_7537579232564574425[2] = 0;
   out_7537579232564574425[3] = 0;
   out_7537579232564574425[4] = 0;
   out_7537579232564574425[5] = 0;
   out_7537579232564574425[6] = 0;
   out_7537579232564574425[7] = 1;
   out_7537579232564574425[8] = 0;
}
void h_27(double *state, double *unused, double *out_7644831426870215727) {
   out_7644831426870215727[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6824091875564097251) {
   out_6824091875564097251[0] = 0;
   out_6824091875564097251[1] = 0;
   out_6824091875564097251[2] = 0;
   out_6824091875564097251[3] = 1;
   out_6824091875564097251[4] = 0;
   out_6824091875564097251[5] = 0;
   out_6824091875564097251[6] = 0;
   out_6824091875564097251[7] = 0;
   out_6824091875564097251[8] = 0;
}
void h_29(double *state, double *unused, double *out_7802018095221344872) {
   out_7802018095221344872[0] = state[1];
}
void H_29(double *state, double *unused, double *out_4139097219449280156) {
   out_4139097219449280156[0] = 0;
   out_4139097219449280156[1] = 1;
   out_4139097219449280156[2] = 0;
   out_4139097219449280156[3] = 0;
   out_4139097219449280156[4] = 0;
   out_4139097219449280156[5] = 0;
   out_4139097219449280156[6] = 0;
   out_4139097219449280156[7] = 0;
   out_4139097219449280156[8] = 0;
}
void h_28(double *state, double *unused, double *out_8731558075226680054) {
   out_8731558075226680054[0] = state[0];
}
void H_28(double *state, double *unused, double *out_9221496236518810730) {
   out_9221496236518810730[0] = 1;
   out_9221496236518810730[1] = 0;
   out_9221496236518810730[2] = 0;
   out_9221496236518810730[3] = 0;
   out_9221496236518810730[4] = 0;
   out_9221496236518810730[5] = 0;
   out_9221496236518810730[6] = 0;
   out_9221496236518810730[7] = 0;
   out_9221496236518810730[8] = 0;
}
void h_31(double *state, double *unused, double *out_7250870771807373411) {
   out_7250870771807373411[0] = state[8];
}
void H_31(double *state, double *unused, double *out_6911371130331222949) {
   out_6911371130331222949[0] = 0;
   out_6911371130331222949[1] = 0;
   out_6911371130331222949[2] = 0;
   out_6911371130331222949[3] = 0;
   out_6911371130331222949[4] = 0;
   out_6911371130331222949[5] = 0;
   out_6911371130331222949[6] = 0;
   out_6911371130331222949[7] = 0;
   out_6911371130331222949[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_7193476806311499376) {
  err_fun(nom_x, delta_x, out_7193476806311499376);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8103533003814411399) {
  inv_err_fun(nom_x, true_x, out_8103533003814411399);
}
void car_H_mod_fun(double *state, double *out_4590973562022669193) {
  H_mod_fun(state, out_4590973562022669193);
}
void car_f_fun(double *state, double dt, double *out_738104619053175027) {
  f_fun(state,  dt, out_738104619053175027);
}
void car_F_fun(double *state, double dt, double *out_7442358059839394283) {
  F_fun(state,  dt, out_7442358059839394283);
}
void car_h_25(double *state, double *unused, double *out_5800239379178490244) {
  h_25(state, unused, out_5800239379178490244);
}
void car_H_25(double *state, double *unused, double *out_7167661522270920967) {
  H_25(state, unused, out_7167661522270920967);
}
void car_h_24(double *state, double *unused, double *out_3363033020538636619) {
  h_24(state, unused, out_3363033020538636619);
}
void car_H_24(double *state, double *unused, double *out_4990447098663770994) {
  H_24(state, unused, out_4990447098663770994);
}
void car_h_30(double *state, double *unused, double *out_366874660792644059) {
  h_30(state, unused, out_366874660792644059);
}
void car_H_30(double *state, double *unused, double *out_4649328563763672340) {
  H_30(state, unused, out_4649328563763672340);
}
void car_h_26(double *state, double *unused, double *out_576635192109018376) {
  h_26(state, unused, out_576635192109018376);
}
void car_H_26(double *state, double *unused, double *out_7537579232564574425) {
  H_26(state, unused, out_7537579232564574425);
}
void car_h_27(double *state, double *unused, double *out_7644831426870215727) {
  h_27(state, unused, out_7644831426870215727);
}
void car_H_27(double *state, double *unused, double *out_6824091875564097251) {
  H_27(state, unused, out_6824091875564097251);
}
void car_h_29(double *state, double *unused, double *out_7802018095221344872) {
  h_29(state, unused, out_7802018095221344872);
}
void car_H_29(double *state, double *unused, double *out_4139097219449280156) {
  H_29(state, unused, out_4139097219449280156);
}
void car_h_28(double *state, double *unused, double *out_8731558075226680054) {
  h_28(state, unused, out_8731558075226680054);
}
void car_H_28(double *state, double *unused, double *out_9221496236518810730) {
  H_28(state, unused, out_9221496236518810730);
}
void car_h_31(double *state, double *unused, double *out_7250870771807373411) {
  h_31(state, unused, out_7250870771807373411);
}
void car_H_31(double *state, double *unused, double *out_6911371130331222949) {
  H_31(state, unused, out_6911371130331222949);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
