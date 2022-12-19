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
void err_fun(double *nom_x, double *delta_x, double *out_4847323219916667687) {
   out_4847323219916667687[0] = delta_x[0] + nom_x[0];
   out_4847323219916667687[1] = delta_x[1] + nom_x[1];
   out_4847323219916667687[2] = delta_x[2] + nom_x[2];
   out_4847323219916667687[3] = delta_x[3] + nom_x[3];
   out_4847323219916667687[4] = delta_x[4] + nom_x[4];
   out_4847323219916667687[5] = delta_x[5] + nom_x[5];
   out_4847323219916667687[6] = delta_x[6] + nom_x[6];
   out_4847323219916667687[7] = delta_x[7] + nom_x[7];
   out_4847323219916667687[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_6388774152708464872) {
   out_6388774152708464872[0] = -nom_x[0] + true_x[0];
   out_6388774152708464872[1] = -nom_x[1] + true_x[1];
   out_6388774152708464872[2] = -nom_x[2] + true_x[2];
   out_6388774152708464872[3] = -nom_x[3] + true_x[3];
   out_6388774152708464872[4] = -nom_x[4] + true_x[4];
   out_6388774152708464872[5] = -nom_x[5] + true_x[5];
   out_6388774152708464872[6] = -nom_x[6] + true_x[6];
   out_6388774152708464872[7] = -nom_x[7] + true_x[7];
   out_6388774152708464872[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_2614378576525301819) {
   out_2614378576525301819[0] = 1.0;
   out_2614378576525301819[1] = 0;
   out_2614378576525301819[2] = 0;
   out_2614378576525301819[3] = 0;
   out_2614378576525301819[4] = 0;
   out_2614378576525301819[5] = 0;
   out_2614378576525301819[6] = 0;
   out_2614378576525301819[7] = 0;
   out_2614378576525301819[8] = 0;
   out_2614378576525301819[9] = 0;
   out_2614378576525301819[10] = 1.0;
   out_2614378576525301819[11] = 0;
   out_2614378576525301819[12] = 0;
   out_2614378576525301819[13] = 0;
   out_2614378576525301819[14] = 0;
   out_2614378576525301819[15] = 0;
   out_2614378576525301819[16] = 0;
   out_2614378576525301819[17] = 0;
   out_2614378576525301819[18] = 0;
   out_2614378576525301819[19] = 0;
   out_2614378576525301819[20] = 1.0;
   out_2614378576525301819[21] = 0;
   out_2614378576525301819[22] = 0;
   out_2614378576525301819[23] = 0;
   out_2614378576525301819[24] = 0;
   out_2614378576525301819[25] = 0;
   out_2614378576525301819[26] = 0;
   out_2614378576525301819[27] = 0;
   out_2614378576525301819[28] = 0;
   out_2614378576525301819[29] = 0;
   out_2614378576525301819[30] = 1.0;
   out_2614378576525301819[31] = 0;
   out_2614378576525301819[32] = 0;
   out_2614378576525301819[33] = 0;
   out_2614378576525301819[34] = 0;
   out_2614378576525301819[35] = 0;
   out_2614378576525301819[36] = 0;
   out_2614378576525301819[37] = 0;
   out_2614378576525301819[38] = 0;
   out_2614378576525301819[39] = 0;
   out_2614378576525301819[40] = 1.0;
   out_2614378576525301819[41] = 0;
   out_2614378576525301819[42] = 0;
   out_2614378576525301819[43] = 0;
   out_2614378576525301819[44] = 0;
   out_2614378576525301819[45] = 0;
   out_2614378576525301819[46] = 0;
   out_2614378576525301819[47] = 0;
   out_2614378576525301819[48] = 0;
   out_2614378576525301819[49] = 0;
   out_2614378576525301819[50] = 1.0;
   out_2614378576525301819[51] = 0;
   out_2614378576525301819[52] = 0;
   out_2614378576525301819[53] = 0;
   out_2614378576525301819[54] = 0;
   out_2614378576525301819[55] = 0;
   out_2614378576525301819[56] = 0;
   out_2614378576525301819[57] = 0;
   out_2614378576525301819[58] = 0;
   out_2614378576525301819[59] = 0;
   out_2614378576525301819[60] = 1.0;
   out_2614378576525301819[61] = 0;
   out_2614378576525301819[62] = 0;
   out_2614378576525301819[63] = 0;
   out_2614378576525301819[64] = 0;
   out_2614378576525301819[65] = 0;
   out_2614378576525301819[66] = 0;
   out_2614378576525301819[67] = 0;
   out_2614378576525301819[68] = 0;
   out_2614378576525301819[69] = 0;
   out_2614378576525301819[70] = 1.0;
   out_2614378576525301819[71] = 0;
   out_2614378576525301819[72] = 0;
   out_2614378576525301819[73] = 0;
   out_2614378576525301819[74] = 0;
   out_2614378576525301819[75] = 0;
   out_2614378576525301819[76] = 0;
   out_2614378576525301819[77] = 0;
   out_2614378576525301819[78] = 0;
   out_2614378576525301819[79] = 0;
   out_2614378576525301819[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_6166246757365277964) {
   out_6166246757365277964[0] = state[0];
   out_6166246757365277964[1] = state[1];
   out_6166246757365277964[2] = state[2];
   out_6166246757365277964[3] = state[3];
   out_6166246757365277964[4] = state[4];
   out_6166246757365277964[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6166246757365277964[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6166246757365277964[7] = state[7];
   out_6166246757365277964[8] = state[8];
}
void F_fun(double *state, double dt, double *out_6257381623778681963) {
   out_6257381623778681963[0] = 1;
   out_6257381623778681963[1] = 0;
   out_6257381623778681963[2] = 0;
   out_6257381623778681963[3] = 0;
   out_6257381623778681963[4] = 0;
   out_6257381623778681963[5] = 0;
   out_6257381623778681963[6] = 0;
   out_6257381623778681963[7] = 0;
   out_6257381623778681963[8] = 0;
   out_6257381623778681963[9] = 0;
   out_6257381623778681963[10] = 1;
   out_6257381623778681963[11] = 0;
   out_6257381623778681963[12] = 0;
   out_6257381623778681963[13] = 0;
   out_6257381623778681963[14] = 0;
   out_6257381623778681963[15] = 0;
   out_6257381623778681963[16] = 0;
   out_6257381623778681963[17] = 0;
   out_6257381623778681963[18] = 0;
   out_6257381623778681963[19] = 0;
   out_6257381623778681963[20] = 1;
   out_6257381623778681963[21] = 0;
   out_6257381623778681963[22] = 0;
   out_6257381623778681963[23] = 0;
   out_6257381623778681963[24] = 0;
   out_6257381623778681963[25] = 0;
   out_6257381623778681963[26] = 0;
   out_6257381623778681963[27] = 0;
   out_6257381623778681963[28] = 0;
   out_6257381623778681963[29] = 0;
   out_6257381623778681963[30] = 1;
   out_6257381623778681963[31] = 0;
   out_6257381623778681963[32] = 0;
   out_6257381623778681963[33] = 0;
   out_6257381623778681963[34] = 0;
   out_6257381623778681963[35] = 0;
   out_6257381623778681963[36] = 0;
   out_6257381623778681963[37] = 0;
   out_6257381623778681963[38] = 0;
   out_6257381623778681963[39] = 0;
   out_6257381623778681963[40] = 1;
   out_6257381623778681963[41] = 0;
   out_6257381623778681963[42] = 0;
   out_6257381623778681963[43] = 0;
   out_6257381623778681963[44] = 0;
   out_6257381623778681963[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_6257381623778681963[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_6257381623778681963[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6257381623778681963[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6257381623778681963[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_6257381623778681963[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_6257381623778681963[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_6257381623778681963[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_6257381623778681963[53] = -9.8000000000000007*dt;
   out_6257381623778681963[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_6257381623778681963[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_6257381623778681963[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6257381623778681963[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6257381623778681963[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_6257381623778681963[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_6257381623778681963[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_6257381623778681963[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6257381623778681963[62] = 0;
   out_6257381623778681963[63] = 0;
   out_6257381623778681963[64] = 0;
   out_6257381623778681963[65] = 0;
   out_6257381623778681963[66] = 0;
   out_6257381623778681963[67] = 0;
   out_6257381623778681963[68] = 0;
   out_6257381623778681963[69] = 0;
   out_6257381623778681963[70] = 1;
   out_6257381623778681963[71] = 0;
   out_6257381623778681963[72] = 0;
   out_6257381623778681963[73] = 0;
   out_6257381623778681963[74] = 0;
   out_6257381623778681963[75] = 0;
   out_6257381623778681963[76] = 0;
   out_6257381623778681963[77] = 0;
   out_6257381623778681963[78] = 0;
   out_6257381623778681963[79] = 0;
   out_6257381623778681963[80] = 1;
}
void h_25(double *state, double *unused, double *out_9107992015900631963) {
   out_9107992015900631963[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3138229866802375607) {
   out_3138229866802375607[0] = 0;
   out_3138229866802375607[1] = 0;
   out_3138229866802375607[2] = 0;
   out_3138229866802375607[3] = 0;
   out_3138229866802375607[4] = 0;
   out_3138229866802375607[5] = 0;
   out_3138229866802375607[6] = 1;
   out_3138229866802375607[7] = 0;
   out_3138229866802375607[8] = 0;
}
void h_24(double *state, double *unused, double *out_736457797204437392) {
   out_736457797204437392[0] = state[4];
   out_736457797204437392[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6080449020837980784) {
   out_6080449020837980784[0] = 0;
   out_6080449020837980784[1] = 0;
   out_6080449020837980784[2] = 0;
   out_6080449020837980784[3] = 0;
   out_6080449020837980784[4] = 1;
   out_6080449020837980784[5] = 0;
   out_6080449020837980784[6] = 0;
   out_6080449020837980784[7] = 0;
   out_6080449020837980784[8] = 0;
   out_6080449020837980784[9] = 0;
   out_6080449020837980784[10] = 0;
   out_6080449020837980784[11] = 0;
   out_6080449020837980784[12] = 0;
   out_6080449020837980784[13] = 0;
   out_6080449020837980784[14] = 1;
   out_6080449020837980784[15] = 0;
   out_6080449020837980784[16] = 0;
   out_6080449020837980784[17] = 0;
}
void h_30(double *state, double *unused, double *out_703670155659848332) {
   out_703670155659848332[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1389466463325232591) {
   out_1389466463325232591[0] = 0;
   out_1389466463325232591[1] = 0;
   out_1389466463325232591[2] = 0;
   out_1389466463325232591[3] = 0;
   out_1389466463325232591[4] = 1;
   out_1389466463325232591[5] = 0;
   out_1389466463325232591[6] = 0;
   out_1389466463325232591[7] = 0;
   out_1389466463325232591[8] = 0;
}
void h_26(double *state, double *unused, double *out_3740456207104181017) {
   out_3740456207104181017[0] = state[7];
}
void H_26(double *state, double *unused, double *out_603273452071680617) {
   out_603273452071680617[0] = 0;
   out_603273452071680617[1] = 0;
   out_603273452071680617[2] = 0;
   out_603273452071680617[3] = 0;
   out_603273452071680617[4] = 0;
   out_603273452071680617[5] = 0;
   out_603273452071680617[6] = 0;
   out_603273452071680617[7] = 1;
   out_603273452071680617[8] = 0;
}
void h_27(double *state, double *unused, double *out_515033759526319992) {
   out_515033759526319992[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3564229775125657502) {
   out_3564229775125657502[0] = 0;
   out_3564229775125657502[1] = 0;
   out_3564229775125657502[2] = 0;
   out_3564229775125657502[3] = 1;
   out_3564229775125657502[4] = 0;
   out_3564229775125657502[5] = 0;
   out_3564229775125657502[6] = 0;
   out_3564229775125657502[7] = 0;
   out_3564229775125657502[8] = 0;
}
void h_29(double *state, double *unused, double *out_239839697241814103) {
   out_239839697241814103[0] = state[1];
}
void H_29(double *state, double *unused, double *out_879235119010840407) {
   out_879235119010840407[0] = 0;
   out_879235119010840407[1] = 1;
   out_879235119010840407[2] = 0;
   out_879235119010840407[3] = 0;
   out_879235119010840407[4] = 0;
   out_879235119010840407[5] = 0;
   out_879235119010840407[6] = 0;
   out_879235119010840407[7] = 0;
   out_879235119010840407[8] = 0;
}
void h_28(double *state, double *unused, double *out_7517796463319385771) {
   out_7517796463319385771[0] = state[0];
}
void H_28(double *state, double *unused, double *out_1084395152554485844) {
   out_1084395152554485844[0] = 1;
   out_1084395152554485844[1] = 0;
   out_1084395152554485844[2] = 0;
   out_1084395152554485844[3] = 0;
   out_1084395152554485844[4] = 0;
   out_1084395152554485844[5] = 0;
   out_1084395152554485844[6] = 0;
   out_1084395152554485844[7] = 0;
   out_1084395152554485844[8] = 0;
}
void h_31(double *state, double *unused, double *out_8464433017542650598) {
   out_8464433017542650598[0] = state[8];
}
void H_31(double *state, double *unused, double *out_1229481554305032093) {
   out_1229481554305032093[0] = 0;
   out_1229481554305032093[1] = 0;
   out_1229481554305032093[2] = 0;
   out_1229481554305032093[3] = 0;
   out_1229481554305032093[4] = 0;
   out_1229481554305032093[5] = 0;
   out_1229481554305032093[6] = 0;
   out_1229481554305032093[7] = 0;
   out_1229481554305032093[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_4847323219916667687) {
  err_fun(nom_x, delta_x, out_4847323219916667687);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6388774152708464872) {
  inv_err_fun(nom_x, true_x, out_6388774152708464872);
}
void car_H_mod_fun(double *state, double *out_2614378576525301819) {
  H_mod_fun(state, out_2614378576525301819);
}
void car_f_fun(double *state, double dt, double *out_6166246757365277964) {
  f_fun(state,  dt, out_6166246757365277964);
}
void car_F_fun(double *state, double dt, double *out_6257381623778681963) {
  F_fun(state,  dt, out_6257381623778681963);
}
void car_h_25(double *state, double *unused, double *out_9107992015900631963) {
  h_25(state, unused, out_9107992015900631963);
}
void car_H_25(double *state, double *unused, double *out_3138229866802375607) {
  H_25(state, unused, out_3138229866802375607);
}
void car_h_24(double *state, double *unused, double *out_736457797204437392) {
  h_24(state, unused, out_736457797204437392);
}
void car_H_24(double *state, double *unused, double *out_6080449020837980784) {
  H_24(state, unused, out_6080449020837980784);
}
void car_h_30(double *state, double *unused, double *out_703670155659848332) {
  h_30(state, unused, out_703670155659848332);
}
void car_H_30(double *state, double *unused, double *out_1389466463325232591) {
  H_30(state, unused, out_1389466463325232591);
}
void car_h_26(double *state, double *unused, double *out_3740456207104181017) {
  h_26(state, unused, out_3740456207104181017);
}
void car_H_26(double *state, double *unused, double *out_603273452071680617) {
  H_26(state, unused, out_603273452071680617);
}
void car_h_27(double *state, double *unused, double *out_515033759526319992) {
  h_27(state, unused, out_515033759526319992);
}
void car_H_27(double *state, double *unused, double *out_3564229775125657502) {
  H_27(state, unused, out_3564229775125657502);
}
void car_h_29(double *state, double *unused, double *out_239839697241814103) {
  h_29(state, unused, out_239839697241814103);
}
void car_H_29(double *state, double *unused, double *out_879235119010840407) {
  H_29(state, unused, out_879235119010840407);
}
void car_h_28(double *state, double *unused, double *out_7517796463319385771) {
  h_28(state, unused, out_7517796463319385771);
}
void car_H_28(double *state, double *unused, double *out_1084395152554485844) {
  H_28(state, unused, out_1084395152554485844);
}
void car_h_31(double *state, double *unused, double *out_8464433017542650598) {
  h_31(state, unused, out_8464433017542650598);
}
void car_H_31(double *state, double *unused, double *out_1229481554305032093) {
  H_31(state, unused, out_1229481554305032093);
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
