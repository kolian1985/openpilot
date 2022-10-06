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
void err_fun(double *nom_x, double *delta_x, double *out_1077158114236301922) {
   out_1077158114236301922[0] = delta_x[0] + nom_x[0];
   out_1077158114236301922[1] = delta_x[1] + nom_x[1];
   out_1077158114236301922[2] = delta_x[2] + nom_x[2];
   out_1077158114236301922[3] = delta_x[3] + nom_x[3];
   out_1077158114236301922[4] = delta_x[4] + nom_x[4];
   out_1077158114236301922[5] = delta_x[5] + nom_x[5];
   out_1077158114236301922[6] = delta_x[6] + nom_x[6];
   out_1077158114236301922[7] = delta_x[7] + nom_x[7];
   out_1077158114236301922[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5038409095361470945) {
   out_5038409095361470945[0] = -nom_x[0] + true_x[0];
   out_5038409095361470945[1] = -nom_x[1] + true_x[1];
   out_5038409095361470945[2] = -nom_x[2] + true_x[2];
   out_5038409095361470945[3] = -nom_x[3] + true_x[3];
   out_5038409095361470945[4] = -nom_x[4] + true_x[4];
   out_5038409095361470945[5] = -nom_x[5] + true_x[5];
   out_5038409095361470945[6] = -nom_x[6] + true_x[6];
   out_5038409095361470945[7] = -nom_x[7] + true_x[7];
   out_5038409095361470945[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_3319339604940750727) {
   out_3319339604940750727[0] = 1.0;
   out_3319339604940750727[1] = 0;
   out_3319339604940750727[2] = 0;
   out_3319339604940750727[3] = 0;
   out_3319339604940750727[4] = 0;
   out_3319339604940750727[5] = 0;
   out_3319339604940750727[6] = 0;
   out_3319339604940750727[7] = 0;
   out_3319339604940750727[8] = 0;
   out_3319339604940750727[9] = 0;
   out_3319339604940750727[10] = 1.0;
   out_3319339604940750727[11] = 0;
   out_3319339604940750727[12] = 0;
   out_3319339604940750727[13] = 0;
   out_3319339604940750727[14] = 0;
   out_3319339604940750727[15] = 0;
   out_3319339604940750727[16] = 0;
   out_3319339604940750727[17] = 0;
   out_3319339604940750727[18] = 0;
   out_3319339604940750727[19] = 0;
   out_3319339604940750727[20] = 1.0;
   out_3319339604940750727[21] = 0;
   out_3319339604940750727[22] = 0;
   out_3319339604940750727[23] = 0;
   out_3319339604940750727[24] = 0;
   out_3319339604940750727[25] = 0;
   out_3319339604940750727[26] = 0;
   out_3319339604940750727[27] = 0;
   out_3319339604940750727[28] = 0;
   out_3319339604940750727[29] = 0;
   out_3319339604940750727[30] = 1.0;
   out_3319339604940750727[31] = 0;
   out_3319339604940750727[32] = 0;
   out_3319339604940750727[33] = 0;
   out_3319339604940750727[34] = 0;
   out_3319339604940750727[35] = 0;
   out_3319339604940750727[36] = 0;
   out_3319339604940750727[37] = 0;
   out_3319339604940750727[38] = 0;
   out_3319339604940750727[39] = 0;
   out_3319339604940750727[40] = 1.0;
   out_3319339604940750727[41] = 0;
   out_3319339604940750727[42] = 0;
   out_3319339604940750727[43] = 0;
   out_3319339604940750727[44] = 0;
   out_3319339604940750727[45] = 0;
   out_3319339604940750727[46] = 0;
   out_3319339604940750727[47] = 0;
   out_3319339604940750727[48] = 0;
   out_3319339604940750727[49] = 0;
   out_3319339604940750727[50] = 1.0;
   out_3319339604940750727[51] = 0;
   out_3319339604940750727[52] = 0;
   out_3319339604940750727[53] = 0;
   out_3319339604940750727[54] = 0;
   out_3319339604940750727[55] = 0;
   out_3319339604940750727[56] = 0;
   out_3319339604940750727[57] = 0;
   out_3319339604940750727[58] = 0;
   out_3319339604940750727[59] = 0;
   out_3319339604940750727[60] = 1.0;
   out_3319339604940750727[61] = 0;
   out_3319339604940750727[62] = 0;
   out_3319339604940750727[63] = 0;
   out_3319339604940750727[64] = 0;
   out_3319339604940750727[65] = 0;
   out_3319339604940750727[66] = 0;
   out_3319339604940750727[67] = 0;
   out_3319339604940750727[68] = 0;
   out_3319339604940750727[69] = 0;
   out_3319339604940750727[70] = 1.0;
   out_3319339604940750727[71] = 0;
   out_3319339604940750727[72] = 0;
   out_3319339604940750727[73] = 0;
   out_3319339604940750727[74] = 0;
   out_3319339604940750727[75] = 0;
   out_3319339604940750727[76] = 0;
   out_3319339604940750727[77] = 0;
   out_3319339604940750727[78] = 0;
   out_3319339604940750727[79] = 0;
   out_3319339604940750727[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_1751646311246354893) {
   out_1751646311246354893[0] = state[0];
   out_1751646311246354893[1] = state[1];
   out_1751646311246354893[2] = state[2];
   out_1751646311246354893[3] = state[3];
   out_1751646311246354893[4] = state[4];
   out_1751646311246354893[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_1751646311246354893[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_1751646311246354893[7] = state[7];
   out_1751646311246354893[8] = state[8];
}
void F_fun(double *state, double dt, double *out_2896000726666627884) {
   out_2896000726666627884[0] = 1;
   out_2896000726666627884[1] = 0;
   out_2896000726666627884[2] = 0;
   out_2896000726666627884[3] = 0;
   out_2896000726666627884[4] = 0;
   out_2896000726666627884[5] = 0;
   out_2896000726666627884[6] = 0;
   out_2896000726666627884[7] = 0;
   out_2896000726666627884[8] = 0;
   out_2896000726666627884[9] = 0;
   out_2896000726666627884[10] = 1;
   out_2896000726666627884[11] = 0;
   out_2896000726666627884[12] = 0;
   out_2896000726666627884[13] = 0;
   out_2896000726666627884[14] = 0;
   out_2896000726666627884[15] = 0;
   out_2896000726666627884[16] = 0;
   out_2896000726666627884[17] = 0;
   out_2896000726666627884[18] = 0;
   out_2896000726666627884[19] = 0;
   out_2896000726666627884[20] = 1;
   out_2896000726666627884[21] = 0;
   out_2896000726666627884[22] = 0;
   out_2896000726666627884[23] = 0;
   out_2896000726666627884[24] = 0;
   out_2896000726666627884[25] = 0;
   out_2896000726666627884[26] = 0;
   out_2896000726666627884[27] = 0;
   out_2896000726666627884[28] = 0;
   out_2896000726666627884[29] = 0;
   out_2896000726666627884[30] = 1;
   out_2896000726666627884[31] = 0;
   out_2896000726666627884[32] = 0;
   out_2896000726666627884[33] = 0;
   out_2896000726666627884[34] = 0;
   out_2896000726666627884[35] = 0;
   out_2896000726666627884[36] = 0;
   out_2896000726666627884[37] = 0;
   out_2896000726666627884[38] = 0;
   out_2896000726666627884[39] = 0;
   out_2896000726666627884[40] = 1;
   out_2896000726666627884[41] = 0;
   out_2896000726666627884[42] = 0;
   out_2896000726666627884[43] = 0;
   out_2896000726666627884[44] = 0;
   out_2896000726666627884[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2896000726666627884[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2896000726666627884[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2896000726666627884[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2896000726666627884[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2896000726666627884[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2896000726666627884[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2896000726666627884[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2896000726666627884[53] = -9.8000000000000007*dt;
   out_2896000726666627884[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2896000726666627884[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2896000726666627884[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2896000726666627884[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2896000726666627884[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2896000726666627884[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2896000726666627884[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2896000726666627884[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2896000726666627884[62] = 0;
   out_2896000726666627884[63] = 0;
   out_2896000726666627884[64] = 0;
   out_2896000726666627884[65] = 0;
   out_2896000726666627884[66] = 0;
   out_2896000726666627884[67] = 0;
   out_2896000726666627884[68] = 0;
   out_2896000726666627884[69] = 0;
   out_2896000726666627884[70] = 1;
   out_2896000726666627884[71] = 0;
   out_2896000726666627884[72] = 0;
   out_2896000726666627884[73] = 0;
   out_2896000726666627884[74] = 0;
   out_2896000726666627884[75] = 0;
   out_2896000726666627884[76] = 0;
   out_2896000726666627884[77] = 0;
   out_2896000726666627884[78] = 0;
   out_2896000726666627884[79] = 0;
   out_2896000726666627884[80] = 1;
}
void h_25(double *state, double *unused, double *out_6359605128438556057) {
   out_6359605128438556057[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2082710691673890726) {
   out_2082710691673890726[0] = 0;
   out_2082710691673890726[1] = 0;
   out_2082710691673890726[2] = 0;
   out_2082710691673890726[3] = 0;
   out_2082710691673890726[4] = 0;
   out_2082710691673890726[5] = 0;
   out_2082710691673890726[6] = 1;
   out_2082710691673890726[7] = 0;
   out_2082710691673890726[8] = 0;
}
void h_24(double *state, double *unused, double *out_5470019761866181810) {
   out_5470019761866181810[0] = state[4];
   out_5470019761866181810[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4492861114917627375) {
   out_4492861114917627375[0] = 0;
   out_4492861114917627375[1] = 0;
   out_4492861114917627375[2] = 0;
   out_4492861114917627375[3] = 0;
   out_4492861114917627375[4] = 1;
   out_4492861114917627375[5] = 0;
   out_4492861114917627375[6] = 0;
   out_4492861114917627375[7] = 0;
   out_4492861114917627375[8] = 0;
   out_4492861114917627375[9] = 0;
   out_4492861114917627375[10] = 0;
   out_4492861114917627375[11] = 0;
   out_4492861114917627375[12] = 0;
   out_4492861114917627375[13] = 0;
   out_4492861114917627375[14] = 1;
   out_4492861114917627375[15] = 0;
   out_4492861114917627375[16] = 0;
   out_4492861114917627375[17] = 0;
}
void h_30(double *state, double *unused, double *out_2722155600184533707) {
   out_2722155600184533707[0] = state[4];
}
void H_30(double *state, double *unused, double *out_4833979649817726029) {
   out_4833979649817726029[0] = 0;
   out_4833979649817726029[1] = 0;
   out_4833979649817726029[2] = 0;
   out_4833979649817726029[3] = 0;
   out_4833979649817726029[4] = 1;
   out_4833979649817726029[5] = 0;
   out_4833979649817726029[6] = 0;
   out_4833979649817726029[7] = 0;
   out_4833979649817726029[8] = 0;
}
void h_26(double *state, double *unused, double *out_3480005745345352517) {
   out_3480005745345352517[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5824214010547946950) {
   out_5824214010547946950[0] = 0;
   out_5824214010547946950[1] = 0;
   out_5824214010547946950[2] = 0;
   out_5824214010547946950[3] = 0;
   out_5824214010547946950[4] = 0;
   out_5824214010547946950[5] = 0;
   out_5824214010547946950[6] = 0;
   out_5824214010547946950[7] = 1;
   out_5824214010547946950[8] = 0;
}
void h_27(double *state, double *unused, double *out_9050352893826662874) {
   out_9050352893826662874[0] = state[3];
}
void H_27(double *state, double *unused, double *out_2659216338017301118) {
   out_2659216338017301118[0] = 0;
   out_2659216338017301118[1] = 0;
   out_2659216338017301118[2] = 0;
   out_2659216338017301118[3] = 1;
   out_2659216338017301118[4] = 0;
   out_2659216338017301118[5] = 0;
   out_2659216338017301118[6] = 0;
   out_2659216338017301118[7] = 0;
   out_2659216338017301118[8] = 0;
}
void h_29(double *state, double *unused, double *out_2279517667476311938) {
   out_2279517667476311938[0] = state[1];
}
void H_29(double *state, double *unused, double *out_5344210994132118213) {
   out_5344210994132118213[0] = 0;
   out_5344210994132118213[1] = 1;
   out_5344210994132118213[2] = 0;
   out_5344210994132118213[3] = 0;
   out_5344210994132118213[4] = 0;
   out_5344210994132118213[5] = 0;
   out_5344210994132118213[6] = 0;
   out_5344210994132118213[7] = 0;
   out_5344210994132118213[8] = 0;
}
void h_28(double *state, double *unused, double *out_3708485657940632311) {
   out_3708485657940632311[0] = state[0];
}
void H_28(double *state, double *unused, double *out_4136545405921780489) {
   out_4136545405921780489[0] = 1;
   out_4136545405921780489[1] = 0;
   out_4136545405921780489[2] = 0;
   out_4136545405921780489[3] = 0;
   out_4136545405921780489[4] = 0;
   out_4136545405921780489[5] = 0;
   out_4136545405921780489[6] = 0;
   out_4136545405921780489[7] = 0;
   out_4136545405921780489[8] = 0;
}
void h_31(double *state, double *unused, double *out_1100953635810332268) {
   out_1100953635810332268[0] = state[8];
}
void H_31(double *state, double *unused, double *out_2052064729796930298) {
   out_2052064729796930298[0] = 0;
   out_2052064729796930298[1] = 0;
   out_2052064729796930298[2] = 0;
   out_2052064729796930298[3] = 0;
   out_2052064729796930298[4] = 0;
   out_2052064729796930298[5] = 0;
   out_2052064729796930298[6] = 0;
   out_2052064729796930298[7] = 0;
   out_2052064729796930298[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_1077158114236301922) {
  err_fun(nom_x, delta_x, out_1077158114236301922);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5038409095361470945) {
  inv_err_fun(nom_x, true_x, out_5038409095361470945);
}
void car_H_mod_fun(double *state, double *out_3319339604940750727) {
  H_mod_fun(state, out_3319339604940750727);
}
void car_f_fun(double *state, double dt, double *out_1751646311246354893) {
  f_fun(state,  dt, out_1751646311246354893);
}
void car_F_fun(double *state, double dt, double *out_2896000726666627884) {
  F_fun(state,  dt, out_2896000726666627884);
}
void car_h_25(double *state, double *unused, double *out_6359605128438556057) {
  h_25(state, unused, out_6359605128438556057);
}
void car_H_25(double *state, double *unused, double *out_2082710691673890726) {
  H_25(state, unused, out_2082710691673890726);
}
void car_h_24(double *state, double *unused, double *out_5470019761866181810) {
  h_24(state, unused, out_5470019761866181810);
}
void car_H_24(double *state, double *unused, double *out_4492861114917627375) {
  H_24(state, unused, out_4492861114917627375);
}
void car_h_30(double *state, double *unused, double *out_2722155600184533707) {
  h_30(state, unused, out_2722155600184533707);
}
void car_H_30(double *state, double *unused, double *out_4833979649817726029) {
  H_30(state, unused, out_4833979649817726029);
}
void car_h_26(double *state, double *unused, double *out_3480005745345352517) {
  h_26(state, unused, out_3480005745345352517);
}
void car_H_26(double *state, double *unused, double *out_5824214010547946950) {
  H_26(state, unused, out_5824214010547946950);
}
void car_h_27(double *state, double *unused, double *out_9050352893826662874) {
  h_27(state, unused, out_9050352893826662874);
}
void car_H_27(double *state, double *unused, double *out_2659216338017301118) {
  H_27(state, unused, out_2659216338017301118);
}
void car_h_29(double *state, double *unused, double *out_2279517667476311938) {
  h_29(state, unused, out_2279517667476311938);
}
void car_H_29(double *state, double *unused, double *out_5344210994132118213) {
  H_29(state, unused, out_5344210994132118213);
}
void car_h_28(double *state, double *unused, double *out_3708485657940632311) {
  h_28(state, unused, out_3708485657940632311);
}
void car_H_28(double *state, double *unused, double *out_4136545405921780489) {
  H_28(state, unused, out_4136545405921780489);
}
void car_h_31(double *state, double *unused, double *out_1100953635810332268) {
  h_31(state, unused, out_1100953635810332268);
}
void car_H_31(double *state, double *unused, double *out_2052064729796930298) {
  H_31(state, unused, out_2052064729796930298);
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
