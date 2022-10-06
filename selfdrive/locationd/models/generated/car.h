#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_1077158114236301922);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5038409095361470945);
void car_H_mod_fun(double *state, double *out_3319339604940750727);
void car_f_fun(double *state, double dt, double *out_1751646311246354893);
void car_F_fun(double *state, double dt, double *out_2896000726666627884);
void car_h_25(double *state, double *unused, double *out_6359605128438556057);
void car_H_25(double *state, double *unused, double *out_2082710691673890726);
void car_h_24(double *state, double *unused, double *out_5470019761866181810);
void car_H_24(double *state, double *unused, double *out_4492861114917627375);
void car_h_30(double *state, double *unused, double *out_2722155600184533707);
void car_H_30(double *state, double *unused, double *out_4833979649817726029);
void car_h_26(double *state, double *unused, double *out_3480005745345352517);
void car_H_26(double *state, double *unused, double *out_5824214010547946950);
void car_h_27(double *state, double *unused, double *out_9050352893826662874);
void car_H_27(double *state, double *unused, double *out_2659216338017301118);
void car_h_29(double *state, double *unused, double *out_2279517667476311938);
void car_H_29(double *state, double *unused, double *out_5344210994132118213);
void car_h_28(double *state, double *unused, double *out_3708485657940632311);
void car_H_28(double *state, double *unused, double *out_4136545405921780489);
void car_h_31(double *state, double *unused, double *out_1100953635810332268);
void car_H_31(double *state, double *unused, double *out_2052064729796930298);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}