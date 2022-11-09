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
void car_err_fun(double *nom_x, double *delta_x, double *out_7193476806311499376);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8103533003814411399);
void car_H_mod_fun(double *state, double *out_4590973562022669193);
void car_f_fun(double *state, double dt, double *out_738104619053175027);
void car_F_fun(double *state, double dt, double *out_7442358059839394283);
void car_h_25(double *state, double *unused, double *out_5800239379178490244);
void car_H_25(double *state, double *unused, double *out_7167661522270920967);
void car_h_24(double *state, double *unused, double *out_3363033020538636619);
void car_H_24(double *state, double *unused, double *out_4990447098663770994);
void car_h_30(double *state, double *unused, double *out_366874660792644059);
void car_H_30(double *state, double *unused, double *out_4649328563763672340);
void car_h_26(double *state, double *unused, double *out_576635192109018376);
void car_H_26(double *state, double *unused, double *out_7537579232564574425);
void car_h_27(double *state, double *unused, double *out_7644831426870215727);
void car_H_27(double *state, double *unused, double *out_6824091875564097251);
void car_h_29(double *state, double *unused, double *out_7802018095221344872);
void car_H_29(double *state, double *unused, double *out_4139097219449280156);
void car_h_28(double *state, double *unused, double *out_8731558075226680054);
void car_H_28(double *state, double *unused, double *out_9221496236518810730);
void car_h_31(double *state, double *unused, double *out_7250870771807373411);
void car_H_31(double *state, double *unused, double *out_6911371130331222949);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}