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
void car_err_fun(double *nom_x, double *delta_x, double *out_985208881199282010);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_2500588504390326563);
void car_H_mod_fun(double *state, double *out_8329055924251883713);
void car_f_fun(double *state, double dt, double *out_3209088959410559628);
void car_F_fun(double *state, double dt, double *out_1203173105859445993);
void car_h_25(double *state, double *unused, double *out_6192974400242693888);
void car_H_25(double *state, double *unused, double *out_8405907936316016165);
void car_h_24(double *state, double *unused, double *out_2380021089705406482);
void car_H_24(double *state, double *unused, double *out_7863621713786385478);
void car_h_30(double *state, double *unused, double *out_4360830178582334369);
void car_H_30(double *state, double *unused, double *out_7522503178886286824);
void car_h_26(double *state, double *unused, double *out_1084982365834877780);
void car_H_26(double *state, double *unused, double *out_4664404617441959941);
void car_h_27(double *state, double *unused, double *out_2630743898039905773);
void car_H_27(double *state, double *unused, double *out_5298909107702343607);
void car_h_29(double *state, double *unused, double *out_2670382343846653409);
void car_H_29(double *state, double *unused, double *out_7012271834571894640);
void car_h_28(double *state, double *unused, double *out_8273473769120862608);
void car_H_28(double *state, double *unused, double *out_6352073222068126402);
void car_h_31(double *state, double *unused, double *out_6468168462527199777);
void car_H_31(double *state, double *unused, double *out_4038196515208608465);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}