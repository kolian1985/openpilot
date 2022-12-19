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
void car_err_fun(double *nom_x, double *delta_x, double *out_4847323219916667687);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6388774152708464872);
void car_H_mod_fun(double *state, double *out_2614378576525301819);
void car_f_fun(double *state, double dt, double *out_6166246757365277964);
void car_F_fun(double *state, double dt, double *out_6257381623778681963);
void car_h_25(double *state, double *unused, double *out_9107992015900631963);
void car_H_25(double *state, double *unused, double *out_3138229866802375607);
void car_h_24(double *state, double *unused, double *out_736457797204437392);
void car_H_24(double *state, double *unused, double *out_6080449020837980784);
void car_h_30(double *state, double *unused, double *out_703670155659848332);
void car_H_30(double *state, double *unused, double *out_1389466463325232591);
void car_h_26(double *state, double *unused, double *out_3740456207104181017);
void car_H_26(double *state, double *unused, double *out_603273452071680617);
void car_h_27(double *state, double *unused, double *out_515033759526319992);
void car_H_27(double *state, double *unused, double *out_3564229775125657502);
void car_h_29(double *state, double *unused, double *out_239839697241814103);
void car_H_29(double *state, double *unused, double *out_879235119010840407);
void car_h_28(double *state, double *unused, double *out_7517796463319385771);
void car_H_28(double *state, double *unused, double *out_1084395152554485844);
void car_h_31(double *state, double *unused, double *out_8464433017542650598);
void car_H_31(double *state, double *unused, double *out_1229481554305032093);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}