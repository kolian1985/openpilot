#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_6759259717540961999);
void live_err_fun(double *nom_x, double *delta_x, double *out_712671743316778553);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_1551098332499718356);
void live_H_mod_fun(double *state, double *out_7646199264428369681);
void live_f_fun(double *state, double dt, double *out_2537679088963151294);
void live_F_fun(double *state, double dt, double *out_2272821566001163028);
void live_h_4(double *state, double *unused, double *out_3772097625145546202);
void live_H_4(double *state, double *unused, double *out_3391770475938062923);
void live_h_9(double *state, double *unused, double *out_1961205282106539496);
void live_H_9(double *state, double *unused, double *out_502908923657983581);
void live_h_10(double *state, double *unused, double *out_1566763082912586440);
void live_H_10(double *state, double *unused, double *out_7997271869951827050);
void live_h_12(double *state, double *unused, double *out_3013028664800893947);
void live_H_12(double *state, double *unused, double *out_2770671450890469256);
void live_h_31(double *state, double *unused, double *out_7147559448430490051);
void live_H_31(double *state, double *unused, double *out_25108418565455547);
void live_h_32(double *state, double *unused, double *out_351136401439433302);
void live_H_32(double *state, double *unused, double *out_5776195953642157011);
void live_h_13(double *state, double *unused, double *out_3516376497702871266);
void live_H_13(double *state, double *unused, double *out_3583990703157424414);
void live_h_14(double *state, double *unused, double *out_1961205282106539496);
void live_H_14(double *state, double *unused, double *out_502908923657983581);
void live_h_33(double *state, double *unused, double *out_1044032644426440076);
void live_H_33(double *state, double *unused, double *out_3125448586073402057);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}