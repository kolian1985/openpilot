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
void live_H(double *in_vec, double *out_6339076308092829214);
void live_err_fun(double *nom_x, double *delta_x, double *out_3402840131243625295);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_7045607237347834889);
void live_H_mod_fun(double *state, double *out_1201729222388576131);
void live_f_fun(double *state, double dt, double *out_4166071926154386066);
void live_F_fun(double *state, double dt, double *out_2554711962185533829);
void live_h_4(double *state, double *unused, double *out_873783986208232369);
void live_H_4(double *state, double *unused, double *out_6553766903658930946);
void live_h_9(double *state, double *unused, double *out_1294967175197148676);
void live_H_9(double *state, double *unused, double *out_733452031605516524);
void live_h_10(double *state, double *unused, double *out_205006460607604927);
void live_H_10(double *state, double *unused, double *out_1112432945254149966);
void live_h_12(double *state, double *unused, double *out_4053631831362776715);
void live_H_12(double *state, double *unused, double *out_5511718793007887674);
void live_h_31(double *state, double *unused, double *out_8269375197882396583);
void live_H_31(double *state, double *unused, double *out_3858924442348533255);
void live_h_32(double *state, double *unused, double *out_2605664927509823398);
void live_H_32(double *state, double *unused, double *out_8786515259153405803);
void live_h_13(double *state, double *unused, double *out_709067257945039309);
void live_H_13(double *state, double *unused, double *out_5411788331758843461);
void live_h_14(double *state, double *unused, double *out_1294967175197148676);
void live_H_14(double *state, double *unused, double *out_733452031605516524);
void live_h_33(double *state, double *unused, double *out_73942277000080270);
void live_H_33(double *state, double *unused, double *out_7009481446987390859);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}