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
void live_H(double *in_vec, double *out_2878163332676000540);
void live_err_fun(double *nom_x, double *delta_x, double *out_1459947820727129062);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_4470484474002255370);
void live_H_mod_fun(double *state, double *out_7668485600411002239);
void live_f_fun(double *state, double dt, double *out_1547214029111368966);
void live_F_fun(double *state, double dt, double *out_288403762753222131);
void live_h_4(double *state, double *unused, double *out_7908376479176728560);
void live_H_4(double *state, double *unused, double *out_799079216110286149);
void live_h_9(double *state, double *unused, double *out_6558772355698767102);
void live_H_9(double *state, double *unused, double *out_557889569480695504);
void live_h_10(double *state, double *unused, double *out_1840225830705905829);
void live_H_10(double *state, double *unused, double *out_6875568155516334644);
void live_h_12(double *state, double *unused, double *out_3938035012682619485);
void live_H_12(double *state, double *unused, double *out_177980191062692482);
void live_h_31(double *state, double *unused, double *out_6069685274528503166);
void live_H_31(double *state, double *unused, double *out_2567582841262321227);
void live_h_32(double *state, double *unused, double *out_7831050498024517196);
void live_H_32(double *state, double *unused, double *out_5351241005680738820);
void live_h_13(double *state, double *unused, double *out_3491309275453140949);
void live_H_13(double *state, double *unused, double *out_1664464166163638855);
void live_h_14(double *state, double *unused, double *out_6558772355698767102);
void live_H_14(double *state, double *unused, double *out_557889569480695504);
void live_h_33(double *state, double *unused, double *out_4130482847495987501);
void live_H_33(double *state, double *unused, double *out_5718139845901178831);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}