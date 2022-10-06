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
void live_H(double *in_vec, double *out_7569511977236190551);
void live_err_fun(double *nom_x, double *delta_x, double *out_6649131745326665324);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_1500566119272722476);
void live_H_mod_fun(double *state, double *out_3395029699315419861);
void live_f_fun(double *state, double dt, double *out_5012703960485696517);
void live_F_fun(double *state, double dt, double *out_3176767779633191742);
void live_h_4(double *state, double *unused, double *out_236166562164894611);
void live_H_4(double *state, double *unused, double *out_7675434606189476881);
void live_h_9(double *state, double *unused, double *out_8283415105117256638);
void live_H_9(double *state, double *unused, double *out_7434244959559886236);
void live_h_10(double *state, double *unused, double *out_7543635371771011316);
void live_H_10(double *state, double *unused, double *out_8600325372552353718);
void live_h_12(double *state, double *unused, double *out_4381482374354067233);
void live_H_12(double *state, double *unused, double *out_2655978198157515086);
void live_h_31(double *state, double *unused, double *out_462198029778151580);
void live_H_31(double *state, double *unused, double *out_89584834167498623);
void live_h_32(double *state, double *unused, double *out_4585011177888638454);
void live_H_32(double *state, double *unused, double *out_2796896588680955873);
void live_h_13(double *state, double *unused, double *out_8483101578250308525);
void live_H_13(double *state, double *unused, double *out_823566997576918527);
void live_h_14(double *state, double *unused, double *out_8283415105117256638);
void live_H_14(double *state, double *unused, double *out_7434244959559886236);
void live_h_33(double *state, double *unused, double *out_7253605793644970757);
void live_H_33(double *state, double *unused, double *out_3240141838806356227);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}