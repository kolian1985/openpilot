Version 8-SpecialAlza (2022-11-01)
=========================
 * Supported Perodua Alza

Version 5-longOverhaul  (2022-08-15)
=========================
 * bukapilot speed display matches stock speedometer
 * Cleaner settings UI
 * Rework of GPS time sync
 * New package feature thanks to @benmasato
 * Toyota: Enable stock Lane Departure Prevention when not engaged
 * Perodua Ativa & MFL: Adjustable lead car follow profile
 * From upstream openpilot v0.8.13:
   *  New driving model from upstream which is now trained on 1 million minutes of driving data
   *  Fixed lead training making lead predictions significantly more accurate
   *  Combined longitudinal planning now happens in a single MPC system
   *  New vision based forward collision warning
   *  New alert sounds
   *  New MPC acceleration lag compensation
   *  Fixed vehicleModelInvalid triggering due to false positives


Version 4-firstbatch  (2022-07-04)
=========================
 * Update UI: Training guide page, font change, add offroad alerts
 * Allow vision to leave power saver mode through ignition button and manual power on
 * Solve intermittent device power off
 * Change speed display to show GPS speed, NOT odometry speed
 * Toyota Alphard: Fix revving and wheelspeed scaling
 * Perodua PSD: Smoother low pump noise braking (lesser jerks)
 * Perodua PSD: Solve pump pressure bleeding issue
 * Perodua PSD: Enable stock Lane Departure Prevention when not engaged
 * Perodua Semi ACC: Add lateral tunes presets through settings
 * Support Toyota Corolla Cross 2021

Version 3-firstbatch  (2022-06-06)
=========================
 * Add full powersaver mode to 13 hours, battery can last at least 2 weeks idle with KommuAssist in power saver
 * Add stock ACC option in settings (only for Perodua)
 * Perodua PSD: Improve longitudinal
 * Perodua PSD: Add stock HUD warnings for front departure, forward collision warning & braking
 * Perodua PSD: Revert cruise speed set logic back to stock behaviour
 * Perodua PSD: Add brake pressure bleed warning
 * Perodua MG3: Fix odometer scaling
 * Increase standstill braking distance for all vehicles
 * QC: Add pre-fulfillment QC test
 * Fix minor UI display problems and release note popup
 * Add file corruption recovery mechanism during scons build

Version 2-firstbatch  (2022-05-11)
=========================
 * First bukapilot prebuilt release
 * Perodua: Remove horseriding
 * Perodua: 1.5x better longitudinal
 * Bug: Remove intermittent controls mismatch

Version 1-firstbatch  (2022-04-26)
=========================
 * Initial release of bukapilot based on openpilot v0.8.6
 * Supported Perodua Ativa 2020 - 2022 (AV)
 * Supported Perodua Myvi FL 2022 (AV)
 * Supported Toyota Alphard/Hybrid 2019-20
 * Supported Toyota Camry/Hybrid 2021-22
 * Supported Toyota Corolla Altis/Hybrid 2020-P
 * Supported Toyota Prius 2021-22
 * Supported Lexus ES/ES Hybrid 2019-21
 * Supported Lexus NX 2020
 * Supported Lexus UX/RX/Hybrid 2020-21
 * Supported Axia with KommuActuator2019 - 2022 (GXtra, Style, SE, AV)
 * Supported Bezza with KommuActuator 2020 - 2022 (X, AV)
 * Supported Myvi Gen 3 Pre-FL with KommuActuator 2017 - 2021 (G, X, H, AV)
