#ifndef __GENERATE_H_INCLUDED__
#define __GENERATE_H_INCLUDED__

#include "filament_ensemble.h"
#include "motor_ensemble.h"

vector<vector<double>> generate_filament_ensemble(
        box *bc, int npolymer, int nbeads_min, int nbeads_extra, int nbeads_extra_prob,
        double dt, double temp, double radius, double l0,
        vector<array<double, 3>> pos_sets, double kb, double seed);

vector<vector<double>> generate_filament_ensemble(
        box *bc, double density, double dt, double temp, double radius, int nbeads, double l0,
        vector<array<double, 3>> pos_sets, double kb, double seed);

vector<vector<double>> generate_motor_ensemble(box *bc, double density, double l0, vector<array<double, 3>> positions);

#endif
