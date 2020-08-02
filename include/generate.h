#ifndef AFINES_GENERATE_H
#define AFINES_GENERATE_H

#include "filament_ensemble.h"
#include "motor_ensemble.h"

vector<vector<double>> generate_filament_ensemble(
        box *bc, int npolymer, int nbeads_min, int nbeads_extra, double nbeads_extra_prob,
        double temp, double radius, double l0,
        vector<array<double, 3>> pos_sets, double kb);

vector<vector<double>> generate_filament_ensemble(
        box *bc, double density, double temp, double radius, int nbeads, double l0,
        vector<array<double, 3>> pos_sets, double kb);

vector<vector<double>> generate_motor_ensemble(box *bc, double density, double l0, vector<array<double, 3>> positions);

vector<vector<double>> spring_spring_intersections(box *bc, vector<vector<double>> beads, double len, double prob);

#endif
