#ifndef AFINES_POTENTIALS_H
#define AFINES_POTENTIALS_H

#include "globals.h"

struct bend_result_type {
    double energy;
    array<double, 2> force1;
    array<double, 2> force2;
};

bend_result_type bend_harmonic(
        double kb, double theta0,
        array<double, 2> delr1, array<double, 2> delr2);

double bend_harmonic_energy(
        double kb, double theta0,
        array<double, 2> delr1, array<double, 2> delr2);

#endif
