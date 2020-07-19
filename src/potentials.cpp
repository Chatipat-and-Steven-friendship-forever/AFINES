#include "potentials.h"
#include "globals.h"

// modified from LAMMPS src/angle_harmonic.cpp
bend_result_type bend_harmonic(
        double kb, double theta0,
        array<double, 2> delr1,
        array<double, 2> delr2)
{
    bend_result_type result;

    // 1st bond
    double rsq1 = delr1[0] * delr1[0] + delr1[1] * delr1[1];
    double r1 = sqrt(rsq1);

    // 2nd bond
    double rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1];
    double r2 = sqrt(rsq2);

    // angle (cos and sin)
    double c = (delr1[0] * delr2[0] + delr1[1] * delr2[1]) / (r1 * r2);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    double s = sqrt(1.0 - c * c);
    if (s < maxSmallAngle) s = maxSmallAngle;

    double theta = acos(c) - theta0;

    // energy
    result.energy = 0.5 * kb * theta * theta;

    // force
    double a = -kb * theta / s;
    double a11 = a * c / rsq1;
    double a12 = -a / (r1 * r2);
    double a22 = a * c / rsq2;

    result.force1[0] = a11 * delr1[0] + a12 * delr2[0];
    result.force1[1] = a11 * delr1[1] + a12 * delr2[1];
    result.force2[0] = a22 * delr2[0] + a12 * delr1[0];
    result.force2[1] = a22 * delr2[1] + a12 * delr1[1];

    return result;
}


double bend_harmonic_energy(
        double kb, double theta0,
        array<double, 2> delr1,
        array<double, 2> delr2)
{
    // 1st bond
    double rsq1 = delr1[0] * delr1[0] + delr1[1] * delr1[1];
    double r1 = sqrt(rsq1);

    // 2nd bond
    double rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1];
    double r2 = sqrt(rsq2);

    // angle (cos and sin)
    double c = (delr1[0] * delr2[0] + delr1[1] * delr2[1]) / (r1 * r2);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    double theta = acos(c) - theta0;

    // energy
    double energy = 0.5 * kb * theta * theta;

    return energy;
}
