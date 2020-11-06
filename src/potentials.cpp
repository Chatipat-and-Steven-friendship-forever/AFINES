#include "potentials.h"

// most of these methods are modified from LAMMPS

// from hsieh, jain, larson, jcp 2006; eqn 5
// adapted by placing a cutoff, analogous to LAMMPS src/bond_fene.cpp
stretch_fene::stretch_fene(double kl_, double l0_, double max_ext_ratio)
{
    kl = kl_;
    l0 = l0_;
    max_ext = max_ext_ratio * l0;
    eps_ext = 0.01 * max_ext;
}

// multiply by direction to get the force
double stretch_fene::tension(double len)
{
    double ext = abs(len - l0);
    double scaled_ext;
    if (max_ext - ext > eps_ext)
        scaled_ext = ext / max_ext;
    else
        scaled_ext = (max_ext - eps_ext) / max_ext;
    return kl / (1.0 - scaled_ext * scaled_ext) * (len - l0);
}

double stretch_fene::energy(double len)
{
    double ext = abs(len - l0);
    if (max_ext - ext > eps_ext)
        return -0.5 * kl * max_ext * max_ext * log(1.0 - (ext / max_ext) * (ext / max_ext));
    else
        return 0.25 * kl * ext * ext * (max_ext / eps_ext);
}


bend_result_type bend_angle(
        vec_type delr1,
        vec_type delr2)
{
    bend_result_type result;

    // 1st bond
    double rsq1 = abs2(delr1);
    double r1 = sqrt(rsq1);

    // 2nd bond
    double rsq2 = abs2(delr2);
    double r2 = sqrt(rsq2);

    // angle (cos and sin)
    double c = dot(delr1, delr2) / (r1 * r2);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    // energy
    result.energy = c;

    // force
    double a11 = c / rsq1;
    double a12 = -1 / (r1 * r2);
    double a22 = c / rsq2;

    result.force1 = a11 * delr1 + a12 * delr2;
    result.force2 = a22 * delr2 + a12 * delr1;

    return result;
}

double angle(
        vec_type delr1,
        vec_type delr2)
{
    // 1st bond
    double rsq1 = abs2(delr1);
    double r1 = sqrt(rsq1);

    // 2nd bond
    double rsq2 = abs2(delr2);
    double r2 = sqrt(rsq2);

    // angle (cos and sin)
    double c = dot(delr1, delr2) / (r1 * r2);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    return c;
}

bend_result_type bend_harmonic(
        double kb, double theta0,
        vec_type delr1,
        vec_type delr2)
{
    bend_result_type result;

    bend_result_type tmp = bend_angle(delr1, delr2);

    double c = tmp.energy;

    double s = sqrt(1.0 - c * c);
    if (s < maxSmallAngle) s = maxSmallAngle;

    double theta = acos(c) - theta0;

    // energy
    result.energy = 0.5 * kb * theta * theta;

    // force
    double a = -kb * theta / s;

    result.force1 = a * tmp.force1;
    result.force2 = a * tmp.force2;

    return result;
}

double bend_harmonic_energy(
        double kb, double theta0,
        vec_type delr1,
        vec_type delr2)
{
    double c = angle(delr1, delr2);

    double theta = acos(c) - theta0;

    // energy
    double energy = 0.5 * kb * theta * theta;

    return energy;
}
