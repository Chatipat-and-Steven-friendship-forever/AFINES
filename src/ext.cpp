#include "ext.h"

// [circle]

ext_circle::ext_circle(double k_, double r_)
{
    k = k_;
    r = r_;
}

ext_result_type ext_circle::compute(vec_type pos)
{
    ext_result_type result;
    double lsq = abs2(pos);
    if (lsq < r * r) return result;
    double l = sqrt(lsq);
    double dl = l - r;
    result.energy = 0.5 * k * dl * dl;
    result.force = (-k * dl / l) * pos;
    result.virial = -0.5 * outer(pos, result.force);
    return result;
}

// [circle linear]

ext_circle_linear::ext_circle_linear(double k0_, double r0_, double dk_, double dr_)
    : ext_circle(k0_, r0_)
{
    k0 = k0_;
    r0 = r0_;
    dk = dk_;
    dr = dr_;
}

void ext_circle_linear::update(double t)
{
    k = k0 + dk * t;
    r = r0 + dr * t;
}

// [ring]

ext_ring::ext_ring(double k_, double r1_, double r2_)
{
    k = k_;
    r1 = r1_;
    r2 = r2_;
}

ext_result_type ext_ring::compute(vec_type pos)
{
    ext_result_type result;
    double lsq = abs2(pos);
    if (lsq < r1 * r1) {
        double l = sqrt(lsq);
        double dl = l - r1;
        result.energy = 0.5 * k * dl * dl;
        result.force = (-k * dl / l) * pos;
        result.virial = -0.5 * outer(pos, result.force);
        return result;
    }
    if (lsq > r2 * r2) {
        double l = sqrt(lsq);
        double dl = l - r2;
        result.energy = 0.5 * k * dl * dl;
        result.force = (-k * dl / l) * pos;
        result.virial = -0.5 * outer(pos, result.force);
        return result;
    }
    return result;
}

// [circle fene]

ext_circle_fene::ext_circle_fene(double k_, double rmin_, double rmax_)
{
    k = k_;
    rmin = rmin_;
    rmax = rmax_;
}

ext_result_type ext_circle_fene::compute(vec_type pos)
{
    ext_result_type result;
    double lsq = abs2(pos);
    if (lsq < rmin * rmin) return result;
    double l = sqrt(lsq);
    double dl = l - rmin;
    double ratio = dl / rmax;
    result.force = -k * dl / (1.0 - ratio * ratio) * pos;
    result.energy = -0.5 * k * rmax * rmax * log1p(-ratio * ratio);
    result.virial = -0.5 * outer(pos, result.force);
    return result;
}

// [circle fene linear]

ext_circle_fene_linear::ext_circle_fene_linear(
        double k0_, double rmin0_, double rmax0_,
        double dk_, double drmin_, double drmax_)
    : ext_circle_fene(k0_, rmin0_, rmax0_)
{
    k0 = k0_;
    rmin0 = rmin0_;
    rmax0 = rmax0_;
    dk = dk_;
    drmin = drmin_;
    drmax = drmax_;
}

void ext_circle_fene_linear::update(double t)
{
    k = k0 + dk * t;
    rmin = rmin0 + drmin * t;
    rmax = rmax0 + drmax * t;
}

// [rectangle]

ext_rectangle::ext_rectangle(double k_, double x_, double y_)
{
    k = k_;
    x = x_;
    y = y_;
}

ext_result_type ext_rectangle::compute(vec_type pos)
{
    ext_result_type result;
    if (pos.x > 0.5 * x) {
        double dx = x - pos.x;
        result.energy += 0.5 * k * dx * dx;
        result.force.x -= k * dx;
    } else if (pos.x < -0.5 * x) {
        double dx = x + pos.x;
        result.energy += 0.5 * k * dx * dx;
        result.force.x -= k * dx;
    }
    if (pos.y > 0.5 * y) {
        double dy = y - pos.y;
        result.energy += 0.5 * k * dy * dy;
        result.force.y -= k * dy;
    } else if (pos.y < -0.5 * y) {
        double dy = y + pos.y;
        result.energy += 0.5 * k * dy * dy;
        result.force.y -= k * dy;
    }
    result.virial = -0.5 * outer(pos, result.force);
    return result;
}

// [xperiodic]

ext_xperiodic::ext_xperiodic(double k_, double y_)
{
    k = k_;
    y = y_;
}

ext_result_type ext_xperiodic::compute(vec_type pos)
{
    ext_result_type result;
    if (pos.y > 0.5 * y) {
        double dy = y - pos.y;
        result.energy += 0.5 * k * dy * dy;
        result.force.y -= k * dy;
    } else if (pos.y < -0.5 * y) {
        double dy = y + pos.y;
        result.energy += 0.5 * k * dy * dy;
        result.force.y -= k * dy;
    }
    result.virial.yy = -0.5 * pos.y * result.force.y;
    return result;
}
