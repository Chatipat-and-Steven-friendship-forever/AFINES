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
