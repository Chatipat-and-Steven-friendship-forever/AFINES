#include "shear.h"
#include "box.h"

// [sequence of shear protocols]

shear_sequence::shear_sequence(vector<shear_protocol *> p, vector<double> t)
{
    if (p.size() != t.size() + 1) {
        throw std::logic_error("times at which to switch the shear protocol should be one less than the number of shear protocols");
    }
    protocols = p;
    times = t;
}

shear_sequence::~shear_sequence()
{
    for (shear_protocol* s : protocols) {
        delete s;
    }
}

bool shear_sequence::need_virial(double t)
{
    size_t i = 0;
    for (; i < times.size(); i++) {
        if (i < times[i]) break;
    }
    return protocols[i]->need_virial(t);
}

double shear_sequence::compute(double t, virial_type virial)
{
    size_t i = 0;
    for (; i < times.size(); i++) {
        if (i < times[i]) break;
    }
    return protocols[i]->compute(t, virial);
}

// [constant strain]

strain_constant::strain_constant(double offset_)
{
    offset = offset_;
}

double strain_constant::compute(double t, virial_type dummy)
{
    return offset;
}

// [sine strain]

strain_sine::strain_sine(double amp_, double freq_, double t0_, double offset_)
{
    amp = amp_;
    freq = freq_;
    t0 = t0_;
    offset = offset_;
}

double strain_sine::compute(double t, virial_type dummy)
{
    return amp * sin(2.0 * pi * freq * (t - t0)) + offset;
}

// [triangle strain]

strain_triangle::strain_triangle(double amp_, double freq_, double t0_, double offset_)
{
    amp = amp_;
    freq = freq_;
    t0 = t0_;
    offset = offset_;
}

double strain_triangle::compute(double t, virial_type dummy)
{
    double phase = std::remainder(freq * (t - t0) + 0.25, 1.0);
    return 4.0 * amp * (std::fabs(phase) - 0.25) + offset;
}

// [linear strain]

strain_linear::strain_linear(double amp_, double freq_, double t0_, double offset_)
{
    amp = amp_;
    freq = freq_;
    t0 = t0_;
    offset = offset_;
}

double strain_linear::compute(double t, virial_type dummy)
{
    return amp * freq * (t - t0) + offset;
}

// [constant stress]

stress_constant::stress_constant(box *bc_, double dt_, double rate_, double value_)
{
    bc = bc_;
    dt = dt_;
    rate = rate_;  // rate at which to relax to value
    value = value_;  // target value of stress
}

bool stress_constant::need_virial(double t)
{
    return true;
}

double stress_constant::compute(double t, virial_type virial)
{
    // if strain = delrx/ybox, then
    // strain += stress_rate (-dU/dstrain / area + stress) dt
    // units:
    // - stress: energy / area
    // - stress_rate: area / (energy time)
    double xbox = bc->get_xbox();
    double ybox = bc->get_ybox();
    double delrx = bc->get_delrx();
    double f = -2.0 * virial.yx / ybox;
    return delrx + rate * (value + f / xbox) * ybox * dt;
}
