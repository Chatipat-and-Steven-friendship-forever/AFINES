#ifndef AFINES_EXT_H
#define AFINES_EXT_H

#include "globals.h"

struct ext_result_type
{
    ext_result_type() { energy = 0.0; }

    double energy;
    vec_type force;
    virial_type virial;
};

class external
{
    public:
        external() {}
        virtual ext_result_type compute(vec_type pos) = 0;
};

class ext_circle : public external
{
    public:
        ext_circle(double k, double r);
        ext_result_type compute(vec_type pos) override;
    protected:
        double k, r;
};

class ext_rectangle : public external
{
    public:
        ext_rectangle(double k, double x, double y);
        ext_result_type compute(vec_type pos) override;
    protected:
        double k, x, y;
};

class ext_xperiodic : public external
{
    public:
        ext_xperiodic(double k, double y);
        ext_result_type compute(vec_type pos) override;
    protected:
        double k, y;
};

#endif
