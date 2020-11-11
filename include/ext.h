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
        virtual ~external() {}
        virtual void update(double t) {}
        virtual ext_result_type compute(vec_type pos)
        {
            throw std::logic_error("not implemented");
        }
};

class ext_circle : public external
{
    public:
        ext_circle(double k, double r);
        ext_result_type compute(vec_type pos) override;
    protected:
        double k, r;
};

class ext_circle_linear : public ext_circle
{
    public:
        ext_circle_linear(double k0, double r0, double dk, double dr);
        void update(double t) override;
    protected:
        double k0, r0, dk, dr;
};

class ext_ring : public external
{
    public:
        ext_ring(double k, double r1, double r2);
        ext_result_type compute(vec_type pos) override;
    protected:
        double k, r1, r2;
};

class ext_circle_fene : public external
{
    public:
        ext_circle_fene(double k, double rmin, double rmax);
        ext_result_type compute(vec_type pos) override;
    protected:
        double k, rmin, rmax;
};

class ext_circle_fene_linear : public ext_circle_fene
{
    public:
        ext_circle_fene_linear(double k0, double rmin0, double rmax0, double dk, double drmin, double drmax);
        void update(double t) override;
    protected:
        double k0, rmin0, rmax0, dk, drmin, drmax;
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
