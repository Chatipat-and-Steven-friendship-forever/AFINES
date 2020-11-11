#ifndef AFINES_SHEAR_H
#define AFINES_SHEAR_H

#include "globals.h"

// interface

class shear_protocol
{
    public:
        virtual ~shear_protocol() {}
        virtual bool need_virial(double t) { return false; }
        virtual double compute(double t, virial_type virial) = 0;
};

// implementations

class shear_sequence : public shear_protocol
{
    public:
        shear_sequence(vector<shear_protocol *> protocols, vector<double> times);
        ~shear_sequence() override;
        bool need_virial(double t) override;
        double compute(double t, virial_type virial) override;
    protected:
        vector<shear_protocol *> protocols;
        vector<double> times;
};

class strain_constant : public shear_protocol
{
    public:
        strain_constant(double offset);
        double compute(double t, virial_type virial) override;
    protected:
        double offset;
};

class strain_sine : public shear_protocol
{
    public:
        strain_sine(double amp, double freq, double t0, double offset);
        double compute(double t, virial_type virial) override;
    protected:
        double amp, freq, t0, offset;
};

class strain_triangle : public shear_protocol
{
    public:
        strain_triangle(double amp, double freq, double t0, double offset);
        double compute(double t, virial_type virial) override;
    protected:
        double amp, freq, t0, offset;
};

class strain_linear : public shear_protocol
{
    public:
        strain_linear(double amp, double freq, double t0, double offset);
        double compute(double t, virial_type virial) override;
    protected:
        double amp, freq, t0, offset;
};

class stress_constant : public shear_protocol
{
    public:
        stress_constant(class box *bc, double dt, double rate, double value);
        bool need_virial(double t) override;
        double compute(double t, virial_type virial) override;
    protected:
        class box *bc;
        double dt, rate, value;
};

#endif
