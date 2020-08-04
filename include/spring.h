/*
 *  spring.h
 *
 *  Created by Simon Freedman on 9/10/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

#ifndef AFINES_SPRING_H
#define AFINES_SPRING_H

#include "globals.h"

class spring
{
    public:
        spring() {}
        spring(class box *bc, double l0, double kl);
        virtual ~spring() {}

        // updates all state
        void step(vec_type h0, vec_type h1);

        // [state]

        vec_type get_h0();
        vec_type get_h1();

        double get_length();
        vec_type get_disp();
        vec_type get_direction();

        // -dU/dh0
        vec_type get_force();

        // [parameters]

        double get_kl();

        double get_l0();
        void set_l0(double myl0);

        // [thermo]

        double get_stretching_energy();

        virial_type get_virial();

        // [misc]

        // compare parameters
        bool operator==(const spring& that);

        vec_type intpoint(vec_type pos);
        bool get_line_intersect(spring *l2); 

    protected:

        // state
        vec_type h0, h1;

        // derived state
        double llen;
        vec_type disp, direc;

        vec_type force;

        // parameters
        double kl, l0;
        class box *bc;
};

#endif
