/*
 *  spring.h
 *
 *  Created by Simon Freedman on 9/10/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

#ifndef AFINES_SPRING_H
#define AFINES_SPRING_H

class filament;

#include "globals.h"
#include "box.h"
#include "motor.h"

class spring {
    public:
        spring() {}
        spring(double l0, double kl, filament* f, array<int, 2> aindex);
        virtual ~spring() {}

        // [dynamics]

        // reads positions from filament
        // updates all state except force
        void step();

        // updates force
        void update_force();

        // applies force to filament
        void filament_update();

        // [state]

        vec_type get_h0();
        vec_type get_h1();

        double get_length();
        vec_type get_disp();
        vec_type get_direction();

        vec_type get_force();

        // [parameters]

        double get_kl();

        double get_l0();
        void set_l0(double myl0);

        array<int, 2> get_aindex();
        void set_aindex(array<int, 2> idx);
        void inc_aindex();

        // [thermo]

        double get_stretching_energy();

        virial_type get_virial();

        // [misc]

        bool operator==(const spring& that);
        bool is_similar(const spring& that);

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
        array<int, 2> aindex;  // bead indices
        double kl, l0;
        box *bc;
        filament *fil;
};

#endif
