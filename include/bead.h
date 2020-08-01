/*
 *  bead.cpp
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#ifndef AFINES_BEAD_H
#define AFINES_BEAD_H

#include "globals.h"
#include "vec.h"

class bead {
    public:
        bead(double xcm, double ycm, double len, double vis); 
        bead(const bead& other);

        // state

        vec_type get_pos();
        void set_pos(vec_type pos);

        vec_type get_force();
        void update_force(vec_type f);
        void reset_force();

        // parameters

        double get_length();
        double get_friction();
        double get_viscosity();

        // misc

        bool operator==(const bead& that);

    private:
        // state
        vec_type pos, force;

        // parameters
        double rad, visc, friction;
};

#endif
