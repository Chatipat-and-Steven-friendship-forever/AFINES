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

class bead {
    public:
        bead();
        bead(double xcm, double ycm, double len, double vis); 
        bead(const bead& other);
        ~bead();

        // state

        array<double, 2> get_pos();
        void set_pos(array<double, 2> pos);
        double get_xcm();
        double get_ycm();
        void set_xcm(double xcm);
        void set_ycm(double ycm);

        double get_fx();
        double get_fy();
        void set_fx(double);
        void set_fy(double);

        array<double, 2> get_force();
        void update_force(array<double, 2> f);
        void reset_force();

        // parameters

        double get_length();
        double get_friction();
        double get_viscosity();

        // output

        vector<double> output();
        string write();
        string to_string();

        // misc

        bool operator==(const bead& that);

    private:
        // state
        double x, y;
        array<double, 2> force;

        // parameters
        double rad, visc, friction;
};

#endif
