/*
 *  bead.cpp
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef AFINES_BEAD_H
#define AFINES_BEAD_H

//=====================================
// forward declared dependencies

//=====================================
//included dependences
#include "globals.h"

//=====================================
//filament bead class
class bead
{
    public:

        bead();
        
        bead(double xcm, double ycm, double len, double vis); 
        
        bead(const bead& other);
        
        ~bead();
    
        void update();

        double get_length();

        double get_xcm();
        
        double get_ycm();
        
        void update_force(double f1, double f2);
        
        void reset_force();

        double get_fx();

        double get_fy();

        void set_fx(double);

        void set_fy(double);

        array<double,2> get_force();

        double get_friction();
        
        double get_viscosity();

        void set_xcm(double xcm);

        void set_ycm(double ycm);

        vector<double> output();

        string write();
        
        string to_string();
        
        bool operator==(const bead& that);    
        

    private:
        
        double x, y, rad, visc, friction;

        array<double, 2> force;
        
};

#endif
