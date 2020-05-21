/*
 *  spring.h
 *
 *  Created by Simon Freedman on 9/10/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __LINK_H_INCLUDED__
#define __LINK_H_INCLUDED__


//=====================================
// forward declared dependencies
class filament;

//=====================================
//included dependences
#include "globals.h"
#include "box.h"

//=====================================
//spring class
class spring
{
    public:
        spring();
        
        spring(double len, double stiffness, double max_ext, filament* f, array<int, 2> aindex);
        
        virtual ~spring();

        array<double, 2> get_hx();
        
        array<double, 2> get_hy();
        
        double get_kl();
        
        double get_length();
        
        double get_length_sq();
       
        double get_l0();
        
        double get_fene_ext();
        
        double get_stretching_energy();

        array<array<double, 2>, 2> get_virial();
        
        double get_xcm();
        
        double get_ycm();
        
        string to_string();

        vector<double> output();

        string write();
        
        void step();
        
        void filament_update();
        
        bool operator==(const spring& that);    
        
        bool is_similar(const spring& that);    

        void update_force();
        
        void update_force_fraenkel_fene();
        
        double get_stretching_energy_fene();
        
        void update_force_marko_siggia(double kToverA);

        array<double,2> get_force();

        void set_aindex1(int i);
        
        double get_distance_sq(double xp, double yp);

        double get_int_angle(double xp, double yp);
        
        array<double,2> get_intpoint();

        void calc_intpoint(double xp, double yp);

        array<double, 2> get_direction();
        
        array<double, 2> get_disp();
        
        array<double, 2> get_neg_disp();

    protected:

        box *bc;

        double xcm, ycm, l0, kl, max_ext, eps_ext, llen, llensq;//, force;
       
        array<double,2> hx, hy;
        array<double, 2> disp, force, intpoint, direc;

        array<int, 2> aindex;
         
        filament *fil;
        
        vector< array<int,2> > quad; //vector of two vectors(x and y quadrants) of integers
};
#endif
