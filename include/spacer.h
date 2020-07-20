/*
 * spacer.h
 *  
 *
 *  Created by Simon Freedmane on 8/8/16.
 *  Copyright 2016 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef AFINES_SPACER_H
#define AFINES_SPACER_H

//=====================================
// forward declared dependencies

//=====================================
//included dependences
#include "motor.h"

//spacer class

class spacer : public motor {
    public:

        spacer(vector<double> svec, double mlen, filament_ensemble* network,
                double delta_t, double temp, double v0, double stiffness, double max_ext_ratio,
                double ron, double roff, double rend,
                double fstall, double rcut, double vis);

        virtual ~spacer();

        void set_bending(double, double);

        double get_kb();
        
        double get_th0();

        virtual void update_force();
        
        void identify();
        
        int get_further_end(int, int, int);

        vec_type disp_from_bead(int, int, int);

        void update_bending(int);

        array<vec_type, 2> get_b_force();

        virtual void brownian_relax(int hd);

        void filament_update();

        virtual double metropolis_prob(int hd, array<int, 2> flidx, vec_type newpos, double maxRate);

        virtual bool allowed_bind(int hd, array<int, 2> flidx);

        array<double,2> get_bending_energy();

    private:
        double kb, th0;
        array<vec_type, 2> b_force;
        array<double, 2> b_eng;
};

#endif
