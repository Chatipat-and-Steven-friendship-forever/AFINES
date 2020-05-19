/*
 * motor.h
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __MOTOR_ENSEMBLE_H_INCLUDED__
#define __MOTOR_ENSEMBLE_H_INCLUDED__

//=====================================
// forward declared dependencies
//class filament_ensemble;

//=====================================
//included dependences
#include "motor.h"

//motor ensemble class

class motor_ensemble
{
    public:

        motor_ensemble(vector<vector<double>> motors, double delta_t, double temp,
                double mlen, filament_ensemble * network, double v0, double stiffness, double max_ext_ratio,
                double ron, double roff, double rend,
                double fstall, double rcut,
                double vis, bool use_attach_opt_);

        ~motor_ensemble();

        int get_nmotors();

        void check_broken_filaments();

        void motor_walk(double t);

        void motor_update();

        void update_d_strain(double g);
        
        void update_energies();
        
        double get_potential_energy();

        array<array<double, 2>, 2> get_virial();

        void motor_write(ostream& fout);

        void print_ensemble_thermo();
        
        void motor_tension(ofstream& fout);

        void add_motor(motor * m);

        void set_shear(double g);
        
        void kill_heads(int i);

        void unbind_all_heads();

        void revive_heads();

    private:

        double mld, gamma, tMove;
        double ke, pe, v;

        array<double, 2> fov;
        filament_ensemble *f_network;
        vector<motor *> n_motors;  
        bool use_attach_opt;
        array<array<double, 2>, 2> virial;
};

#endif
