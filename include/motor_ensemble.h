/*
 * motor.h
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#ifndef AFINES_MOTOR_ENSEMBLE_H
#define AFINES_MOTOR_ENSEMBLE_H

#include "motor.h"

class motor_ensemble
{
    public:

        motor_ensemble(vector<vector<double>> motors, double delta_t, double temp,
                double mlen, filament_ensemble * network, double v0, double stiffness, double max_ext_ratio,
                double ron, double roff, double rend, double fstall, double rcut, double vis);

        ~motor_ensemble();

        // flags
        void use_attach_opt(bool flag);
        void use_static(bool flag);
        void use_shear(bool flag);

        // settings
        void kill_heads(int i);
        void set_binding_two(double, double, double);
        void set_bending(double kb, double th0);

        // thermo
        double get_potential_energy();
        double get_kinetic_energy(); 
        virial_type get_virial();

        // output
        vector<vector<double>> output();
        void motor_write(ostream& fout);
        void motor_write_doubly_bound(ostream& fout);
        void print_ensemble_thermo();
        void motor_tension(ofstream& fout);

        // state
        int get_nmotors();
        motor *get_motor(int);

        // shear
        void update_d_strain(double g);

        // update
        void integrate();
        void montecarlo();
        void update_energies();
        void compute_forces();
        void unbind_all_heads();
        void revive_heads();

        // misc
        void add_motor(motor *m);

    protected:
        filament_ensemble *f_network;
        vector<motor *> n_motors;
        double mld;

        // flags
        bool attach_opt_flag, shear_flag, static_flag;

        // thermo
        double ke_vel, ke_vir, pe;
        virial_type virial;
};

#endif
