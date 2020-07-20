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
#ifndef AFINES_MOTOR_ENSEMBLE_H
#define AFINES_MOTOR_ENSEMBLE_H

//=====================================
// forward declared dependencies
//class filament_ensemble;

//=====================================
//included dependences
#include "motor.h"
#include "spacer.h"

template <class motor_type>
class motor_ensemble
{
    public:

        motor_ensemble(vector<vector<double>> motors, double delta_t, double temp,
                double mlen, filament_ensemble * network, double v0, double stiffness, double max_ext_ratio,
                double ron, double roff, double rend, double fstall, double rcut, double vis);

        ~motor_ensemble();

        void use_attach_opt(bool flag);

        void use_shear(bool flag);

        void use_static(bool flag);

        int get_nmotors();

        motor *get_motor(int);

        void check_broken_filaments();

        void motor_update();

        void update_d_strain(double g);

        void update_energies();

        double get_potential_energy();

        double get_kinetic_energy(); 

        virial_type get_virial();

        vector<vector<double>> output();

        void motor_write(ostream& fout);

        void motor_write_doubly_bound(ostream& fout);

        void print_ensemble_thermo();

        void motor_tension(ofstream& fout);

        void add_motor(motor_type *m);

        void set_shear(double g);

        void kill_heads(int i);

        void set_binding_two(double, double, double);

        void unbind_all_heads();

        void revive_heads();

    protected:
        double mld;
        double ke_vel, ke_vir, pe;
        filament_ensemble *f_network;
        vector<motor_type *> n_motors;
        bool attach_opt_flag, shear_flag, static_flag;
        virial_type virial;
};

class spacer_ensemble : public motor_ensemble<spacer>
{
    public:
        using motor_ensemble<spacer>::motor_ensemble;    
        void set_bending(double, double);

};

class xlink_ensemble : public motor_ensemble<motor>
{
    public:
        using motor_ensemble<motor>::motor_ensemble;    

};

#endif
