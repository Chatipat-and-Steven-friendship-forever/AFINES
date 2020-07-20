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
#ifndef AFINES_MOTOR_H
#define AFINES_MOTOR_H

//=====================================
// forward declared dependencies
class filament_ensemble;

//=====================================
//included dependences
#include "globals.h"
#include "box.h"

enum class motor_state {
    free = 0,
    bound = 1,
    dead = -1,
    inactive = -2
};

class motor
{
    public:

        motor() {}

        motor(vector<double> mvec,
                double mlen, filament_ensemble *network,
                double delta_t, double v0, double temp, double stiffness, double max_ext_ratio,
                double ron, double roff, double rend,
                double fstall, double rcut,
                double vis);

        virtual ~motor() {}

        virtual bool allowed_bind( int hd, array<int, 2> fl_idx);

        void attach_head(int hd, vec_type intpoint, array<int, 2> fl);

        bool attach( int hd);

        bool attach_opt(int hd);

        void relax_head( int hd);

        void kill_head( int hd);

        void identify();

        void update_position_attached(int hd);

        virtual void update_force();

        virtual void update_force_fraenkel_fene();

        void update_angle();

        virtual void brownian_relax(int hd);

        void step_onehead( int hd);

        void filament_update_hd(int hd, vec_type f);

        void filament_update();

        void update_d_strain(double g);

        array<int,2> get_f_index();

        array<int,2> get_l_index();

        void init_l_index(int hd, int idx);

        void set_l_index(int hd, int idx);

        void set_f_index(int hd, int idx);

        array<double,2> get_pos_a_end();

        void set_pos_a_end(int hd, double pos);

        double get_pos_a_end(int hd);

        vec_type get_force();

        array<motor_state, 2> get_states();

        vec_type get_h0();

        vec_type get_h1();

        void detach_head(int hd);

        void detach_head(int hd, vec_type pos);

        void detach_head_without_moving(int hd);

        void revive_head(int hd);

        void deactivate_head(int hd);

        void update_pos_a_end(int hd, double pos);

        double get_stretching_energy();

        virial_type get_virial();

        double get_stretching_energy_fene();

        double get_kinetic_energy_vel();

        double get_kinetic_energy_vir();

        virtual double metropolis_prob(int hd, array<int, 2> fl_idx, vec_type newpos, double maxprob);

        vec_type generate_off_pos(int hd);

        string to_string();

        vector<double> output();

        string write();

        void remove_from_spring(int hd);

        void add_to_spring(int hd);

        void inc_l_index(int hd);

        void set_binding_two(double ron2, double roff2, double rend2);

    protected:

        box *bc;
        filament_ensemble* filament_network;

        double mld, vs, stall_force, max_bind_dist, max_bind_dist_sq, mk, kon, koff, kend, dt, temperature, 
               damp, max_ext, eps_ext, ke_vel, ke_vir, bd_prefactor, tension, len, kon2, koff2, kend2;

        vec_type force, disp, direc;

        // head state
        array<motor_state, 2> state;
        array<int,2> f_index, l_index, spring_mot_idx;
        array<double, 2> pos_a_end;
        array<bool, 2> at_barbed_end;
        array<vec_type, 2> h, prv_rnd;
        array<vec_type, 2> ldir_bind, bind_disp;
};

#endif
