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

        void set_binding_two(double ron2, double roff2, double rend2);

        void set_bending(double kb, double th0);
        double get_kb();
        double get_th0();

        // get state
        array<motor_state, 2> get_states();
        vec_type get_h0();
        vec_type get_h1();
        array<int, 2> get_f_index();
        array<int, 2> get_l_index();
        vec_type get_force();
        array<vec_type, 2> get_b_force();

        // attach/detach
        bool allowed_bind( int hd, array<int, 2> fl_idx);
        bool try_attach(int head, bool opt, mc_prob &p);
        bool attach(int hd, mc_prob &p);
        bool attach_opt(int hd, mc_prob &p);
        void attach_head(int hd, vec_type intpoint, array<int, 2> fl);
        bool try_detach(int head, mc_prob &p);
        void detach_head(int hd);
        void detach_head(int hd, vec_type pos);
        void detach_head_without_moving(int hd);
        virtual double metropolis_prob(int hd, array<int, 2> fl_idx, vec_type newpos, double maxprob);
        vec_type generate_off_pos(int hd);

        // step/walk
        void update_angle();
        void walk(int hd);
        void update_pos_a_end(int hd, double pos);
        void brownian_relax(int hd);
        void set_pos_a_end(int hd, double pos);

        // from filament
        void update_position_attached(int hd);

        // to filament
        void filament_update_hd(int hd, vec_type f);
        void filament_update();

        // potential
        void update_force();
        void update_force_fraenkel_fene();
        void update_bending(int hd);

        // stress
        void update_d_strain(double g);

        // bound spring
        void init_l_index(int hd, int idx);
        void inc_l_index(int hd);
        void add_to_spring(int hd);
        void remove_from_spring(int hd);
        void set_f_index(int hd, int idx);
        void set_l_index(int hd, int idx);

        // thermo
        double get_stretching_energy();
        virial_type get_virial();
        double get_stretching_energy_fene();
        array<double, 2> get_bending_energy();
        double get_kinetic_energy_vel();
        double get_kinetic_energy_vir();

        // output
        string to_string();
        vector<double> output();
        string write();

        // misc
        void relax_head( int hd);
        void kill_head( int hd);
        void revive_head(int hd);
        void deactivate_head(int hd);

    protected:

        box *bc;
        filament_ensemble* filament_network;

        // parameters
        double dt, temperature, damp;

        // attach/detach parameters
        double kon, koff, kend;
        double kon2, koff2, kend2;

        // thermo
        double ke_vel, ke_vir;
        array<double, 2> b_eng;

        // potential parameters
        double mk, mld;
        double max_ext, eps_ext;
        double max_bind_dist;
        double kb, th0;

        // walking parameters
        double vs, stall_force;

        // precompute
        double bd_prefactor;
        double max_bind_dist_sq;

        // derived state
        double len;
        vec_type disp, direc;

        // forces
        double tension;
        vec_type force;
        array<vec_type, 2> b_force;

        // head state
        array<motor_state, 2> state;
        array<vec_type, 2> h;
        // unbound only
        array<vec_type, 2> prv_rnd;
        // bound only
        array<fp_index_type, 2> fp_index;
        array<vec_type, 2> ldir_bind, bind_disp;
};

#endif
