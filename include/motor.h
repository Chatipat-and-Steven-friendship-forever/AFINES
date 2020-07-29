/*
 * motor.h
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#ifndef AFINES_MOTOR_H
#define AFINES_MOTOR_H

class filament_ensemble;

#include "globals.h"
#include "box.h"
#include "ext.h"

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

        // [settings]

        void set_binding_two(double ron2, double roff2, double rend2);

        void set_bending(double kb, double th0);
        double get_kb();
        double get_th0();

        void set_par(double k);
        void set_antipar(double k);
        void set_align(double k);

        void set_external(external *ext);

        void relax_head( int hd);
        void kill_head( int hd);
        void revive_head(int hd);
        void deactivate_head(int hd);

        // get [state]
        array<motor_state, 2> get_states();
        vec_type get_h0();
        vec_type get_h1();
        array<int, 2> get_f_index();
        array<int, 2> get_l_index();

        // calculations done by update_force
        vec_type get_force();  // spring forces
        array<vec_type, 2> get_b_force();  // bending forces

        // update [forces]
        void update_force();  // updates all forces

        // helper functions, used by update_force
        void update_force_fraenkel_fene();
        void update_bending(int hd);
        void update_alignment();
        void update_external(int hd);
        void update_force_proj(int hd);
        void filament_update();

        // [dynamics]
        void update_d_strain(double g);  // stress
        void brownian_relax(int hd);
        void walk(int hd);
        void step();  // update derived state

        // [attach/detach]
        double metropolis_prob(int hd, array<int, 2> fl_idx, vec_type newpos);
        double alignment_penalty(vec_type a, vec_type b);
        // attachment
        bool try_attach(int head, mc_prob &p);
        bool allowed_bind( int hd, array<int, 2> fl_idx);
        void attach_head(int hd, vec_type intpoint, array<int, 2> fl);
        // detachment
        bool try_detach(int head, mc_prob &p);
        vec_type generate_off_pos(int hd);
        void detach_head_without_moving(int hd);
        void detach_head(int hd, vec_type pos);

        // thermo
        double get_stretching_energy();
        array<double, 2> get_bending_energy();
        virial_type get_virial();
        double get_stretching_energy_fene();
        double get_kinetic_energy_vel();
        double get_kinetic_energy_vir();

        // output
        string to_string();
        vector<double> output();
        string write();

    protected:

        box *bc;
        filament_ensemble *filament_network;

        // [parameters]

        double dt, temperature, damp;
        double bd_prefactor;

        // attach/detach
        double kon, koff, kend;
        double kon2, koff2, kend2;
        double max_bind_dist, max_bind_dist_sq;

        // walk
        double vs, stall_force;

        // stretch
        double mk, mld;
        double max_ext, eps_ext; // fene

        // bend
        double kb, th0;

        // align
        double kalign;
        bool par_flag;

        // external
        external *ext;

        // [state]

        // head state
        array<motor_state, 2> state;
        array<vec_type, 2> h;  // for bound, updated by step
        // unbound only
        array<vec_type, 2> prv_rnd;
        // bound only
        array<fp_index_type, 2> fp_index;  // location bound
        array<vec_type, 2> ldir_bind, bind_disp;  // for unbinding

        // [derived] from state
        double len;
        vec_type disp, direc;

        // [forces] computed from state and derived
        double tension;
        vec_type force;
        array<vec_type, 2> b_force;
        array<vec_type, 2> ext_force;
        array<double, 2> f_proj;

        // [thermo] computed along with forces
        double ke_vel, ke_vir;
        array<double, 2> b_eng, ext_eng;
        double align_eng;
};

#endif
