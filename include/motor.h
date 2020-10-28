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
                double l0, filament_ensemble *network,
                double delta_t, double v0, double temp, double kl,
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
        void set_align(double k);
        void set_antipar(double k);
        double get_kalign();
        int get_align();

        void set_external(external *ext);

        void relax_head( int hd);
        void kill_head( int hd);
        void revive_head(int hd);
        void deactivate_head(int hd);

        void set_velocity(double v1, double v2);
        void set_stall_force(double f1, double f2);

        void set_occ(double occ);

        // get [state]
        array<motor_state, 2> get_states();
        vec_type get_h0();
        vec_type get_h1();
        array<int, 2> get_f_index();
        array<int, 2> get_l_index();

        // calculations done by update_force
        array<vec_type, 2> get_force();  // total forces
        array<vec_type, 2> get_s_force();  // spring forces
        array<vec_type, 2> get_b_force();  // bending forces
        array<vec_type, 2> get_ext_force();  // external forces
        array<double, 2> get_force_proj();  // projected forces for walking

        // update [forces]
        void clear_forces();
        void update_force();  // updates all forces

        // helper functions, used by update_force
        void update_bending(int hd);  // compute and partially apply bending forces
        void update_alignment();  // compute and apply alignment forces
        void update_external(int hd);  // compute external forces
        void update_force_proj(int hd);  // compute projected forces for walking
        void filament_update();  // apply remaining forces to filaments

        // [dynamics]
        void update_d_strain(double g);  // apply shear (called automatically by box)
        void brownian_relax(int hd);  // Brownian dynamics for free heads
        void walk(int hd);  // walking for bound heads
        void step();  // compute derived state (incl. bound head positions)

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

        // [thermo]

        double get_stretching_energy();
        array<double, 2> get_bending_energy();
        double get_alignment_energy();
        array<double, 2> get_external_energy();

        virial_type get_stretching_virial();
        virial_type get_bending_virial();
        virial_type get_alignment_virial();
        virial_type get_external_virial();

        // output
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
        double occ;

        // walk
        array<double, 2> vs, stall_force;

        // stretch
        double mk, mld;

        // bend
        double kb, th0;

        // align
        double kalign;
        int par_flag;

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
        array<vec_type, 2> force;
        array<vec_type, 2> s_force;
        array<vec_type, 2> b_force;
        array<vec_type, 2> ext_force;
        array<double, 2> f_proj;

        // [thermo] computed along with forces
        double s_eng;
        array<double, 2> b_eng;
        double align_eng;
        array<double, 2> ext_eng;
        virial_type vir_stretch, vir_bend, vir_align, vir_ext;
};

#endif
