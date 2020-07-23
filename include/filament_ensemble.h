/*
 *  filament_ensemble.h
 *  
 *
 *  Authors : Shiladitya Banerjee, Simon Freedman
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef AFINES_FILAMENT_ENSEMBLE_H
#define AFINES_FILAMENT_ENSEMBLE_H

//=====================================
// forward declared dependencies

//=====================================
//included dependences
#include "filament.h"
#include "box.h"
#include "quadrants.h"
//=====================================
//filament network class

enum class external_force_type
{
    none,
    circle
};

class filament_ensemble
{
    public:

        filament_ensemble(box *bc, vector< vector<double> > beads, array<int, 2> mynq, double delta_t, double temp,
                double vis, double spring_len, double stretching, double ext, double bending, double frac_force, double RMAX, double A);

        ~filament_ensemble();

        // settings
        void set_circle_wall(double radius, double spring_constant);
        void set_fene_dist_pct(double);

        // quadrants
        quadrants *get_quads();
        void quad_update_serial();
        vector<array<int, 2>> *get_attach_list(vec_type pos);

        // state

        int get_nbeads();
        int get_nsprings();
        int get_nfilaments();

        box *get_box();
        vector<filament *> * get_network();
        filament * get_filament(int index);

        vec_type get_force(int fil, int spring);
        vec_type get_direction(int fil, int spring);
        vec_type get_start(int fil, int spring);
        vec_type get_end(int fil, int spring);

        double get_int_direction(int fil, int spring, double xp, double yp);
        double get_llength(int fil, int spring);

        bool is_polymer_start(int f, int a);

        // thermo
        double get_stretching_energy();
        virial_type get_stretching_virial();
        double get_bending_energy();
        virial_type get_bending_virial();
        double get_exv_energy();
        double get_kinetic_energy_vel();
        double get_kinetic_energy_vir();

        // dynamics
        void update();
        void update_positions();
        void update_d_strain(double);

        // update forces/energies

        void update_stretching();
        void update_filament_stretching(int); // fragments also
        void update_bending();
        void update_int_forces();
        void update_forces(int fil, int bead, vec_type f);
        vec_type external_force(vec_type pos);
        void update_energies();

        // excluded volume
        void update_spring_forces(int f); 
        void update_spring_forces_from_quads(); 
        void update_force_between_filaments(double n1, double l1, double n2, double l2); 
        void update_excluded_volume(int f); 

        // output

        vector<vector<double>> output_beads();
        vector<vector<double>> output_springs();
        vector<vector<double>> output_thermo();

        void write_beads(ofstream& fout);
        void write_springs(ofstream& fout);
        void write_thermo(ofstream& fout);

        void print_filament_thermo();
        void print_network_thermo();
        void print_filament_lengths();

        // fracture
        vector<int> get_broken();
        void clear_broken();

        // growing
        void set_growing(double, double, double, double, int);

    protected:
        box *bc;
        quadrants *quads;
        vector<filament *> network;

        // parameters
        double dt, temperature, spring_rest_len, visc;

        // thermo
        double pe_stretch, pe_bend, pe_exv, ke_vel, ke_vir;
        virial_type vir_stretch, vir_bend;

        // excluded volume
        double rmax, kexv;
        int nsprings_per_fil_max;

        // external forces
        external_force_type external_force_flag;
        double circle_wall_radius, circle_wall_spring_constant;

        // fracture
        vector<int> broken_filaments;
};

#endif
