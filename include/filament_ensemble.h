/*
 *  filament_ensemble.h
 *  
 *
 *  Authors : Shiladitya Banerjee, Simon Freedman
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#ifndef AFINES_FILAMENT_ENSEMBLE_H
#define AFINES_FILAMENT_ENSEMBLE_H

#include "filament.h"
#include "box.h"
#include "quadrants.h"
#include "exv.h"

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

        // monte carlo
        void set_growing(double, double, double, double, int);
        void try_grow();
        void try_fracture();
        void montecarlo();

        // dynamics
        void integrate();
        void update_positions();
        void update_d_strain(double);

        // update forces/energies
        void compute_forces();
        void update_stretching();
        void update_bending();
        void update_external();
        void update_excluded_volume();
        void update_forces(int fil, int bead, vec_type f);
        vec_type external_force(vec_type pos);
        void update_energies();

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

    protected:
        box *bc;
        quadrants *quads;
        excluded_volume *exv;
        vector<filament *> network;

        int nsprings_per_fil_max;

        // parameters
        double dt, temperature, spring_rest_len, visc;

        // thermo
        double pe_stretch, pe_bend, pe_exv, ke_vel, ke_vir;
        virial_type vir_stretch, vir_bend;

        // external forces
        external_force_type external_force_flag;
        double circle_wall_radius, circle_wall_spring_constant;
};

#endif
