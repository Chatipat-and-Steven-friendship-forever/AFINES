/*
 *  filament.cpp
 *  
 *
 *  Created by Simon Freedman on 12/22/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

#ifndef AFINES_FILAMENT_H
#define AFINES_FILAMENT_H

class filament_ensemble;

#include "spring.h"
#include "globals.h"

class filament
{
    public:
        filament(filament_ensemble *net, vector<vector<double>> beadvec, double spring_length,
                double stretching_stiffnes, double ext, double bending_stiffness,
                double deltat, double temp, double fracture);
        ~filament();

        // output

        vector<vector<double>> output_beads(int fil);
        vector<vector<double>> output_springs(int fil);
        vector<double> output_thermo(int fil);

        string write_beads(int fil);
        string write_springs(int fil);
        string write_thermo(int fil);

        void print_thermo();
        string to_string();

        // thermo

        double get_stretching_energy();
        double get_bending_energy();
        double get_potential_energy();
        double get_kinetic_energy_vir();
        double get_kinetic_energy_vel();
        double get_total_energy();

        virial_type get_stretching_virial();
        virial_type get_bending_virial();

        // state

        int get_nbeads();
        vec_type get_bead_position(int bead);

        int get_nsprings();
        spring *get_spring(int i);

        double get_end2end();
        vec_type get_force(int i);

        // dynamics

        void update_d_strain(double);
        void update_positions();
        vector<filament *> fracture(int node);
        vector<filament *> try_fracture();
        void detach_all_motors();

        // internal forces

        void update_stretching();
        void update_bending();

        // external forces

        void pull_on_ends(double f);
        void affine_pull(double f);
        void update_forces(int index, vec_type f);

        // misc

        box *get_box();

        void add_bead(vector<double> a, double l0, double kl, double me);

        void init_ubend(); // computes bending energies only
        inline double angle_between_springs(int i, int j);

        vector<vector<double>> get_beads(size_t first, size_t last); // used for fractures

        bool operator==(const filament& that);

        // growing

        void set_l0_max(double);
        void set_nsprings_max(int);
        void set_l0_min(double);
        void set_kgrow(double);
        void set_lgrow(double);

        void update_length();
        void grow(double);

    protected:
        box *bc;
        filament_ensemble *filament_network;

        // state
        vector<vec_type> prv_rnds;
        vector<class bead *> beads;
        vector<spring *> springs;

        // thermo
        double ke_vir, ke_vel;
        double ubend;
        virial_type bending_virial;

        // parameters
        double kb, temperature, dt, fracture_force, damp;

        // growing parameters
        int nsprings_max;
        double spring_l0, l0_max, l0_min, kgrow, lgrow;

        // precompute
        double bd_prefactor;
        double fracture_force_sq;
};

#endif
