/*
 *  filament.cpp
 *  
 *
 *  Created by Simon Freedman on 12/22/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef AFINES_FILAMENT_H
#define AFINES_FILAMENT_H

//=====================================
// forward declared dependencies
class filament_ensemble;

//=====================================
//included dependences
#include "spring.h"

//=====================================
//bead filament class
class filament
{
    public:

        filament(filament_ensemble *net, vector<vector<double>> beadvec, double spring_length,
                double stretching_stiffnes, double ext, double bending_stiffness,
                double deltat, double temp, double fracture);

        ~filament();
    
        void update_d_strain(double);

        box *get_box();

        void pull_on_ends(double f);
        
        void affine_pull(double f);
        
        vector<filament *> update_stretching();
        
        void update_bending();
        
        void update_positions();

        spring *get_spring(int i);

        int get_nsprings();

        vector<vector<double>> output_beads(int fil);
        vector<vector<double>> output_springs(int fil);
        vector<double> output_thermo(int fil);

        string write_beads(int fil);
        
        string write_springs(int fil);
        
        string to_string();
        
        string write_thermo(int fil);
        
        double get_end2end();

        vector<vector<double>> get_beads(unsigned int first, unsigned int last);
        
        vector<filament *> fracture(int node);
        
        void update_forces(int index, vec_type f);
        
        bool operator==(const filament& that);
        
        void add_bead(vector<double> a, double l0, double kl, double me);

        inline double angle_between_springs(int i, int j);

        int get_nbeads();

        double get_bending_energy();

        virial_type get_bending_virial();

        double get_stretching_energy();

        double get_kinetic_energy_vel();

        double get_kinetic_energy_vir();

        virial_type get_stretching_virial();
        
        double get_potential_energy();
        
        double get_total_energy();
        
        void init_ubend();
    
        void print_thermo();
        
        void set_l0_max(double);
        
        void set_nsprings_max(int);
        
        void set_l0_min(double);
        
        void set_kgrow(double);
        
        void set_lgrow(double);

        vec_type get_bead_position(int bead);

        void update_length();

        void grow(double);

        void shrink(double);

    protected:
        box *bc;
        filament_ensemble *filament_network;

        double kb, temperature, dt, fracture_force, fracture_force_sq, ke_vir, ke_vel, damp, bd_prefactor, ubend;
        double spring_l0, l0_max, l0_min, kgrow, lgrow;
        int nsprings_max;

        virial_type bending_virial;
        vector<vec_type> prv_rnds;
        vector<class bead *> beads;
        vector<spring *> springs;
};

#endif
