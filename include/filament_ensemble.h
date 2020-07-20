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

        quadrants *get_quads();

        void quad_update_serial();

        vector<array<int, 2>> *get_attach_list(double, double);

        vector<filament *> * get_network();

        filament * get_filament(int index);

        array<double,2> get_direction(int fil, int spring);

        array<double,2> get_start(int fil, int spring);
        
        array<double,2> get_end(int fil, int spring);

        double get_int_direction(int fil, int spring, double xp, double yp);

        double get_llength(int fil, int spring);

        box *get_box();

        double get_stretching_energy();

        array<array<double, 2>, 2> get_stretching_virial();
        
        double get_bending_energy();

        array<array<double, 2>, 2> get_bending_virial();

        double get_exv_energy();

        double get_kinetic_energy_vel();

        double get_kinetic_energy_vir();

        int get_nbeads();

        int get_nsprings();

        int get_nfilaments();

        void update_d_strain(double);

        void update_stretching();

        void update_filament_stretching(int);

        void update_bending();

        void update_int_forces();

        void update_positions();

        void update_forces(int fil, int bead, double f2, double f3);

        vector<vector<double>> output_beads();

        vector<vector<double>> output_springs();

        vector<vector<double>> output_thermo();

        void update_spring_forces(int f); 

        void update_spring_forces_from_quads(); 

        void update_force_between_filaments(double n1, double l1, double n2, double l2); 

        void update_excluded_volume(int f); 

        void write_beads(ofstream& fout);
        
        void write_springs(ofstream& fout);
        
        void write_thermo(ofstream& fout);

        void set_fene_dist_pct(double);

        bool is_polymer_start(int f, int a);

        vector<int> get_broken();

        void clear_broken();
        
        void print_filament_thermo();

        void print_network_thermo();

        void print_filament_lengths();
        
        void update();
        
        void update_energies();
        
        void set_circle_wall(double radius, double spring_constant);

        array<double, 2> external_force(array<double, 2> pos);

        void set_growing(double, double, double, double, int);

    protected:

        box *bc;

        external_force_type external_force_flag;
        double circle_wall_radius, circle_wall_spring_constant;

        double t, dt;
        double temperature, spring_rest_len, visc, min_time;
        double rmax; 
        double kexv;  
        int nsprings_per_fil_max;

        double pe_stretch, pe_bend, pe_exv, ke_vel, ke_vir;
        array<array<double, 2>, 2> vir_stretch, vir_bend;

        vector<int> broken_filaments;

        quadrants *quads;

        vector<filament *> network;

};

#endif
