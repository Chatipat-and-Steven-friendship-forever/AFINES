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
#ifndef __FILAMENT_ENSEMBLE_H_INCLUDED__
#define __FILAMENT_ENSEMBLE_H_INCLUDED__

//=====================================
// forward declared dependencies

//=====================================
//included dependences
#include "filament.h"
#include "box.h"
#include "quadrants.h"
//=====================================
//filament network class

int const CIRCLE = 1;

class filament_ensemble
{
    public:

        filament_ensemble(box *bc, vector< vector<double> > beads, array<int,2> mynq, double delta_t, double temp,
                double vis, double spring_len, double stretching, double ext, double bending, double frac_force);
        
        ~filament_ensemble();

        // quadrants
        quadrants *get_quads();
        void quad_update_serial();
        vector<array<int, 2>> *get_attach_list(double, double);

        vector<filament *> * get_network();

        filament * get_filament(int index);

        array<double,2> get_direction(int fil, int spring);

        array<double,2> get_start(int fil, int spring);
        
        array<double,2> get_end(int fil, int spring);
        
        array<double,2> get_force(int fil, int bead);
        
        double get_int_direction(int fil, int spring, double xp, double yp);

        double get_xcm(int fil, int spring);
       
        double get_ycm(int fil, int spring);

        double get_llength(int fil, int spring);
       
        double get_bead_friction();
        
        box *get_box();

        double get_stretching_energy();

        array<array<double, 2>, 2> get_stretching_virial();
        
        double get_bending_energy();

        array<array<double, 2>, 2> get_bending_virial();

        int get_nbeads();
        
        int get_nsprings();
        
        int get_nfilaments();

        void update_shear();
        
        void update_d_strain(double);
        
        void update_delrx(double);
        
        void update_stretching();
        
        void update_filament_stretching(int);
        
        void update_bending();
        
        void update_int_forces();

        void update_positions();

        void update_positions_range(int lo, int hi);
        
        void update_forces(int fil, int bead, double f2, double f3);

        vector<vector<double>> output_beads();
        vector<vector<double>> output_springs();
        vector<vector<double>> output_thermo();

        void write_beads(ofstream& fout);
        
        void write_springs(ofstream& fout);
        
        void write_thermo(ofstream& fout);
        
        void set_straight_filaments(bool is_straight);

        void set_y_thresh(double);
        
        void set_fene_dist_pct(double);
        
        void set_shear_rate(double);
        
        void set_shear_stop(double);

        void set_shear_dt(double);
        
        bool is_polymer_start(int f, int a);

        void set_visc(double v);

        vector<int> get_broken();

        void clear_broken();
        
        void print_filament_thermo();

        void print_network_thermo();

        void print_filament_lengths();
        
        void update();
        
        void update_energies();
        
        void set_circle_wall(double radius, double spring_constant);

        array<double, 2> external_force(array<double, 2> pos);

    protected:

        box *bc;

        int external_force_flag;
        double circle_wall_radius, circle_wall_spring_constant;

        double t, dt, temperature, visc;
        double gamma, shear_stop, shear_dt, shear_speed;
        double pe_stretch, pe_bend, ke;

        array<array<double, 2>, 2> vir_stretch, vir_bend;

        array<int, 2> nq;
        vector<int> broken_filaments;

        quadrants *quads;

        vector<filament *> network;
};

#endif
