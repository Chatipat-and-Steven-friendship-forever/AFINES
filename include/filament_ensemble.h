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
//=====================================
//filament network class

int const CIRCLE = 1;

class filament_ensemble
{
    public:

        filament_ensemble(box *bc, vector< vector<double> > beads, array<int,2> mynq, double delta_t, double temp,
                double vis, double spring_len, double stretching, double ext, double bending, double frac_force, bool check_dup_in_quad_);
        
        ~filament_ensemble();
        
        void nlist_init();
        
        void nlist_init_serial();
        
        void quad_update();
        
        void quad_update_serial();
        
        void consolidate_quads();

        void update_quads_per_filament(int);

        void reset_n_springs(int);

        void update_dist_map(set<pair<double, array<int, 2>>>& t_map, const array<int, 2>& mquad, double x, double y);

        vector<array<int, 2>> *get_attach_list(double, double);
        
        vector<filament *> * get_network();

        filament * get_filament(int index);

        set<pair<double, array<int,2>>> get_dist(double x, double y);
        
        set<pair<double, array<int,2>>> get_dist_all(double x, double y);
        
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
        array<int, 2> get_nq();
        
        double get_stretching_energy();

        array<array<double, 2>, 2> get_stretching_virial();
        
        double get_bending_energy();

        array<array<double, 2>, 2> get_bending_virial();

        int get_nbeads();
        
        int get_nsprings();
        
        int get_nfilaments();

        vector<vector<double> > spring_spring_intersections(double cllen, double prob);

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

        void set_fov(double x, double y);

        void set_nq(double x, double y);

        void set_visc(double v);

        vector<int> get_broken();

        void clear_broken();
        
        void print_filament_thermo();

        void print_network_thermo();

        void print_filament_lengths();
        
        void update();
        
        void update_energies();
        
        void turn_quads_off();
        
        void set_circle_wall(double radius, double spring_constant);

        array<double, 2> external_force(array<double, 2> pos);

    protected:

        box *bc;

        int external_force_flag;
        double circle_wall_radius, circle_wall_spring_constant;

        double t, dt, temperature, spring_rest_len, visc, min_time;
        double gamma, shear_stop, shear_dt, shear_speed;
        double max_springs_per_quad_per_filament, max_springs_per_quad; 
        bool straight_filaments = false, quad_off_flag;
        double pe_stretch, pe_bend, ke;

        array<array<double, 2>, 2> vir_stretch, vir_bend;

        array<double,2> view;
        array<int, 2> nq;
        vector<int> broken_filaments, empty_vector;
        
        vector< vector < vector< array<int, 2 > >* > * > springs_per_quad;
        vector< vector < int >* > n_springs_per_quad;
        vector<array<int, 2>> all_springs;

        bool check_dup_in_quad;

        vector<array<int, 2>* > all_quads;
        vector<filament *> network;
        unordered_set<array<int, 2>, boost::hash<array<int,2>>> fls;
};

#endif
