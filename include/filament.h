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
#ifndef __FILAMENT_H_INCLUDED__
#define __FILAMENT_H_INCLUDED__

//=====================================
// forward declared dependencies
class filament_ensemble;

//=====================================
//included dependences
#include "bead.h"
#include "spring.h"

//=====================================
//bead filament class
class filament
{
    public:

        filament(array<double, 2> myfov, array<int, 2> mynq, double deltat, double temp, double shear, 
		 double frac, double bending_stiffness, string bndcnd, double drx);

        filament(filament_ensemble *net, array<double, 3> startpos, int nbead, array<int,2> mynq,
                double vis, double deltat, double temp, bool isStraight,
                double beadLength, double spring_length, double stretching, double ext, double bending, double fracture);

        filament(filament_ensemble *net, vector<bead *> beadvec, array<int, 2> mynq, double spring_length, 
                double stretching_stiffnes, double ext, double bending_stiffness, 
                double deltat, double temp, double fracture, double gamma);

        filament();
        
        ~filament();
    
        void set_y_thresh(double y);
        
        void set_shear(double g);

        void update_d_strain(double);

        void update_shear(double t); 

        void pull_on_ends(double f);
        
        void affine_pull(double f);
        
        vector<filament *> update_stretching(double t);
        
        void update_bending(double);
        
        void update_positions();
        
        void update_positions_range(int lo, int hi);
        
        array<double, 2> boundary_check(int i, double x, double y);
        
        void update(double t);
        
        bead * get_bead(int i);
        
        spring * get_spring(int i);

        int get_nsprings();

        vector<vector<array<int, 2> > > get_quadrants();
        //multimap<int, array<int, 2> > > get_quadrants();

        string write_beads(int fil);
        
        string write_springs(int fil);
        
        string to_string();
        
        string write_thermo(int fil);
        
        double get_end2end();

        vector<bead *> get_beads(unsigned int first, unsigned int last);
        
        vector<filament *> fracture(int node);
        
        void update_forces(int index, double f1, double f2);
        
        bool operator==(const filament& that);
        
        void add_bead(bead * a, double l0, double kl, double me);
        
        void set_BC(string s);

        string get_BC();

        array<double, 2> get_fov();

        double get_delrx();

        inline double angle_between_springs(int i, int j);

        void fwd_bwd_bending_update();
        
        void lammps_bending_update();
        
        int get_nbeads();

        double get_bending_energy();

        array<array<double, 2>, 2> get_bending_virial();

        double get_stretching_energy();

        double get_kinetic_energy();

        array<array<double, 2>, 2> get_stretching_virial();
        
        double get_potential_energy();
        
        double get_total_energy();
        
        void init_ubend();
    
        void print_thermo();
        
        array<double,2> get_bead_position(int bead);

    protected:

        filament_ensemble *filament_network;

        double kb, temperature, dt, fracture_force, fracture_force_sq, kinetic_energy, damp, kToverLp, bd_prefactor, ubend;
        double gamma, max_shear, y_thresh;

        array<array<double, 2>, 2> bending_virial;
        array<int,2> nq;
        vector<array<double, 2> > prv_rnds;
        vector<bead *> beads;
        vector<spring *> springs;
};

#endif
