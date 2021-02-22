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
#include "ext.h"

class filament_ensemble
{
    public:

        filament_ensemble(box *bc, vector<vector<double>> beads, array<int, 2> nq,
                double dt, double temp, double visc, double l0, double kl, double kb,
                double frac_force, double rmax, double kexv, double shake_tol, int max_shake);

        ~filament_ensemble();

        // settings
        void set_external(external *);

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

        // attached locations
        fp_index_type new_attached(motor *m, int hd, int f_index, int l_index, vec_type pos);
        void del_attached(fp_index_type i);
        array<int, 2> get_attached_fl(fp_index_type i);
        vec_type get_attached_pos(fp_index_type i);
        void add_attached_force(fp_index_type i, vec_type f);
        void add_attached_pos(fp_index_type i, double dist);
        vec_type get_attached_direction(fp_index_type i);
        bool at_pointed_end(fp_index_type i);
        bool at_barbed_end(fp_index_type i);

        // thermo
        double get_potential_energy();
        double get_stretching_energy();
        double get_bending_energy();
        double get_excluded_energy();
        double get_external_energy();

        virial_type get_potential_virial();
        virial_type get_stretching_virial();
        virial_type get_bending_virial();
        virial_type get_excluded_virial();
        virial_type get_external_virial();

        // monte carlo
        void montecarlo();

        void set_growing(double, double, double, double, int);
        void try_grow();
        void try_fracture();

        // dynamics
        void integrate();
        void update_d_strain(double);

        // update forces/energies
        void compute_forces();

        void update_stretching();
        void update_bending();
        void update_excluded_volume();
        void update_external();
        void update_energies();

        void update_forces(int fil, int bead, vec_type f);

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
        external *ext;
        vector<filament *> network;

        // thermo
        double pe_stretch, pe_bend, pe_exv, pe_ext;
        virial_type vir_stretch, vir_bend, vir_exv, vir_ext;
};

#endif
