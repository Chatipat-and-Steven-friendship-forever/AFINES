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
#include "motor.h"
#include "globals.h"

class filament
{
    public:
        filament(filament_ensemble *net, vector<vec_type> positions, vector<vec_type> randoms,
                double bead_radius, double visc,
                double l0, double kl, double kb, double dt, double temp, double fracture);
        ~filament();

        // [output]

        vector<vector<double>> output_beads(int fil);
        vector<vector<double>> output_springs(int fil);
        vector<double> output_thermo(int fil);

        string write_beads(int fil);
        string write_springs(int fil);
        string write_thermo(int fil);

        void print_thermo();

        // [thermo]

        double get_stretching_energy();
        double get_bending_energy();

        virial_type get_stretching_virial();
        virial_type get_bending_virial();

        // [state]

        int get_nbeads();
        vec_type get_bead_position(int bead);
        vec_type get_force(int i);

        int get_nsprings();
        spring *get_spring(int i);

        double get_end2end();

        // [dynamics]

        // updates bead and spring positions
        // computes kinetic energies
        // clears forces, but doesn't compute them
        void update_positions();

        // shears beads and springs
        void update_d_strain(double);

        // attempts to fracture the filament
        // failure: returns no filaments
        // success: returns two filaments, split at the first fracture site
        vector<filament *> try_fracture();
        vector<filament *> fracture(int node);  // helper method

        void detach_all_motors();

        void update_length();
        void grow(double);  // helper method

        // [internal forces]

        // updates spring forces
        // and applies them to beads
        void update_stretching();

        // computes and applies bending forces to beads
        // also computes bending energy and virial
        void update_bending();

        // [external forces]

        // applies forces to beads
        void pull_on_ends(double f);
        void affine_pull(double f);
        void update_forces(int index, vec_type f);

        // [misc]

        box *get_box();

        void add_bead(vec_type a, double l0, double kl);

        bool operator==(const filament& that);

        // [attached positions]

        int new_attached(motor *m, int hd, int l, vec_type pos);
        void del_attached(int i);

        int get_attached_l(int i);
        vec_type get_attached_pos(int i);
        void add_attached_force(int i, vec_type f);
        void add_attached_pos(int i, double dist);

        bool at_barbed_end(int i);
        bool at_pointed_end(int i);

        // [growing settings]
        void set_l0_max(double);
        void set_nsprings_max(int);
        void set_l0_min(double);
        void set_kgrow(double);
        void set_lgrow(double);

    protected:
        box *bc;
        filament_ensemble *filament_network;

        // state

        vector<vec_type> positions;
        vector<vec_type> forces;
        vector<vec_type> prv_rnds;

        vector<spring *> springs;

        struct attached_type { class motor *m; int hd; int l; double pos; };
        vector<attached_type> attached;

        // thermo
        double ubend;
        virial_type bending_virial;

        // parameters
        double kb, temperature, dt, fracture_force, damp, bd_prefactor;
        double bead_radius, visc;

        // growing parameters
        int nsprings_max;
        double spring_l0, l0_max, l0_min, kgrow, lgrow;
};

#endif
