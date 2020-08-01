/*
 * motor.h
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#ifndef AFINES_MOTOR_ENSEMBLE_H
#define AFINES_MOTOR_ENSEMBLE_H

#include "motor.h"

class motor_ensemble
{
    public:

        motor_ensemble(vector<vector<double>> motors, double dt, double temp,
                double l0, filament_ensemble *network, double v0, double kl,
                double ron, double roff, double rend,
                double fstall, double rcut, double vis);
        ~motor_ensemble();

        void add_motor(motor *m);
        int get_nmotors();
        motor *get_motor(int);

        // [settings]
        void set_binding_two(double, double, double);
        void set_bending(double kb, double th0);
        void use_shear(bool flag);
        void use_static(bool flag);
        void set_par(double k);
        void set_antipar(double k);
        void set_align(double k);
        void kill_heads(int i);
        void unbind_all_heads();
        void revive_heads();
        void set_external(external *ext);

        // [dynamics]
        void montecarlo();  // attach/detach
        void integrate();  // brownian/walk
        void update_d_strain(double g);  // shear
        void compute_forces();  // compute force/energy/virial
        void update_energies();  // compute energy/virial

        // [thermo]
        // calculated by update_energies

        double get_potential_energy();
        double get_stretching_energy();
        double get_bending_energy();
        double get_alignment_energy();
        double get_external_energy();
        double get_binding_energy();

        virial_type get_potential_virial();
        virial_type get_stretching_virial();
        virial_type get_bending_virial();
        virial_type get_alignment_virial();
        virial_type get_external_virial();

        // [output]
        void print_ensemble_thermo();
        vector<vector<double>> output();
        void motor_write(ostream &fout);
        // void motor_write_doubly_bound(ostream &fout);
        // void motor_tension(ofstream& fout);

    protected:
        filament_ensemble *f_network;
        vector<motor *> n_motors;

        double mld;  // saved for set_bending
        external *ext;  // saved for deletion in destructor

        // flags
        bool shear_flag, static_flag;

        // thermo
        double pe_stretch, pe_bend, pe_align, pe_ext, pe_bind;
        virial_type vir_stretch, vir_bend, vir_align, vir_ext;
};

#endif
