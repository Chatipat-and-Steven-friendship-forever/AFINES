/*------------------------------------------------------------------
 motor_ensemble.cpp : container class for motors

 Copyright (C) 2016
 Created by: Simon Freedman, Shiladitya Banerjee, Glen Hocky, Aaron Dinner
 Contact: dinner@uchicago.edu

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version. See ../LICENSE for details.
-------------------------------------------------------------------*/

#include "filament_ensemble.h"
#include "globals.h"
#include "motor_ensemble.h"
#include "motor.h"

motor_ensemble::motor_ensemble(vector<vector<double>> motors, double delta_t, double temp,
        double mlen, filament_ensemble *network, double v0, double stiffness, double max_ext_ratio,
        double ron, double roff, double rend, double fstall, double rcut, double vis)
{
    shear_flag = false;
    static_flag = false;
    f_network = network;
    network->get_box()->add_callback([this](double g) { this->update_d_strain(g); });
    mld = mlen;

    cout << "\nDEBUG: Number of motors:" << motors.size() << "\n";

    for (vector<double> mvec : motors) {
        n_motors.push_back(new motor(mvec, mlen, f_network, delta_t, temp,
                    v0, stiffness, max_ext_ratio, ron, roff, rend, fstall, rcut, vis));
    }

    this->update_energies();
}

motor_ensemble::~motor_ensemble()
{
    for (motor *m : n_motors) delete m;
}

void motor_ensemble::add_motor(motor *m)
{
    n_motors.push_back(m);
}

int motor_ensemble::get_nmotors()
{
    return n_motors.size();
}

motor *motor_ensemble::get_motor(int i)
{
    return n_motors[i];
}

// begin [settings]

void motor_ensemble::set_binding_two(double kon2, double koff2, double kend2){
    for(unsigned int i = 0; i < n_motors.size(); i++)
        n_motors[i]->set_binding_two(kon2, koff2, kend2);
}

void motor_ensemble::set_bending(double modulus, double ang){
    double kb  = modulus/mld;
    for(unsigned int i = 0; i < n_motors.size(); i++)
        n_motors[i]->set_bending(kb, ang);
}

void motor_ensemble::use_shear(bool flag)
{
    shear_flag = flag;
}

void motor_ensemble::use_static(bool flag)
{
    static_flag = flag;
}

void motor_ensemble::set_par(double k)
{
    for (motor *m : n_motors)
        m->set_par(k);
}

void motor_ensemble::set_antipar(double k)
{
    for (motor *m : n_motors)
        m->set_antipar(k);
}

void motor_ensemble::kill_heads(int hd)
{
    for (motor *m : n_motors) {
        m->kill_head(hd);
    }
}

void motor_ensemble::unbind_all_heads()
{
    for (motor *m : n_motors) {
        m->detach_head(0);
        m->detach_head(1);
        m->deactivate_head(0);
        m->deactivate_head(1);
    }
}

void motor_ensemble::revive_heads()
{
    for (motor *m : n_motors) {
        m->revive_head(0);
        m->revive_head(1);
    }
}

// end [settings]

// begin [dynamics]

void motor_ensemble::montecarlo()
{
    for (motor *m : n_motors) {
        array<motor_state, 2> s = m->get_states();

        mc_prob p;

        if (s[0] == motor_state::free) {
            m->try_attach(0, p);
        } else if (s[0] != motor_state::inactive) {
            m->try_detach(0, p);
        }

        if (s[1] == motor_state::free) {
            m->try_attach(1, p);
        } else if (s[1] != motor_state::inactive) {
            m->try_detach(1, p);
        }
    }
}

void motor_ensemble::integrate()
{
    for (motor *m : n_motors) {
        array<motor_state, 2> s = m->get_states();
        if (s[0] == motor_state::free || s[0] == motor_state::inactive) {
            m->brownian_relax(0);
        } else if (s[0] == motor_state::bound) {
            if (!static_flag) m->walk(0);
            m->update_position_attached(0);
        }
        if (s[1] == motor_state::free || s[1] == motor_state::inactive) {
            m->brownian_relax(1);
        } else if (s[1] == motor_state::bound) {
            if (!static_flag) m->walk(1);
            m->update_position_attached(1);
        }
        m->update_angle();
    }
}

void motor_ensemble::update_d_strain(double g)
{
    if (shear_flag) {
        for (motor *m : n_motors) {
            m->update_d_strain(g);
        }
    }
}

void motor_ensemble::compute_forces()
{
    for (motor *m : n_motors) {
        m->update_force();
        m->filament_update();
    }
    update_energies();
}

void motor_ensemble::update_energies()
{
    ke_vel = 0.0;
    ke_vir = 0.0;
    pe = 0.0;
    virial.zero();

    for (motor *m : n_motors) {
        ke_vel += m->get_kinetic_energy_vel();
        ke_vir += m->get_kinetic_energy_vir();
        pe += m->get_stretching_energy();
        virial += m->get_virial();
    }
}

// end [dynamics]

// begin [thermo]

double motor_ensemble::get_potential_energy()
{
    return pe;
}

virial_type motor_ensemble::get_virial()
{
    return virial;
}

// end [thermo]

// begin [output]

void motor_ensemble::print_ensemble_thermo()
{
    fmt::print(
            "\n"
            "All Motors\t:\t"
            "KEvel = {}\t"
            "KEvir = {}\t"
            "PEs = {}\t"
            "PEb = {}\t"
            "TE = {}",
            ke_vel, ke_vir, pe, 0, ke_vel + pe);
}

vector<vector<double>> motor_ensemble::output()
{
    vector<vector<double>> out;
    for (motor *m : n_motors) {
        out.push_back(m->output());
    }
    return out;
}

void motor_ensemble::motor_write(ostream& fout)
{
    for (motor *m : n_motors) {
        fout << m->write();
    }
}

void motor_ensemble::motor_write_doubly_bound(ostream& fout)
{
    array<motor_state, 2> doubly_bound = {motor_state::bound, motor_state::bound};
    for (size_t i = 0; i < n_motors.size(); i++) {
        if (n_motors[i]->get_states() == doubly_bound) {
            fmt::print(fout, "{}\t{}", n_motors[i]->write(), i);
        }
    }
}

// end [output]
