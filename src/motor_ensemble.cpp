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
#include "spacer.h"

template <class motor_type>
motor_ensemble<motor_type>::motor_ensemble(vector<vector<double>> motors, double delta_t, double temp,
        double mlen, filament_ensemble *network, double v0, double stiffness, double max_ext_ratio,
        double ron, double roff, double rend, double fstall, double rcut, double vis)
{
    attach_opt_flag = false;
    shear_flag = false;
    static_flag = false;
    f_network = network;
    network->get_box()->add_callback([this](double g) { this->update_d_strain(g); });
    mld = mlen;

    cout << "\nDEBUG: Number of motors:" << motors.size() << "\n";

    for (vector<double> mvec : motors) {
        n_motors.push_back(new motor_type(mvec, mlen, f_network, delta_t, temp,
                    v0, stiffness, max_ext_ratio, ron, roff, rend, fstall, rcut, vis));
    }

    this->update_energies();
}

template <class motor_type>
motor_ensemble<motor_type>::~motor_ensemble()
{
    cout << "DELETING MOTOR ENSEMBLE\n";
    for (motor *m : n_motors) delete m;
}

template <class motor_type>
int motor_ensemble<motor_type>::get_nmotors()
{
    return n_motors.size();
}

template <class motor_type>
motor *motor_ensemble<motor_type>::get_motor(int i)
{
    return n_motors[i];
}

template <class motor_type>
void motor_ensemble<motor_type>::kill_heads(int hd)
{
    for (motor *m : n_motors) {
        m->kill_head(hd);
    }
}

//check if any motors attached to filaments that no longer exist;
// if they do, detach them
// Worst case scenario, this is a O(m*n*p) function,
// where m is the number of filaments
//       n is the number of rods per filament
//       p is the number of motors
// However: we don't expect to fracture often,
// so this loop should rarely if ever be accessed.
template <class motor_type>
void motor_ensemble<motor_type>::check_broken_filaments()
{
    for (int i : f_network->get_broken()) {
        for (motor *m : n_motors) {
            array<int, 2> f_index = m->get_f_index();

            if (f_index[0] == i) {
                m->detach_head_without_moving(0);
                //cout<<"\nDEBUG: detaching head 0 of motor "<<j<<" from broken filament "<<broken_filaments[i];
            }

            if (f_index[1] == i) {
                m->detach_head_without_moving(1);
                //cout<<"\nDEBUG: detaching head 1 of motor "<<j<<" from broken filament "<<broken_filaments[i];
            }
        }
    }
}

template <class motor_type>
void motor_ensemble<motor_type>::motor_update()
{
    this->check_broken_filaments();

    for (motor *m : n_motors) {
        if (static_flag) {
            // Used for static, constantly attached, motors -- ASSUMES both heads are ALWAYS attached
            m->update_position_attached(0);
            m->update_position_attached(1);
            m->update_angle();
            m->update_force();
            m->filament_update();

        } else {
            array<motor_state, 2> s = m->get_states();

            // Dynamics
            if (s[0] == motor_state::free || s[0] == motor_state::inactive) {
                m->brownian_relax(0);
            }
            if (s[1] == motor_state::free || s[1] == motor_state::inactive) {
                m->brownian_relax(1);
            }

            m->update_angle();
            m->update_force();
            m->filament_update();

            // Attachment or Movement + Detachment
            if (s[0] == motor_state::free) {
                if (attach_opt_flag) {
                    m->attach_opt(0);
                } else {
                    m->attach(0);
                }
            } else if (s[0] != motor_state::inactive) {
                m->step_onehead(0);
            }

            if (s[1] == motor_state::free) {
                if (attach_opt_flag) {
                    m->attach_opt(1);
                } else {
                    m->attach(1);
                }
            } else if (s[1] != motor_state::inactive) {
                m->step_onehead(1);
            }

            /*
            if (s[0] == motor_state::bound) {
                m->step_onehead(0);
            } else if (s[0] == motor_state::free) {
                bool attached = m->attach(0);
                if (!attached) m->brownian_relax(0);
            }

            if (s[1] == motor_state::bound) {
                m->step_onehead(1);
            } else if (s[1] == motor_state::free) {
                bool attached = m->attach(1);
                if (!attached) m->brownian_relax(1);
            }

            m->update_angle();
            m->update_force();
            m->filament_update();
            */

        }
    }
    this->update_energies();
}

template <class motor_type>
vector<vector<double>> motor_ensemble<motor_type>::output()
{
    vector<vector<double>> out;
    for (motor *m : n_motors) {
        out.push_back(m->output());
    }
    return out;
}

template <class motor_type>
void motor_ensemble<motor_type>::motor_write(ostream& fout)
{
    for (motor *m : n_motors) {
        fout << m->write();
    }
}

template <class motor_type>
void motor_ensemble<motor_type>::add_motor(motor_type *m)
{
    n_motors.push_back(m);
}

template <class motor_type>
void motor_ensemble<motor_type>::update_energies()
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

template <class motor_type>
double motor_ensemble<motor_type>::get_potential_energy()
{
    return pe;
}

template <class motor_type>
virial_type motor_ensemble<motor_type>::get_virial()
{
    return virial;
}

template <class motor_type>
void motor_ensemble<motor_type>::print_ensemble_thermo()
{
    cout << "\nAll Motors\t:\tKE = " << ke_vel << "\tKEvir" << ke_vir << "\tPEs = " << pe << "\tPEb = " << 0 << "\tTE = " << (ke_vel + pe);
}

template <class motor_type>
void motor_ensemble<motor_type>::unbind_all_heads()
{
    for (motor *m : n_motors) {
        m->detach_head(0);
        m->detach_head(1);
        m->deactivate_head(0);
        m->deactivate_head(1);
    }
}

template <class motor_type>
void motor_ensemble<motor_type>::revive_heads()
{
    for (motor *m : n_motors) {
        m->revive_head(0);
        m->revive_head(1);
    }
}

template <class motor_type>
void motor_ensemble<motor_type>::motor_write_doubly_bound(ostream& fout)
{
    array<motor_state, 2> doubly_bound = {motor_state::bound, motor_state::bound};

    for (unsigned int i=0; i<n_motors.size(); i++) {
        if (n_motors[i]->get_states() == doubly_bound){
            fout<<n_motors[i]->write();
            fout<<"\t"<<i;
        }
    }
}

template <class motor_type>
void motor_ensemble<motor_type>::set_binding_two(double kon2, double koff2, double kend2){
    for(unsigned int i = 0; i < n_motors.size(); i++)
        n_motors[i]->set_binding_two(kon2, koff2, kend2);
}

void spacer_ensemble::set_bending(double modulus, double ang){
    double kb  = modulus/mld;
    for(unsigned int i = 0; i < n_motors.size(); i++)
        n_motors[i]->set_bending(kb, ang);
}

template <class motor_type>
void motor_ensemble<motor_type>::use_attach_opt(bool flag)
{
    attach_opt_flag = flag;
}

template <class motor_type>
void motor_ensemble<motor_type>::use_shear(bool flag)
{
    shear_flag = flag;
}

template <class motor_type>
void motor_ensemble<motor_type>::use_static(bool flag)
{
    static_flag = flag;
}

template <class motor_type>
void motor_ensemble<motor_type>::update_d_strain(double g)
{
    if (shear_flag) {
        for (motor *m : n_motors) {
            m->update_d_strain(g);
        }
    }
}

template class motor_ensemble<motor>;
template class motor_ensemble<spacer>;
