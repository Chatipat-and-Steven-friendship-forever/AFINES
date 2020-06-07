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

#include "globals.h"
#include "motor.h"
#include "motor_ensemble.h"
#include "filament_ensemble.h"

motor_ensemble::motor_ensemble(vector<vector<double>> motors, double delta_t, double temp,
        double mlen, filament_ensemble *network, double v0, double stiffness, double max_ext_ratio,
        double ron, double roff, double rend,
        double fstall, double rcut,
        double vis, bool use_attach_opt_)
{
    use_attach_opt = use_attach_opt_;
    f_network=network;

    ke = 0;
    pe = 0;

    cout << "\nDEBUG: Number of motors:" << motors.size() << "\n";

    for (vector<double> mvec : motors) {
        n_motors.push_back(new motor(mvec, mlen, f_network, delta_t, temp,
                    v0, stiffness, max_ext_ratio, ron, roff, rend, fstall, rcut, vis));
    }

    this->update_energies();
}


motor_ensemble::~motor_ensemble()
{
    cout << "DELETING MOTOR ENSEMBLE\n";
    for (motor *m : n_motors) delete m;
}

int motor_ensemble::get_nmotors()
{
    return n_motors.size();
}

motor *motor_ensemble::get_motor(int i)
{
    return n_motors[i];
}

void motor_ensemble::kill_heads(int hd){
    for (unsigned int i = 0; i < n_motors.size(); i++)
        n_motors[i]->kill_head(hd);
}

//check if any motors attached to filaments that no longer exist; 
// if they do, detach them
// Worst case scenario, this is a O(m*n*p) function,
// where m is the number of filaments
//       n is the number of rods per filament
//       p is the number of motors
// However: we don't expect to fracture often, 
// so this loop should rarely if ever be accessed.
    


void motor_ensemble::check_broken_filaments()
{
    vector<int> broken_filaments = f_network->get_broken();
    array<int, 2> f_index;
    
    for (unsigned int i = 0; i < broken_filaments.size(); i++){
        
        for(unsigned int j = 0; j < n_motors.size(); j++){
            
            f_index = n_motors[j]->get_f_index();

            if(f_index[0] == broken_filaments[i]){
                n_motors[j]->detach_head_without_moving(0);
                //cout<<"\nDEBUG: detaching head 0 of motor "<<j<<" from broken filament "<<broken_filaments[i];
            }

            if(f_index[1] == broken_filaments[i]){
                n_motors[j]->detach_head_without_moving(1);
                //cout<<"\nDEBUG: detaching head 1 of motor "<<j<<" from broken filament "<<broken_filaments[i];
            }
        }
    }


}


void motor_ensemble::motor_walk(double t)
{

    this->check_broken_filaments();
    int nmotors_sz = int(n_motors.size());
    //#pragma omp parallel for

    for (int i=0; i<nmotors_sz; i++) {

        array<motor_state, 2> s = n_motors[i]->get_states();

        if (t >= 0.0) {

            //Dynamics
            if (s[0] == motor_state::free || s[0] == motor_state::inactive)
                n_motors[i]->brownian_relax(0);
            if (s[1] == motor_state::free || s[1] == motor_state::inactive)
                n_motors[i]->brownian_relax(1);

            n_motors[i]->update_angle();
            n_motors[i]->update_force();
            n_motors[i]->filament_update();

            //Attachment or Movement + Detachment
            if (s[0] == motor_state::free) {
                if (use_attach_opt) {
                    n_motors[i]->attach_opt(0);
                } else {
                    n_motors[i]->attach(0);
                }
            } else if (s[0] != motor_state::inactive) {
                n_motors[i]->step_onehead(0);
            }

            if (s[1] == motor_state::free) {
                if (use_attach_opt) {
                    n_motors[i]->attach_opt(1);
                } else {
                    n_motors[i]->attach(1);
                }
            } else  if (s[1] != motor_state::inactive) {
                n_motors[i]->step_onehead(1);
            }

        }
    
    }
    this->update_energies();
    
}

/* Used for static, contstantly attached, motors -- ASSUMES both heads are ALWAYS attached */

void motor_ensemble::motor_update()
{

    this->check_broken_filaments();
    int nmotors_sz = int(n_motors.size());
    //#pragma omp parallel for
    
    for (int i=0; i<nmotors_sz; i++) {
       
            n_motors[i]->update_position_attached(0);
            n_motors[i]->update_position_attached(1);
            n_motors[i]->update_angle();
            n_motors[i]->update_force();
            //n_motors[i]->update_force_fraenkel_fene();
            n_motors[i]->filament_update();
    
    }
    this->update_energies();
    
}

vector<vector<double>> motor_ensemble::output()
{
    vector<vector<double>> out;
    for (size_t i = 0; i < n_motors.size(); i++) {
        out.push_back(n_motors[i]->output());
    }
    return out;
}

void motor_ensemble::motor_write(ostream& fout)
{
    for (unsigned int i=0; i<n_motors.size(); i++) {
        fout<<n_motors[i]->write();
    } 
}


void motor_ensemble::add_motor(motor * m)
{
    n_motors.push_back(m);
}

void motor_ensemble::update_energies()
{
    ke = 0;
    pe = 0;
    virial[0][0] = virial[0][1] = 0.0;
    virial[1][0] = virial[1][1] = 0.0;
   
    for (unsigned int m = 0; m < n_motors.size(); m++)
    {
        ke += n_motors[m]->get_kinetic_energy();
        pe += n_motors[m]->get_stretching_energy();
        //pe += n_motors[m]->get_stretching_energy_fene();
        array<array<double, 2>, 2> vir = n_motors[m]->get_virial();
        virial[0][0] += vir[0][0]; virial[0][1] += vir[0][1];
        virial[1][0] += vir[1][0]; virial[1][1] += vir[1][1];
    }
}

 
double motor_ensemble::get_potential_energy(){
    return pe;
}

array<array<double, 2>, 2> motor_ensemble::get_virial() {
    return virial;
}

void motor_ensemble::print_ensemble_thermo(){
    cout<<"\nAll Motors\t:\tKE = "<<ke<<"\tPEs = "<<pe<<"\tPEb = "<<0<<"\tTE = "<<(ke+pe);
}


void motor_ensemble::unbind_all_heads()
{
    for (int i = 0; i < n_motors.size(); i++) {
        n_motors[i]->detach_head(0);
        n_motors[i]->detach_head(1);
        n_motors[i]->deactivate_head(0);
        n_motors[i]->deactivate_head(1);
    }
}

void motor_ensemble::revive_heads()
{
    for (int i = 0; i < n_motors.size(); i++) {
        n_motors[i]->revive_head(0);
        n_motors[i]->revive_head(1);
    }
}

void motor_ensemble::update_d_strain(double g)
{
    for (motor *m : n_motors) {
        m->update_d_strain(g);
    }
}
