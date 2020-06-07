/*------------------------------------------------------------------
 filament_ensemble.cpp : container class for filaments
 
 Copyright (C) 2016 
 Created by: Simon Freedman, Shiladitya Banerjee, Glen Hocky, Aaron Dinner
 Contact: dinner@uchicago.edu

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version. See ../LICENSE for details. 
-------------------------------------------------------------------*/

#include "globals.h"
//#include "Link.h"
#include "filament_ensemble.h"

filament_ensemble::~filament_ensemble()
{
    cout << "DELETING FILAMENT_ENSEMBLE\n";
    delete quads;
    for (filament *f : network) delete f;
}

vector<filament *>* filament_ensemble::get_network()
{
    return &network;
}


filament * filament_ensemble::get_filament(int index)
{
    return network[index];
}

quadrants *filament_ensemble::get_quads()
{
    return quads;
}

void filament_ensemble::quad_update_serial()
{
    quads->clear();
    for (int f = 0; f < int(network.size()); f++) {
        for (int l = 0; l < network[f]->get_nsprings(); l++) {
            spring *s = network[f]->get_spring(l);
            quads->add_spring(s, {f, l});
        }
    }
}

vector<array<int, 2>> *filament_ensemble::get_attach_list(double x, double y)
{
    return quads->get_attach_list({x, y});
}

double filament_ensemble::get_llength(int fil, int spring)
{
    return network[fil]->get_spring(spring)->get_length();
}


array<double,2> filament_ensemble::get_start(int fil, int spring)
{
    return {{network[fil]->get_spring(spring)->get_hx()[0] , network[fil]->get_spring(spring)->get_hy()[0]}};
}


array<double,2> filament_ensemble::get_end(int fil, int spring)
{
    return {{network[fil]->get_spring(spring)->get_hx()[1] , network[fil]->get_spring(spring)->get_hy()[1]}};
}


array<double,2> filament_ensemble::get_force(int fil, int bead)
{
    return network[fil]->get_bead(bead)->get_force();
}


array<double,2> filament_ensemble::get_direction(int fil, int spring)
{
    return network[fil]->get_spring(spring)->get_direction();
}

void filament_ensemble::update_positions()
{
    int net_sz = int(network.size());
    for (int f = 0; f < net_sz; f++)
    {
        //if (f==0) cout<<"\nDEBUG: update_positions: using "<<omp_get_num_threads()<<" cores";  
        network[f]->update_positions();
    }

}

 
void filament_ensemble::update_positions_range(int lo, int hi)
{
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_positions_range(lo, hi);
    }

}

vector<vector<double>> filament_ensemble::output_beads()
{
    vector<vector<double>> out;
    for (size_t i = 0; i < network.size(); i++) {
        vector<vector<double>> tmp = network[i]->output_beads(i);
        out.insert(out.end(), tmp.begin(), tmp.end());
    }
    return out;
}

vector<vector<double>> filament_ensemble::output_springs()
{
    vector<vector<double>> out;
    for (size_t i = 0; i < network.size(); i++) {
        vector<vector<double>> tmp = network[i]->output_springs(i);
        out.insert(out.end(), tmp.begin(), tmp.end());
    }
    return out;
}

vector<vector<double>> filament_ensemble::output_thermo()
{
    vector<vector<double>> out;
    for (size_t i = 0; i < network.size(); i++) {
        out.push_back(network[i]->output_thermo(i));
    }
    return out;
}

void filament_ensemble::write_beads(ofstream& fout)
{
    for (unsigned int i=0; i<network.size(); i++) {
        fout<<network[i]->write_beads(i);
    } 
}

 
void filament_ensemble::write_springs(ofstream& fout)
{
    for (unsigned int i=0; i<network.size(); i++) {
        fout<<network[i]->write_springs(i);
    } 
}

 
void filament_ensemble::write_thermo(ofstream& fout){
    for (unsigned int f = 0; f < network.size(); f++)
        fout<<network[f]->write_thermo(f);
    
}

box *filament_ensemble::get_box()
{
    return bc;
}

void filament_ensemble::set_y_thresh(double y)
{
    for (unsigned int f = 0; f < network.size(); f++) network[f]->set_y_thresh(y);
}

void filament_ensemble::update_d_strain(double g)
{
    bc->update_d_strain(g);
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_d_strain(g);
    }
}

void filament_ensemble::print_filament_thermo(){
    
    for (unsigned int f = 0; f < network.size(); f++)
    {
        cout<<"\nF"<<f<<"\t:";
        network[f]->print_thermo();
    }

}

void filament_ensemble::update_energies()
{
    pe_stretch = 0.0;
    pe_bend = 0.0;
    ke = 0.0;
    virial_clear(vir_stretch);
    virial_clear(vir_bend);
    for (filament *f : network) {
        ke += f->get_kinetic_energy();
        pe_bend += f->get_bending_energy();
        pe_stretch += f->get_stretching_energy();
        virial_add(vir_stretch, f->get_stretching_virial());
        virial_add(vir_bend, f->get_bending_virial());
    }
}


double filament_ensemble::get_stretching_energy(){
    return pe_stretch;
}

 
double filament_ensemble::get_bending_energy(){
    return pe_bend;
}

array<array<double, 2>, 2> filament_ensemble::get_stretching_virial() {
    return vir_stretch;
}

array<array<double, 2>, 2> filament_ensemble::get_bending_virial() {
    return vir_bend;
}

void filament_ensemble::print_network_thermo(){
    cout<<"\nAll Fs\t:\tKE = "<<ke<<"\tPEs = "<<pe_stretch<<"\tPEb = "<<pe_bend<<"\tTE = "<<(ke+pe_stretch+pe_bend);
}

 
void filament_ensemble::print_filament_lengths(){
    for (unsigned int f = 0; f < network.size(); f++)
    {
        cout<<"\nF"<<f<<" : "<<network[f]->get_end2end()<<" um";
    }
}


 
bool filament_ensemble::is_polymer_start(int fil, int bead){

    return !(bead);

}

void filament_ensemble::update_forces(int f_index, int a_index, double f1, double f2){
    network[f_index]->update_forces(a_index, f1,f2);
}

 
vector<int> filament_ensemble::get_broken(){
    return broken_filaments;
}

 
void filament_ensemble::clear_broken(){
    broken_filaments.clear();
}

 
int filament_ensemble::get_nbeads(){
    int tot = 0;
    for (unsigned int f = 0; f < network.size(); f++)
        tot += network[f]->get_nbeads();
    return tot;
}

 
int filament_ensemble::get_nsprings(){
    return this->get_nbeads() - network.size();
}

 
int filament_ensemble::get_nfilaments(){
    return network.size();
}

double filament_ensemble::get_bead_friction(){
    
    if (network.size() > 0)
        if (network[0]->get_nbeads() > 0)
            return network[0]->get_bead(0)->get_friction();
    
    return 0;
}

// Update bending forces between monomers

void filament_ensemble::update_bending()
{
    int net_sz = int(network.size());
    
    for (int f = 0; f < net_sz; f++)
    {
        //if (f==0) cout<<"\nDEBUG: update_bending: using "<<omp_get_num_threads()<<" cores";  
        network[f]->update_bending(t);
    }
}


void filament_ensemble::update_stretching(){
    
//    vector<filament *> newfilaments;
    int s = network.size(); //keep it to one fracture per filament per timestep, or things get messy
    for (int f = 0; f < s; f++)
    {
        //if (f==0) cout<<"\nDEBUG: update_stretching: using "<<omp_get_num_threads()<<" cores";  
        this->update_filament_stretching(f);
    }
}


void filament_ensemble::update_filament_stretching(int f){
    vector<filament *> newfilaments = network[f]->update_stretching(t);

    if (newfilaments.size() > 0){ //fracture event occured

        cout<<"\n\tDEBUG: fracturing filament : "<<f;
        filament * broken = network[f];     //store a pointer to the broken filament to delete it with
        network[f] = newfilaments[0];       //replace that pointer with one of the new filaments

        if (newfilaments.size() == 2) network.push_back(newfilaments[1]); //add the second filament to the top of the stack

        broken_filaments.push_back(f);      // record the index, for automatic motor detachment
        delete broken;                      // delete the old filament

    }
}

void filament_ensemble::update_int_forces()
{
    this->update_stretching();
    this->update_bending();
}

/* Overdamped Langevin Dynamics Integrator (Leimkuhler, 2013) */

void filament_ensemble::update()
{
    for (int f = 0; f < int(network.size()); f++) {
        update_filament_stretching(f);
        network[f]->update_bending(t);
        if (external_force_flag != external_force_type::none) {
            for (int i = 0; i < network[f]->get_nbeads(); i++) {
                array<double, 2> pos = network[f]->get_bead_position(i);
                array<double, 2> force = external_force(pos);
                update_forces(f, i, force[0], force[1]);
            }
        }
        network[f]->update_positions();
    }
    update_energies();
    t += dt;
}

array<double, 2> filament_ensemble::external_force(array<double, 2> pos)
{
    if (external_force_flag == external_force_type::circle) {
        double x = pos[0];
        double y = pos[1];
        double rsq = x * x + y * y;
        if (rsq < circle_wall_radius * circle_wall_radius) {
            return {0, 0};
        }
        double r = sqrt(rsq);
        double k = -circle_wall_spring_constant * (1.0 - circle_wall_radius / r);
        return {k * x, k * y};
    } else {
        throw std::logic_error("External force flag not recognized.");
    }
}

void filament_ensemble::set_circle_wall(double radius, double spring_constant)
{
    external_force_flag = external_force_type::circle;
    circle_wall_radius = radius;
    circle_wall_spring_constant = spring_constant;
}

filament_ensemble::filament_ensemble(box *bc_, vector<vector<double> > beads, array<int,2> mynq, double delta_t, double temp,
        double vis, double spring_len, double stretching, double ext, double bending, double frac_force)
{
    external_force_flag = external_force_type::none;
    bc = bc_;

    dt = delta_t;
    t = 0;

    int fil_idx = 0;
    vector<bead *> avec;

    for (int i=0; i < int(beads.size()); i++){

        if (beads[i][3] != fil_idx && avec.size() > 0){

            network.push_back(new filament(this, avec, spring_len, stretching, ext, bending, delta_t, temp, frac_force));
            for (bead *b : avec) delete b;
            avec.clear();
            fil_idx = beads[i][3];
        }
        avec.push_back(new bead(beads[i][0], beads[i][1], beads[i][2], vis));
    }

    if (avec.size() > 0)
      network.push_back(new filament(this, avec, spring_len, stretching, ext, bending, delta_t, temp, frac_force));

    for (bead *b : avec) delete b;
    avec.clear();

    quads = new quadrants(bc, mynq);
    this->update_energies();
}
