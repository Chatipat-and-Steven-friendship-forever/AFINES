/*-------------------------------------------------------------------
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

vector<array<int, 2>> *filament_ensemble::get_attach_list(vec_type pos)
{
    return quads->get_attach_list(pos);
}

double filament_ensemble::get_llength(int fil, int spring)
{
    return network[fil]->get_spring(spring)->get_length();
}


array<double,2> filament_ensemble::get_start(int fil, int spring)
{
    return network[fil]->get_spring(spring)->get_h0();
}


array<double,2> filament_ensemble::get_end(int fil, int spring)
{
    return network[fil]->get_spring(spring)->get_h1();
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

void filament_ensemble::update_d_strain(double g)
{
    for (filament *f : network) {
        f->update_d_strain(g);
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
    ke_vel = 0.0;
    ke_vir = 0.0;
    virial_clear(vir_stretch);
    virial_clear(vir_bend);
    for (filament *f : network) {
        ke_vel += f->get_kinetic_energy_vel();
        ke_vir += f->get_kinetic_energy_vir();
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

double filament_ensemble::get_kinetic_energy_vir(){
    return ke_vir;
}

void filament_ensemble::print_network_thermo(){
    cout<<"\nAll Fs\t:\tKE = "<<ke_vir<<"\tPEs = "<<pe_stretch<<"\tPEb = "<<pe_bend<<"\tPEexv = "<<pe_exv<<"\tTE = "<<(ke_vir+pe_stretch+pe_bend+pe_exv);
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

void filament_ensemble::update_forces(int f_index, int a_index, vec_type f){
    network[f_index]->update_forces(a_index, f);
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
    pe_exv = 0;

    if ( kexv > 0)
        this->update_spring_forces_from_quads();
    for (int f = 0; f < int(network.size()); f++) {
        network[f]->update_length();
        update_filament_stretching(f);
        network[f]->update_bending(t);
        if (external_force_flag != external_force_type::none) {
            for (int i = 0; i < network[f]->get_nbeads(); i++) {
                array<double, 2> pos = network[f]->get_bead_position(i);
                array<double, 2> force = external_force(pos);
                update_forces(f, i, force);
            }
        }
        network[f]->update_positions();
    }
    this->update_energies();
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

void filament_ensemble::update_spring_forces_from_quads()
{
    //This function loops through the quads of the systems and then loops through the filaments and springs described by the neighbor list in every quad.
    //Upon this looping, the pair interactions will be calculated according to the neighbor list.
    //The loop will go through thee values of spring_per_quad() [x][y][i], where x: 0-nq[0], y: 0-nq[1], i: 0-m_springs_per_quad
    //The pairs found in he nieghbor list are then saved in a vector array.
    //On subsequent loops, this value will be searched for in order to ensure no repeats in the force calculation.

    int max_nsprings = network.size()*nsprings_per_fil_max;

    // track which excluded volume forces have already been calculated
    vector<vector<int>> int_lks(max_nsprings, vector<int>(max_nsprings, 0));

    // for each quadrant
    array<int, 2> nq = quads->get_nq();
    for (int x = 0; x < nq[0]; x++) {
        for (int y = 0; y < nq[1]; y++) {
            vector<array<int, 2>> *q = quads->get_quad({x, y});

            // for each pair of springs in the quadrant
            int nsprings_at_quad = q->size();
            for (int i = 0; i < nsprings_at_quad; i++) {
                array<int, 2> spring_1 = q->at(i);
                int f1 = spring_1[0];
                int l1 = spring_1[1];
                for (int j = i+1; j < nsprings_at_quad; j++) {
                    array<int, 2> spring_2 = q->at(j);
                    int f2 = spring_2[0];
                    int l2 = spring_2[1];

                    // adjacent springs would yield excluded volume interactions between the same bead
                    if (f1 == f2 && abs(l1 - l2) < 2) continue;

                    // flattened index of springs
                    int par1 = f1*nsprings_per_fil_max + l1;
                    int par2 = f2*nsprings_per_fil_max + l2;

                    // only compute excluded volume forces if it hasn't been computed already
                    if (!int_lks[par1][par2]) {
                        int_lks[par1][par2] = 1;
                        int_lks[par2][par1] = 1;
                        update_force_between_filaments(f1, l1, f2, l2);
                    }
                }
            }
        }
    }
}

void filament_ensemble::update_spring_forces(int f)
{
    //This function loops through every filament and spring in the network and applies the force calulation under certain limits
    int net_sz = network.size();

    // for every spring in filament f
    int lks_sz = network[f]->get_nsprings();
    for (int i = 0; i < lks_sz; i++) {

        // for every spring in filaments g > f
        for (int g = f+1; g < net_sz; g++) {
            int oth_lks_sz = network[g]->get_nsprings();
            for (int j = 0; j < oth_lks_sz; j++) {

                update_force_between_filaments(f, i, g, j);
            }
        }
    }
}

void filament_ensemble::update_force_between_filaments(double n1, double l1, double n2, double l2)
{
    //This function calculates the forces applied to the actin beads of a pair of filaments under certain limits.
    //Here, we use distance of closest approach to describe the direction and magnitude of the forces.

    double b = 1/rmax;

    spring *s1 = network[n1]->get_spring(l1);
    spring *s2 = network[n2]->get_spring(l2);

    vec_type h0_1 = s1->get_h0();
    vec_type h1_1 = s1->get_h1();

    vec_type h0_2 = s2->get_h0();
    vec_type h1_2 = s2->get_h1();

    array<double, 2> len;
    len[0] = s1->get_length();
    len[1] = s2->get_length();

    // compute the nearest point on the other spring
    vec_type p1 = s1->intpoint(h0_2);
    vec_type p2 = s1->intpoint(h1_2);
    vec_type p3 = s2->intpoint(h0_1);
    vec_type p4 = s2->intpoint(h1_1);

    // compute the distance to the nearest point
    array<double, 4> r_c;
    r_c[0] = bc->dist_bc(h0_2 - p1);
    r_c[1] = bc->dist_bc(h1_2 - p2);
    r_c[2] = bc->dist_bc(h0_1 - p3);
    r_c[3] = bc->dist_bc(h1_1 - p4);

    // find the minimum distance between the filaments
    // assuming that they don't intersect
    double r = r_c[0];
    int index = 0;
    for(int k = 1; k < 4; k++){
        if(r_c[k] < r){
            r = r_c[k];
            index = k;
        }
    }

    bool intersect = s1->get_line_intersect(s2);

    // assume spring length l0 < 2 rmax ?
    if (r < rmax) {

        if (!intersect) {
            // doesn't intersect

            double length=0, len1=0;
            vec_type dist;
            if (index == 0) {
                r = r_c[0];
                len1 = bc->dist_bc(h0_1 - p1);
                length = len[0];
                dist = bc->rij_bc(p1 - h0_2);
            } else if (index == 1) {
                r = r_c[1];
                len1 = bc->dist_bc(h0_1 - p2);
                length = len[0];
                dist = bc->rij_bc(p2 - h1_2);
            } else if (index == 2) {
                r = r_c[2];
                len1 = bc->dist_bc(h0_2 - p3);
                length = len[1];
                dist = bc->rij_bc(p3 - h0_1);
            } else if (index == 3) {
                r = r_c[3];
                len1 = bc->dist_bc(h0_2 - p4);
                length = len[1];
                dist = bc->rij_bc(p4 - h1_1);
            }

            double len2 = length - len1;
            double r_1 = len2/length;
            double r_2 = len1/length;

            vec_type F = 2*kexv*dist*b*((1/r) - b);

            pe_exv += kexv*pow((1-r*b),2);

            if (index == 0) {
                network[n1]->update_forces(l1, F*r_1);
                network[n1]->update_forces(l1+1, F*r_2);
                network[n2]->update_forces(l2, -F);
            } else if (index == 1) {
                network[n1]->update_forces(l1, F*r_1);
                network[n1]->update_forces(l1+1, F*r_2);
                network[n2]->update_forces(l2+1, -F);
            } else if (index == 2) {
                network[n2]->update_forces(l2, F*r_1);
                network[n2]->update_forces(l2+1, F*r_2);
                network[n1]->update_forces(l1, -F);
            } else if (index == 3) {
                network[n2]->update_forces(l2, F*r_1);
                network[n2]->update_forces(l2+1, F*r_2);
                network[n1]->update_forces(l1+1, -F);
            }

        } else {
            // intersects

            // apply constant force?
            // this doesn't look right
            vec_type F = {2*kexv/(rmax*sqrt(2)), 2*kexv/(rmax*sqrt(2))};

            pe_exv += kexv*pow((1-r*b),2);

            network[n1]->update_forces(l1, F);
            network[n1]->update_forces(l1+1, F);
            network[n2]->update_forces(l2, -F);
            network[n2]->update_forces(l2+1, -F);

        }
    }
}

double filament_ensemble::get_exv_energy()
{
    return pe_exv;
}

void filament_ensemble::update_excluded_volume(int f)
{
    //For every filament bead on f, for every bead not on f, calculate the force between the two bead using the Jones potential, and update them ( maybe divide by half due to overcaluclations).

    int net_sz = network.size();
    //10^6 included to account for m to microm conversion
    double a = 0.004;
    double b = 1/rmax;

    // for every bead in filament f
    int act_sz = network[f]->get_nbeads();
    for (int i = 0; i < act_sz; i++) {
        array<double, 2> h1 = network[f]->get_bead_position(i);
        double x1 = h1[0];
        double y1 = h1[1];

        // for every bead in filament g > f
        for (int g = f+1; g < net_sz; g++) {
            int act_sz_other = network[g]->get_nbeads();
            for (int j = 0; j < act_sz_other; j++) {
                array<double, 2> h2 = network[g]->get_bead_position(j);
                double x2 = h2[0];
                double y2 = h2[1];

                // compute forces for the potential
                // U(r) = a (b r - 1)^2 if r < 1 / b
                double dx = x1 - x2;
                double dy = y1 - y2;

                double r = bc->dist_bc({dx, dy});
                if (r > 0.0 && r <= rmax) {
                    double Fx = 2*dx*a*b*((1/r)-b);
                    double Fy = 2*dy*a*b*((1/r)-b);
                    vec_type F = {Fx, Fy};

                    network[f]->update_forces(i,F);
                    network[g]->update_forces(j,-F);
                }
            }
        }
    }
}

filament_ensemble::filament_ensemble(box *bc_, vector<vector<double> > beads, array<int,2> mynq, double delta_t, double temp,
        double vis, double spring_len, double stretching, double ext, double bending, double frac_force, double RMAX, double A)
{
    external_force_flag = external_force_type::none;
    bc = bc_;

    bc->add_callback([this](double g) { this->update_d_strain(g); });

    rmax = RMAX;
    kexv = A;
    dt = delta_t;
    t = 0;

    int fil_idx = 0;
    vector<vector<double>> avec;

    nsprings_per_fil_max = 0;
    for (int i=0; i < int(beads.size()); i++){

        if (beads[i][3] != fil_idx && avec.size() > 0){

            network.push_back(new filament(this, avec, spring_len, stretching, ext, bending, delta_t, temp, frac_force));
            nsprings_per_fil_max = max(nsprings_per_fil_max, int(avec.size() - 1));
            avec.clear();
            fil_idx = beads[i][3];
        }
        avec.push_back({beads[i][0], beads[i][1], beads[i][2], vis});
    }

    if (avec.size() > 0)
      network.push_back(new filament(this, avec, spring_len, stretching, ext, bending, delta_t, temp, frac_force));

    avec.clear();

    quads = new quadrants(bc, mynq);
    this->update_energies();

    pe_stretch = 0;
    pe_bend = 0;
    pe_exv = 0;
    ke_vir = 0;
}

void filament_ensemble::set_growing(double kgrow, double lgrow, double l0min, double l0max, int nsprings_max)
{
    nsprings_per_fil_max = nsprings_max;
    for (int i = 0; i < int(network.size()); i++){
        network[i]->set_kgrow(kgrow);
        network[i]->set_lgrow(lgrow);
        network[i]->set_l0_min(l0min);
        network[i]->set_l0_max(l0max);
        network[i]->set_nsprings_max(nsprings_max);
    }
}
