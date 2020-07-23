/*------------------------------------------------------------------
 filament.cpp : object describing a worm-like chain filament

 Copyright (C) 2016
 Created by: Simon Freedman, Shiladitya Banerjee, Glen Hocky, Aaron Dinner
 Contact: dinner@uchicago.edu

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version. See ../LICENSE for details.
-------------------------------------------------------------------*/

#include "filament.h"
#include "filament_ensemble.h"
#include "bead.h"
#include "globals.h"
#include "potentials.h"

filament::filament(filament_ensemble *net, vector<vector<double>> beadvec, double spring_length,
        double stretching_stiffness, double max_ext_ratio, double bending_stiffness,
        double deltat, double temp, double frac_force)
{
    filament_network = net;

    bc = net->get_box();
    dt = deltat;
    temperature = temp;
    fracture_force = frac_force;
    fracture_force_sq = fracture_force*fracture_force;
    kb = bending_stiffness;

    ke_vel = 0.0;
    ke_vir = 0.0;

    spring_l0 = spring_length;

    nsprings_max = 0;
    l0_max = 0.0;
    l0_min = 0.0;
    kgrow = 0.0;
    lgrow = 0.0;

    if (beadvec.size() > 0)
    {
        vector<double> &entry = beadvec[0];
        if (entry.size() != 4) throw runtime_error("Wrong number of arguments in beadvec.");
        beads.push_back(new bead(entry[0], entry[1], entry[2], entry[3]));
        prv_rnds.push_back({});
        damp = beads[0]->get_friction();
    }

    //spring em up
    if (beadvec.size() > 1){
        for (unsigned int j = 1; j < beadvec.size(); j++) {

            vector<double> &entry = beadvec[j];
            if (entry.size() != 4) throw runtime_error("Wrong number of arguments in beadvec.");
            beads.push_back(new bead(entry[0], entry[1], entry[2], entry[3]));
            springs.push_back(new spring(spring_length, stretching_stiffness, max_ext_ratio, this, {(int)j-1, (int)j}));
            springs[j-1]->step();
            springs[j-1]->update_force();
            prv_rnds.push_back({});

        }
    }

    bd_prefactor = sqrt(temperature/(2*dt*damp));

    this->init_ubend();
}

filament::~filament()
{
    for (bead *b : beads) delete b;
    for (spring *s : springs) delete s;
}

void filament::add_bead(vector<double> a, double spring_length, double stretching_stiffness, double max_ext_ratio){
    if (a.size() != 4) throw runtime_error("Wrong number of arguments in bead.");
    beads.push_back(new bead({a[0], a[1], a[2], a[3]}));
    prv_rnds.push_back({});
    if (beads.size() > 1){
        int j = (int) beads.size() - 1;
        springs.push_back(new spring(spring_length, stretching_stiffness, max_ext_ratio, this, {j-1,  j}));
        springs[j-1]->step();
    }
    if (damp == infty)
        damp = beads[0]->get_friction();
}

void filament::update_positions()
{
    ke_vel = 0.0;
    ke_vir = 0.0;
    size_t sa = beads.size();
    for (size_t i = 0; i < sa; i++) {
        vec_type new_rnds = vec_randn();
        vec_type f = beads[i]->get_force() + bd_prefactor*damp*(new_rnds + prv_rnds[i]);
        vec_type v = f / damp;
        vec_type pos = beads[i]->get_pos();
        prv_rnds[i] = new_rnds;
        ke_vel += abs2(v);
        ke_vir += -0.5 * dot(f, pos);
        vec_type new_pos = bc->pos_bc(pos + v*dt);
        beads[i]->set_pos(new_pos);
        beads[i]->reset_force();
    }
    for (spring *s : springs) s->step();
}

void filament::update_stretching()
{
    for (spring *s : springs) {
        s->update_force();
        s->filament_update();
    }
}

spring *filament::get_spring(int i)
{
    return springs[i];
}

void filament::update_d_strain(double g)
{
    for (bead *b : beads) {
        vec_type pos = b->get_pos();
        b->set_pos({pos.x + g * pos.y / bc->get_ybox(), pos.y});
    }
}

box *filament::get_box()
{
    return bc;
}

vec_type filament::get_force(int i)
{
    return beads[i]->get_force();
}

void filament::update_forces(int index, vec_type f)
{
    beads[index]->update_force(f);
}

void filament::pull_on_ends(double f)
{
    if (beads.size() < 2) return;
    int last = beads.size() - 1;
    vec_type dr = bc->rij_bc(beads[last]->get_pos() - beads[0]->get_pos());
    double len = abs(dr);

    beads[ 0  ]->update_force(-0.5*f*dr/len);
    beads[last]->update_force( 0.5*f*dr/len);
}

void filament::affine_pull(double f)
{
    if (beads.size() < 2) return;
    int last = beads.size() - 1;
    vec_type dr = bc->rij_bc(beads[last]->get_pos() - beads[0]->get_pos());
    double len = abs(dr);
    vec_type fcs = f * dr / len;

    for (int i = 0; i <= last; i++){
        double frac = (double(i)/double(last)-0.5);
        beads[i]->update_force(frac * fcs);
    }
}

vector<vector<double>> filament::output_beads(int fil)
{
    vector<vector<double>> out;
    for (size_t i = 0; i < beads.size(); i++) {
        out.push_back(beads[i]->output());
        out[i].push_back(double(fil));
    }
    return out;
}

vector<vector<double>> filament::output_springs(int fil)
{
    vector<vector<double>> out;
    for (size_t i = 0; i < springs.size(); i++) {
        out.push_back(springs[i]->output());
        out[i].push_back(double(fil));
    }
    return out;
}

vector<double> filament::output_thermo(int fil)
{
    return {get_kinetic_energy_vel(), get_kinetic_energy_vir(), get_potential_energy(), get_total_energy(), double(fil)};
}

string filament::write_beads(int fil)
{
    string all_beads;
    for (bead *b : beads) {
        all_beads += fmt::format("{}\t{}", b->write(), fil);
    }
    return all_beads;
}

string filament::write_springs(int fil)
{
    string all_springs;
    for (spring *s : springs) {
        all_springs += fmt::format("{}\t{}", s->write(), fil);
    }
    return all_springs;
}

string filament::write_thermo(int fil)
{
    return fmt::format(
            "\n{}\t{}\t{}\t{}",
            get_kinetic_energy_vel(), get_kinetic_energy_vir(),
            get_potential_energy(), get_total_energy());
}

vector<vector<double>> filament::get_beads(size_t first, size_t last)
{
    vector<vector<double>> newbeads;
    for (size_t i = first; i < last; i++) {
        if (i >= beads.size()) {
            break;
        } else {
            bead *b = beads[i];
            vec_type pos = b->get_pos();
            newbeads.push_back({pos.x, pos.y, b->get_length(), b->get_viscosity()});
        }
    }
    return newbeads;
}

vector<filament *> filament::try_fracture()
{
    for (size_t i = 0; i < springs.size(); i++) {
        vec_type f = springs[i]->get_force();
        if (abs2(f) > fracture_force_sq) {
            return fracture(i);
        }
    }
    return {};
}

vector<filament *> filament::fracture(int node){

    vector<filament *> newfilaments;
    cout<<"\n\tDEBUG: fracturing at node "<<node;

    if(springs.size() == 0)
        return newfilaments;

    vector<vector<double>> lower_half = this->get_beads(0, node+1);
    vector<vector<double>> upper_half = this->get_beads(node+1, beads.size());

    if (lower_half.size() > 0)
        newfilaments.push_back(
                new filament(filament_network, lower_half, springs[0]->get_l0(), springs[0]->get_kl(), springs[0]->get_fene_ext(), kb,
                    dt, temperature, fracture_force));
    if (upper_half.size() > 0)
        newfilaments.push_back(
                new filament(filament_network, upper_half, springs[0]->get_l0(), springs[0]->get_kl(), springs[0]->get_fene_ext(), kb,
                    dt, temperature, fracture_force));

    return newfilaments;

}

void filament::detach_all_motors()
{
    for (spring *s : springs) {
        for (auto it : s->get_mots()) {
            motor *m = it.first;
            int hd = it.second;
            m->detach_head_without_moving(hd);
        }
    }
}

bool filament::operator==(const filament& that){

    if (beads.size() != that.beads.size() || springs.size() != that.springs.size())
        return false;

    for (unsigned int i = 0; i < beads.size(); i++)
        if (!(*(beads[i]) == *(that.beads[i])))
            return false;

    for (unsigned int i = 0; i < springs.size(); i++)
        if (!(springs[i]->is_similar(*(that.springs[i]))))
            return false;

    return (this->temperature == that.temperature &&
            this->dt == that.dt && this->fracture_force == that.fracture_force);

}

string filament::to_string()
{
    // Note: not including springs in to_string, because spring's to_string includes filament's to_string
    string out = "\n";

    for (bead *b : beads) {
        out += b->to_string();
    }

    out += fmt::format(
            "temperature = {}\t"
            "dt = {}\t"
            "fracture_force = {}\n",
            temperature, dt, fracture_force);

    return out;
}

inline double filament::angle_between_springs(int i, int j){

    // 1st bond
    vec_type delr1 = springs[i]->get_disp();
    double r1 = springs[i]->get_length();

    // 2nd bond
    vec_type delr2 = springs[j]->get_disp();
    double r2 = springs[j]->get_length();

    // cos angle
    double c = dot(delr1, delr2) / (r1 * r2);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    return acos(c);
}


void filament::update_bending()
{
    if (springs.size() <= 1 || kb == 0) return;

    bending_virial.zero();
    ubend = 0.0;
    for (int n = 0; n < int(springs.size())-1; n++) {

        vec_type delr1 = springs[n]->get_disp();
        vec_type delr2 = springs[n+1]->get_disp();

        bend_result_type result = bend_harmonic(kb, 0.0, delr1, delr2);

        ubend += result.energy;

        // apply force to each of 3 atoms
        beads[n  ]->update_force(-result.force1);
        beads[n+1]->update_force(result.force1);

        beads[n+1]->update_force(-result.force2);
        beads[n+2]->update_force(result.force2);

        bending_virial += outer(delr1, result.force1);
        bending_virial += outer(delr2, result.force2);
    }
}


int filament::get_nbeads(){
    return beads.size();
}

int filament::get_nsprings(){
    return springs.size();
}

double filament::get_bending_energy(){

    return ubend;

}

virial_type filament::get_bending_virial()
{
    return bending_virial;
}

void filament::init_ubend()
{
    ubend = 0.0;
    for (int i = 0; i < int(springs.size()) - 1; i++) {
        double theta = angle_between_springs(i + 1, i);
        ubend += 0.5 * kb * theta * theta;
    }
}

double filament::get_stretching_energy()
{
    double u = 0.0;
    for (spring *s : springs) {
        u += s ->get_stretching_energy();
    }
    return u;
}

double filament::get_kinetic_energy_vel()
{
    return ke_vel;
}

double filament::get_kinetic_energy_vir()
{
    return ke_vir;
}

virial_type filament::get_stretching_virial()
{
    virial_type vir;
    vir.zero();
    for (spring *s : springs) {
        vir += s->get_virial();
    }
    return vir;
}

double filament::get_potential_energy()
{
    return this->get_stretching_energy() + this->get_bending_energy();
}

double filament::get_total_energy()
{
    return this->get_potential_energy() + (this->get_kinetic_energy_vel());
}

vec_type filament::get_bead_position(int n)
{
    return beads[n]->get_pos();
}

void filament::print_thermo()
{
    fmt::print("\tKEvel = {}\tKEvir = {}\tPE = {}\tTE = {}",
            get_kinetic_energy_vel(), get_kinetic_energy_vir(),
            get_potential_energy(), get_total_energy());
}

double filament::get_end2end()
{
    if (beads.size() < 2) {
        return 0;
    } else {
        return bc->dist_bc(beads[beads.size() - 1]->get_pos() - beads[0]->get_pos());
    }
}

// GROWING

void filament::set_l0_max(double lmax)
{
    l0_max = lmax;
}

void filament::set_nsprings_max(int nmax)
{
    nsprings_max = nmax;
}

void filament::set_l0_min(double lmin)
{
    l0_min = lmin;
}

void filament::grow(double dL)
{
    double lb = springs[0]->get_l0();
    if ( lb + dL < l0_max ){
        springs[0]->set_l0(lb + dL);
        springs[0]->step();
    }
    else{
        vec_type dir = springs[0]->get_direction();
        vec_type p2 = beads[1]->get_pos();
        //add a bead
        vec_type newpos = bc->pos_bc(p2-spring_l0*dir);
        beads.insert(beads.begin()+1, new bead(newpos.x, newpos.y, beads[0]->get_length(), beads[0]->get_viscosity()));
        prv_rnds.insert(prv_rnds.begin()+1, {});
        //shift all springs from "1" onward forward
        //move backward; otherwise i'll just keep pushing all the motors to the pointed end, i think
        for (int i = int(springs.size()-1); i > 0; i--){
            springs[i]->inc_aindex();
            springs[i]->step();
            //shift all xsprings on these springs forward
            //lmots = springs[i]->get_mots();
            for (auto &it : springs[i]->get_mots()) {
                it.first->inc_l_index(it.second);
            }

        }
        //add spring "1"
        springs.insert(springs.begin()+1, new spring(spring_l0, springs[0]->get_kl(), springs[0]->get_max_ext(), this, {{1, 2}}));
        springs[1]->step();
        //reset l0 at barbed end
        springs[0]->set_l0(spring_l0);
        springs[0]->step();

        //adjust motors and xsprings on first spring
        map<motor *, int> mots0 = springs[0]->get_mots();
        vector<motor *> mots1;
        for (auto &it : mots0) {
            motor *m = it.first;
            int hd = it.second;
            double pos = m->get_pos_a_end()[hd];
            if (pos < spring_l0) {
                mots1.push_back(m);
            } else {
                m->set_pos_a_end(hd, pos - spring_l0);
            }
        }
        for (motor *m : mots1) {
            int hd = mots0[m];
            m->set_l_index(hd, 1);
        }
    }
}

void filament::update_length()
{
    if ( kgrow*lgrow > 0 && this->get_nsprings() + 1 <= nsprings_max && rng_u() < kgrow*dt){
        grow(lgrow);
    }
}

void filament::set_kgrow(double k){
    kgrow = k;
}

void filament::set_lgrow(double dl){
    lgrow = dl;
}
