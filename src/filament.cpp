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
        double stretching_stiffness, double bending_stiffness,
        double deltat, double temp, double frac_force)
{
    filament_network = net;

    bc = net->get_box();
    dt = deltat;
    temperature = temp;
    fracture_force = frac_force;
    fracture_force_sq = fracture_force*fracture_force;
    kb = bending_stiffness;

    ubend = 0.0;

    spring_l0 = spring_length;

    nsprings_max = 0;
    l0_max = 0.0;
    l0_min = 0.0;
    kgrow = 0.0;
    lgrow = 0.0;

    damp = infty;
    if (beadvec.size() > 0)
    {
        vector<double> &entry = beadvec[0];
        if (entry.size() != 4) throw runtime_error("Wrong number of arguments in beadvec.");
        beads.push_back(new bead(entry[0], entry[1], entry[2], entry[3]));
        prv_rnds.push_back(vec_randn());
        damp = beads[0]->get_friction();
    }

    //spring em up
    if (beadvec.size() > 1){
        for (unsigned int j = 1; j < beadvec.size(); j++) {

            vector<double> &entry = beadvec[j];
            if (entry.size() != 4) throw runtime_error("Wrong number of arguments in beadvec.");
            beads.push_back(new bead(entry[0], entry[1], entry[2], entry[3]));
            springs.push_back(new spring(spring_length, stretching_stiffness, this, {(int)j-1, (int)j}));
            springs[j-1]->step();
            springs[j-1]->update_force();
            prv_rnds.push_back(vec_randn());
        }
    }

    bd_prefactor = sqrt(temperature/(2*dt*damp));
}

filament::~filament()
{
    for (bead *b : beads) delete b;
    for (spring *s : springs) delete s;
}

void filament::add_bead(vector<double> a, double spring_length, double stretching_stiffness)
{
    if (a.size() != 4) throw runtime_error("Wrong number of arguments in bead.");
    beads.push_back(new bead({a[0], a[1], a[2], a[3]}));
    prv_rnds.push_back(vec_randn());
    if (beads.size() > 1){
        int j = (int) beads.size() - 1;
        springs.push_back(new spring(spring_length, stretching_stiffness, this, {j-1,  j}));
        springs[j-1]->step();
    }
    if (damp == infty) {
        damp = beads[0]->get_friction();
        bd_prefactor = sqrt(temperature/(2*dt*damp));
    }
}

void filament::update_positions()
{
    size_t sa = beads.size();
    for (size_t i = 0; i < sa; i++) {
        vec_type new_rnds = vec_randn();
        vec_type v = beads[i]->get_force() / damp + bd_prefactor * (new_rnds + prv_rnds[i]);
        vec_type pos = beads[i]->get_pos();
        prv_rnds[i] = new_rnds;
        vec_type new_pos = bc->pos_bc(pos + v * dt);
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

// TODO: energy, virial
void filament::pull_on_ends(double f)
{
    if (beads.size() < 2) return;
    int last = beads.size() - 1;
    vec_type dr = bc->rij_bc(beads[last]->get_pos() - beads[0]->get_pos());
    double len = abs(dr);

    beads[ 0  ]->update_force(-0.5*f*dr/len);
    beads[last]->update_force( 0.5*f*dr/len);
}

// TODO: energy, virial
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
    return {
        this->get_stretching_energy(),
        this->get_bending_energy(),
        double(fil)
    };
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
            "\n{}\t{}",
            this->get_stretching_energy(), this->get_bending_energy());
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
        springs[i]->update_force();
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
                new filament(filament_network, lower_half,
                    springs[0]->get_l0(), springs[0]->get_kl(), kb,
                    dt, temperature, fracture_force));
    if (upper_half.size() > 0)
        newfilaments.push_back(
                new filament(filament_network, upper_half,
                    springs[0]->get_l0(), springs[0]->get_kl(), kb,
                    dt, temperature, fracture_force));

    return newfilaments;

}

void filament::detach_all_motors()
{
    for (size_t i = 0; i < attached.size(); i++) {
        motor *m = attached[i].m;
        int hd = attached[i].hd;
        m->detach_head_without_moving(hd);
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

void filament::update_bending()
{
    if (springs.size() <= 1 || kb == 0) return;

    bending_virial.zero();
    ubend = 0.0;
    for (int n = 0; n < int(springs.size())-1; n++) {

        vec_type delr1 = springs[n+0]->get_disp();
        vec_type delr2 = springs[n+1]->get_disp();

        bend_result_type result = bend_harmonic(kb, 0.0, delr1, delr2);

        ubend += result.energy;

        // apply force to each of 3 atoms
        beads[n+0]->update_force(-result.force1);
        beads[n+1]->update_force(result.force1);

        beads[n+1]->update_force(-result.force2);
        beads[n+2]->update_force(result.force2);

        bending_virial += -0.5 * outer(delr1, result.force1);
        bending_virial += -0.5 * outer(delr2, result.force2);
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

double filament::get_stretching_energy()
{
    double u = 0.0;
    for (spring *s : springs) {
        u += s ->get_stretching_energy();
    }
    return u;
}

virial_type filament::get_stretching_virial()
{
    virial_type vir;
    for (spring *s : springs) {
        vir += s->get_virial();
    }
    return vir;
}

vec_type filament::get_bead_position(int n)
{
    return beads[n]->get_pos();
}

void filament::print_thermo()
{
    fmt::print("\tPEs = {}\tPEb = {}",
            this->get_stretching_energy(),
            this->get_bending_energy());
}

double filament::get_end2end()
{
    if (beads.size() < 2) {
        return 0;
    } else {
        return bc->dist_bc(beads[beads.size() - 1]->get_pos() - beads[0]->get_pos());
    }
}

// begin [growing]

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

    if (lb + dL < l0_max) {
        // make spring "0" longer
        springs[0]->set_l0(lb + dL);
        springs[0]->step();

    } else {

        vec_type dir = springs[0]->get_direction();
        vec_type p2 = beads[1]->get_pos();

        // split spring "0" into two

        // add a new bead "1" to split spring "0"
        vec_type newpos = bc->pos_bc(p2 - spring_l0 * dir);
        bead *b = new bead(
                newpos.x, newpos.y,
                beads[0]->get_length(), beads[0]->get_viscosity());
        beads.insert(beads.begin() + 1, b);
        prv_rnds.insert(prv_rnds.begin() + 1, vec_randn());

        // shift all springs forward, except the first one
        for (size_t i = 1; i < springs.size(); i++) {
            springs[i]->inc_aindex();
        }

        // add new spring "1" with length l0
        spring *s = new spring(spring_l0, springs[0]->get_kl(), this, {1, 2});
        springs.insert(springs.begin() + 1, s);

        // set spring at barbed end to remaining length
        springs[0]->set_l0(lb + dL - spring_l0);

        // update springs (I think only springs "0" and "1" need updating)
        for (spring *s : springs) {
            s->step();
        }

        // fix attached points
        for (auto &a : attached) {
            if (a.l > 0) {
                // shift attachment points forward if they aren't on the first spring
                a.l += 1;
            } else {
                // adjust attachment points on the first spring
                double pos = a.pos;
                if (pos < spring_l0) {
                    a.l = 1;
                } else {
                    a.pos = pos - spring_l0;
                }
            }
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

// end [growing]

// begin [attached]

// TODO: use an in-place linked list instead
// pos = distance from pointed end, relative to bead[l + 1]

int filament::new_attached(motor *m, int hd, int l, vec_type intpoint)
{
    double pos = bc->dist_bc(springs[l]->get_h1() - intpoint);
    // length of attached is usually short,
    // so this isn't very expensive
    for (size_t i = 0; i < attached.size(); i++) {
        if (!attached[i].m) {
            attached[i] = {m, hd, l, pos};
            return int(i);
        }
    }
    size_t j = attached.size();
    attached.push_back({m, hd, l, pos});
    return int(j);
}

void filament::del_attached(int i)
{
    attached[i] = {nullptr, -1, -1, NAN};
}

int filament::get_attached_l(int i)
{
    return attached[i].l;
}

vec_type filament::get_attached_pos(int i)
{
    int l = attached[i].l;
    double pos = attached[i].pos;
    vec_type h1 = springs[l]->get_h1();
    vec_type dir = springs[l]->get_direction();
    return bc->pos_bc(h1 - pos * dir);
}

void filament::add_attached_force(int i, vec_type f)
{
    int l = attached[i].l;
    double pos = attached[i].pos;
    double ratio = pos / springs[l]->get_length();
    beads[l + 0]->update_force(f * ratio);
    beads[l + 1]->update_force(f * (1.0 - ratio));
}

void filament::add_attached_pos(int i, double dist)
{
    int &l = attached[i].l;
    double &pos = attached[i].pos;
    pos += dist;

    double len = springs[l]->get_l0();
    if (pos >= len) {
        // pos is after spring
        if (l == 0) {
            // at barbed end
            pos = len;
        } else {
            // move to next spring
            l -= 1;
            // subtract CURRENT spring length
            pos -= len;
        }
    } else if (pos < 0.0) {
        // pos is before spring
        if (l + 1 == springs.size()) {
            // at pointed end
            pos = 0.0;
        } else {
            // move to previous spring
            l += 1;
            // add NEW spring length
            pos += springs[l]->get_l0();
        }
    }
}

bool filament::at_barbed_end(int i)
{
    return attached[i].l == 0 && attached[i].pos == springs[attached[i].l]->get_l0();
}

bool filament::at_pointed_end(int i)
{
    return attached[i].l + 1 == springs.size() && attached[i].pos == 0.0;
}

// end [attached]
