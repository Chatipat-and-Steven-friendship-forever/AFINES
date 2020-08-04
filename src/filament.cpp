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

#include "box.h"
#include "motor.h"
#include "potentials.h"
#include "spring.h"

filament::filament(box *bc_, vector<vec_type> beadvec, vector<vec_type> randvec,
        double bead_radius_, double visc_,
        double spring_length,
        double stretching_stiffness, double bending_stiffness,
        double deltat, double temp, double frac_force)
{
    bc = bc_;
    dt = deltat;
    temperature = temp;
    fracture_force = frac_force;
    kb = bending_stiffness;
    bead_radius = bead_radius_;
    visc = visc_;
    damp = 6 * pi * visc * bead_radius;
    bd_prefactor = sqrt(temperature/(2*dt*damp));

    ubend = 0.0;

    spring_l0 = spring_length;

    nsprings_max = 0;
    l0_max = 0.0;
    l0_min = 0.0;
    kgrow = 0.0;
    lgrow = 0.0;

    positions = beadvec;
    forces.resize(positions.size());

    if (randvec.size() == beadvec.size()) {
        prv_rnds = randvec;
    } else if (randvec.size() == 0) {
        // initialize random numbers
        for (size_t i = 0; i < positions.size(); i++) {
            prv_rnds.push_back(vec_randn());
        }
    } else {
        throw std::runtime_error("Positions and random numbers have different sizes.");
    }

    // spring em up
    for (size_t j = 0; j + 1 < beadvec.size(); j++) {
        springs.push_back(new spring(bc, spring_length, stretching_stiffness));
        springs[j]->step(positions[j], positions[j + 1]);
    }
}

filament::~filament()
{
    for (spring *s : springs) delete s;
}

void filament::add_bead(vec_type a, double spring_length, double stretching_stiffness)
{
    positions.push_back(a);
    forces.push_back({});
    prv_rnds.push_back(vec_randn());
    if (positions.size() > 1) {
        int j = int(positions.size()) - 1;
        springs.push_back(new spring(bc, spring_length, stretching_stiffness));
        springs[j - 1]->step(positions[j - 1], positions[j]);
    }
}

void filament::update_positions()
{
    size_t sa = positions.size();
    for (size_t i = 0; i < sa; i++) {
        vec_type new_rnds = vec_randn();
        vec_type v = forces[i] / damp + bd_prefactor * (new_rnds + prv_rnds[i]);
        prv_rnds[i] = new_rnds;
        positions[i] = bc->pos_bc(positions[i] + v * dt);
        forces[i].zero();
    }
    for (size_t i = 0; i < springs.size(); i++) {
        springs[i]->step(positions[i], positions[i + 1]);
    }
}

void filament::update_stretching()
{
    for (size_t i = 0; i < springs.size(); i++) {
        vec_type f = springs[i]->get_force();
        forces[i + 0] += f;
        forces[i + 1] -= f;
    }
}

spring *filament::get_spring(int i)
{
    return springs[i];
}

void filament::update_d_strain(double g)
{
    for (vec_type &pos : positions) {
        pos.x += g * pos.y / bc->get_ybox();
    }
}

box *filament::get_box()
{
    return bc;
}

vec_type filament::get_force(int i)
{
    return forces[i];
}

void filament::update_forces(int index, vec_type f)
{
    forces[index] += f;
}

// TODO: energy, virial
void filament::pull_on_ends(double f)
{
    if (positions.size() < 2) return;
    int l = positions.size() - 1;
    vec_type dr = bc->rij_bc(positions[l] - positions[0]);
    double len = abs(dr);

    forces[0] -= 0.5 * f * dr / len;
    forces[l] += 0.5 * f * dr / len;
}

// TODO: energy, virial
void filament::affine_pull(double f)
{
    if (positions.size() < 2) return;
    int last = positions.size() - 1;
    vec_type dr = bc->rij_bc(positions[last] - positions[0]);
    double len = abs(dr);
    vec_type fcs = f * dr / len;

    for (int i = 0; i <= last; i++){
        double frac = (double(i)/double(last)-0.5);
        forces[i] += frac * fcs;
    }
}

vector<vector<double>> filament::output_beads(int fil)
{
    vector<vector<double>> out;
    for (vec_type pos : positions) {
        out.push_back({pos.x, pos.y, bead_radius, double(fil)});
    }
    return out;
}

vector<vector<double>> filament::output_springs(int fil)
{
    vector<vector<double>> out;
    for (spring *s : springs) {
        vec_type h0 = s->get_h0();
        vec_type disp = s->get_disp();
        out.push_back({h0.x, h0.y, disp.x, disp.y, double(fil)});
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
    for (vec_type pos : positions) {
        all_beads += fmt::format("{}\t{}\t{}\t{}\n", pos.x, pos.y, bead_radius, fil);
    }
    return all_beads;
}

string filament::write_springs(int fil)
{
    string all_springs;
    for (spring *s : springs) {
        vec_type h0 = s->get_h0();
        vec_type disp = s->get_disp();
        all_springs += fmt::format("{}\t{}\t{}\t{}\t{}\n", h0.x, h0.y, disp.x, disp.y, fil);
    }
    return all_springs;
}

string filament::write_thermo(int fil)
{
    return fmt::format("{}\t{}\t{}\n",
            this->get_stretching_energy(), this->get_bending_energy(), fil);
}

vector<filament *> filament::try_fracture()
{
    if (!std::isfinite(fracture_force)) return {};
    for (size_t i = 0; i < springs.size(); i++) {
        double kl = springs[i]->get_kl();
        double l0 = springs[i]->get_l0();
        double len = springs[i]->get_length();
        if (kl * abs(len - l0) > fracture_force) {
            return fracture(i);
        }
    }
    return {};
}

// TODO: keep motors attached at the correct positions
vector<filament *> filament::fracture(int node)
{
    if (node < 0 || node >= int(springs.size()))
        throw std::runtime_error(fmt::format("Spring index {} outside of range [0, {}).", node, springs.size()));

    fmt::print("DEBUG: fracturing at node {}\n", node);

    vector<vec_type> lower_pos(positions.begin(), positions.begin() + node + 1);
    vector<vec_type> upper_pos(positions.begin() + node + 1, positions.end());

    vector<vec_type> lower_ran(prv_rnds.begin(), prv_rnds.begin() + node + 1);
    vector<vec_type> upper_ran(prv_rnds.begin() + node + 1, prv_rnds.end());

    assert(lower_pos.size() > 0);
    assert(upper_pos.size() > 0);
    assert(lower_ran.size() > 0);
    assert(upper_ran.size() > 0);

    // create filaments with the same parameters
    // forces are zeroed, but they should be already zeroed anyway

    filament *upper = new filament(bc,
            lower_pos, lower_ran,
            bead_radius, visc,
            springs[0]->get_l0(), springs[0]->get_kl(), kb,
            dt, temperature, fracture_force);

    filament *lower = new filament(bc,
            upper_pos, upper_ran,
            bead_radius, visc,
            springs[0]->get_l0(), springs[0]->get_kl(), kb,
            dt, temperature, fracture_force);

    // set grow parameters

    upper->set_nsprings_max(nsprings_max);
    upper->set_l0_max(l0_max);
    upper->set_l0_min(l0_min);
    upper->set_kgrow(kgrow);
    upper->set_lgrow(lgrow);

    lower->set_nsprings_max(nsprings_max);
    lower->set_l0_max(l0_max);
    lower->set_l0_min(l0_min);
    lower->set_kgrow(kgrow);
    lower->set_lgrow(lgrow);

    return {upper, lower};
}

void filament::detach_all_motors()
{
    for (size_t i = 0; i < attached.size(); i++) {
        motor *m = attached[i].m;
        int hd = attached[i].hd;
        m->detach_head_without_moving(hd);
    }
}

// TODO: update with all fields
bool filament::operator==(const filament& that){

    if (positions.size() != that.positions.size() || springs.size() != that.springs.size())
        return false;

    if (positions != that.positions) return false;
    if (forces != that.forces) return false;
    if (prv_rnds != that.prv_rnds) return false;

    for (size_t i = 0; i < springs.size(); i++)
        if (!(*springs[i] == *that.springs[i])) return false;

    return temperature == that.temperature
        && dt == that.dt
        && fracture_force == that.fracture_force;
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
        forces[n+0] -= result.force1;
        forces[n+1] += result.force1;

        forces[n+1] -= result.force2;
        forces[n+2] += result.force2;

        bending_virial += -0.5 * outer(delr1, result.force1);
        bending_virial += -0.5 * outer(delr2, result.force2);
    }
}


int filament::get_nbeads()
{
    return positions.size();
}

int filament::get_nsprings()
{
    return springs.size();
}

double filament::get_bending_energy()
{
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
    return positions[n];
}

double filament::get_end2end()
{
    if (positions.size() < 2) {
        return 0;
    } else {
        return bc->dist_bc(positions[positions.size() - 1] - positions[0]);
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
        springs[0]->step(positions[0], positions[1]);

    } else {

        vec_type dir = springs[0]->get_direction();
        vec_type p2 = positions[1];

        // split spring "0" into two

        // add a new bead "1" to split spring "0"
        vec_type newpos = bc->pos_bc(p2 - spring_l0 * dir);
        positions.insert(positions.begin() + 1, newpos);
        forces.insert(forces.begin() + 1, {});
        prv_rnds.insert(prv_rnds.begin() + 1, vec_randn());

        // add new spring "1" with length l0
        spring *s = new spring(bc, spring_l0, springs[0]->get_kl());
        springs.insert(springs.begin() + 1, s);  // (positions[1], positions[2])

        // set spring at barbed end to remaining length
        springs[0]->set_l0(lb + dL - spring_l0);

        // update springs (I think only springs "0" and "1" need updating)
        for (size_t i = 0; i < springs.size(); i++) {
            springs[i]->step(positions[i], positions[i + 1]);
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
    forces[l + 0] += f * ratio;
    forces[l + 1] += f * (1.0 - ratio);
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
        if (l + 1 == int(springs.size())) {
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
    return attached[i].l + 1 == int(springs.size()) && attached[i].pos == 0.0;
}

// end [attached]
