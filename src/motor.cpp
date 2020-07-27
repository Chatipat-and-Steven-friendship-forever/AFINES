/*------------------------------------------------------------------
 motor.cpp : object describing a motor or crosslinker

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
#include "filament_ensemble.h"
#include "potentials.h"

motor::motor(vector<double> mvec,
        double mlen, filament_ensemble * network,
        double delta_t,
        double temp,
        double v0,
        double stiffness,
        double max_ext_ratio,
        double ron, double roff, double rend,
        double fstall, double rcut,
        double vis)
{
    bc = network->get_box();
    vs          = v0;
    mk          = stiffness;

    stall_force = fstall;
    temperature = temp;

    max_bind_dist    = rcut;
    max_bind_dist_sq = rcut*rcut;

    mld         = mlen;
    dt          = delta_t;
    kon         = ron*dt;
    koff        = roff*dt;
    kend        = rend*dt;
    kon2        = ron*dt;
    koff2       = roff*dt;
    kend2       = rend*dt;

    kb = 0;
    th0 = 0;

    kalign = 0;
    par_flag = true;

    // filament and spring indices for each head
    array<int, 2> f_index = {int(mvec[4]), int(mvec[5])};
    array<int, 2> l_index = {int(mvec[6]), int(mvec[7])};

    state[0] = (f_index[0] == -1 && l_index[0] == -1) ? motor_state::free : motor_state::bound;
    state[1] = (f_index[1] == -1 && l_index[1] == -1) ? motor_state::free : motor_state::bound;

    filament_network = network;
    damp        =(6*pi*vis*mld);
    bd_prefactor= sqrt(temperature/(2*damp*dt));

    /********for FENE springs*********/
    max_ext     = max_ext_ratio*mlen;
    eps_ext     = 0.01*max_ext;
    /********************************/

    tension = 0;
    ke_vel = 0; //assume m = 1
    ke_vir = 0;
    b_eng = {0, 0};

    h[0] = bc->pos_bc({mvec[0], mvec[1]});
    h[1] = bc->pos_bc({mvec[0] + mvec[2], mvec[1] + mvec[3]});
    //force can be non-zero and angle is determined from disp vector
    this->update_angle();
    this->update_force();

    if (state[0] == motor_state::bound){
        fp_index[0] = filament_network->new_attached(this, 0, f_index[0], l_index[0], h[0]);
        ldir_bind[0] = filament_network->get_attached_direction(fp_index[0]);

    }
    if (state[1] == motor_state::bound){
        fp_index[1] = filament_network->new_attached(this, 1, f_index[1], l_index[1], h[1]);
        ldir_bind[1] = filament_network->get_attached_direction(fp_index[1]);
    }
}

void motor::set_bending(double kb_, double th0_)
{
    kb = kb_;
    th0 = th0_;
}

double motor::get_kb()
{
    return kb;
}

double motor::get_th0()
{
    return th0;
}

void motor::set_par(double k)
{
    kalign = k;
    par_flag = true;
}

void motor::set_antipar(double k)
{
    kalign = k;
    par_flag = false;
}

//return motor state with a given head number
array<motor_state, 2> motor::get_states()
{
    return state;
}

vec_type motor::get_h0()
{
    return h[0];
}


vec_type motor::get_h1()
{
    return h[1];
}

// metropolis algorithm
// probability is in range [0, 1]
double motor::metropolis_prob(int hd, array<int, 2> fl_idx, vec_type newpos)
{
    double len_old = bc->dist_bc(h[hd] - h[pr(hd)]) - mld;
    double len_new = bc->dist_bc(newpos - h[pr(hd)]) - mld;
    double dE = 0.5 * mk * (len_new * len_new - len_old * len_old);

    if (kb > 0.0 && state[pr(hd)] == motor_state::bound) {
        array<int, 2> fl;
        vec_type delr1, delr2;
        if (state[hd] == motor_state::free) {
            // attach: singly bound -> doubly bound

            // new bending energy of attaching head
            fl = fl_idx;
            delr1 = filament_network->get_filament(fl[0])->get_spring(fl[1])->get_disp();
            delr2 = pow(-1, hd) * disp;
            dE += bend_harmonic_energy(kb, th0, delr1, delr2);

            // new bending energy of other head
            fl = filament_network->get_attached_fl(fp_index[pr(hd)]);
            delr1 = filament_network->get_filament(fl[0])->get_spring(fl[1])->get_disp();
            delr2 = pow(-1, pr(hd)) * disp;
            dE += bend_harmonic_energy(kb, th0, delr1, delr2);

        }
        if (state[hd] == motor_state::bound) {
            // detach: doubly bound -> singly bound

            // old bending energy of detaching head
            fl = filament_network->get_attached_fl(fp_index[hd]);
            delr1 = filament_network->get_filament(fl[0])->get_spring(fl[1])->get_disp();
            delr2 = pow(-1, hd) * disp;
            dE -= bend_harmonic_energy(kb, th0, delr1, delr2);

            // old bending energy of other head
            fl = filament_network->get_attached_fl(fp_index[pr(hd)]);
            delr1 = filament_network->get_filament(fl[0])->get_spring(fl[1])->get_disp();
            delr2 = pow(-1, pr(hd)) * disp;
            dE -= bend_harmonic_energy(kb, th0, delr1, delr2);

        }
    }

    if (kalign != 0.0 && state[pr(hd)] == motor_state::bound) {
        array<int, 2> fl;
        if (state[hd] == motor_state::free) {
            // attach: singly bound -> doubly bound

            // spring to be attached
            fl = fl_idx;
            vec_type delr1 = filament_network->get_filament(fl[0])->get_spring(fl[1])->get_disp();

            // other attached spring
            fl = filament_network->get_attached_fl(fp_index[pr(hd)]);
            vec_type delr2 = filament_network->get_filament(fl[0])->get_spring(fl[1])->get_disp();

            dE += alignment_penalty(delr1, delr2);
        }
        if (state[hd] == motor_state::bound) {
            // detach: doubly bound -> singly bound

            // spring to be detached
            fl = filament_network->get_attached_fl(fp_index[hd]);
            vec_type delr1 = filament_network->get_filament(fl[0])->get_spring(fl[1])->get_disp();

            // other attached spring
            fl = filament_network->get_attached_fl(fp_index[pr(hd)]);
            vec_type delr2 = filament_network->get_filament(fl[0])->get_spring(fl[1])->get_disp();

            dE -= alignment_penalty(delr1, delr2);
        }
    }

    return (dE <= 0.0) ? 1.0 : exp(-dE / temperature);
}

// compute the alignment penalty for doubly bound heads
// values are in range [0, kalign], where
// 0 is for parallel/antiparallel and kalign is for antiparallel/parallel
double motor::alignment_penalty(vec_type delr1, vec_type delr2)
{
    // cosine of angle between vectors
    // if parallel, 1; if antiparallel, -1
    double c = dot(delr1, delr2) / (abs(delr1) * abs(delr2));

    // penalize NOT being parallel/antiparallel by at most kalign
    if (par_flag) {
        return kalign * 0.5 * (1.0 - c);
    } else {
        return kalign * 0.5 * (1.0 + c);
    }
}

bool motor::allowed_bind(int hd, array<int, 2> fl_idx){
    array<int, 2> fl = filament_network->get_attached_fl(fp_index[pr(hd)]);
    if (kb > 0.0) return fl_idx[0] != fl[0];
    return fl[0] != fl_idx[0] || fl[1] != fl_idx[1];
}

// attaches head to spring {f, l} at intpoint
// saves information for detachment
// assumes that intpoint is on the spring
void motor::attach_head(int hd, vec_type intpoint, array<int, 2> fl)
{
    // update state
    state[hd] = motor_state::bound;
    fp_index[hd] = filament_network->new_attached(this, hd, fl[0], fl[1], intpoint);

    // record displacement of head and orientation of spring for future unbinding move
    ldir_bind[hd] = filament_network->get_attached_direction(fp_index[hd]);
    bind_disp[hd] = bc->rij_bc(intpoint - h[hd]);

    // update head position
    // TODO: check that intpoint and get_attached_pos coincide
    h[hd] = intpoint;
}

// attempt to attach unbound head to a filament
// does NOT check that the head in unbound
//check for attachment of unbound heads given head index (0 for head 1, and 1 for head 2)
bool motor::try_attach(int hd, mc_prob &p)
{
    double onrate = (state[pr(hd)] == motor_state::bound) ? kon2 : kon;
    vector<array<int, 2>> *attach_list = filament_network->get_attach_list(h[hd]);

    int count = attach_list->size();
    double needprob = onrate * count;
    boost::optional<double> opt_p = p(needprob);

    if (opt_p) {
        double mf_rand = *opt_p;
        int i = floor(mf_rand / onrate);
        double remprob = mf_rand - onrate * i;

        if (i < 0) throw std::logic_error("attach list index < 0");
        if (i >= count) throw std::logic_error("attach list index >= count");
        if (remprob < 0 || remprob > onrate) throw std::logic_error("invalid remaining probability");

        // compute and get attachment point
        array<int, 2> fl = attach_list->at(i);
        filament *f = filament_network->get_filament(fl[0]);
        spring *s = f->get_spring(fl[1]);
        vec_type intpoint = s->intpoint(h[hd]);

        // don't bind if binding site is further away than the cutoff
        vec_type dr = bc->rij_bc(intpoint - h[hd]);
        double dist_sq = abs2(dr);
        if (dist_sq > max_bind_dist_sq || !allowed_bind(hd, fl)) {
            return false;
        }

        double prob = onrate * metropolis_prob(hd, fl, intpoint);

        if (remprob < prob) {
            attach_head(hd, intpoint, fl);
            return true;
        }
    }

    return false;
}

void motor::update_force()
{
    if (kb > 0.0) {
        b_force[0].zero();
        b_force[1].zero();
        if (state[0] == motor_state::bound && state[1] == motor_state::bound) {
            this->update_bending(0);
            this->update_bending(1);
        }
    }

    this->update_angle();

    tension = mk*(len - mld);
    force = tension * direc;
}

void motor::update_bending(int hd)
{
    array<int, 2> fl = filament_network->get_attached_fl(fp_index[hd]);
    vec_type delr1 = filament_network->get_filament(fl[0])->get_spring(fl[1])->get_disp();
    vec_type delr2 = pow(-1, hd) * disp;

    bend_result_type result = bend_harmonic(kb, th0, delr1, delr2);

    filament_network->update_forces(fl[0], fl[1], -result.force1);
    filament_network->update_forces(fl[0], fl[1] + 1, result.force1);

    b_force[pr(hd)] += result.force2;
    b_force[hd] -= result.force2;

    b_eng[hd] = result.energy;
}

/* Taken from hsieh, jain, larson, jcp 2006; eqn (5)
 * Adapted by placing a cutoff, similar to how it's done in LAMMPS src/bond_fene.cpp*/
void motor::update_force_fraenkel_fene()
{
    double ext = abs(mld - len);

    double scaled_ext;
    if (max_ext - ext > eps_ext )
        scaled_ext = ext/max_ext;
    else
        scaled_ext = (max_ext - eps_ext)/max_ext;

    double mkp = mk/(1-scaled_ext*scaled_ext)*(len-mld);
    force = mkp * direc;
}

// Brownian dynamics for unbound particles
// does NOT check that particles are unbound
void motor::brownian_relax(int hd)
{
    vec_type new_rnd = vec_randn();
    vec_type f = pow(-1, hd) * force + b_force[hd];
    vec_type v = f / damp + bd_prefactor * (new_rnd + prv_rnd[hd]);
    ke_vel = abs2(v);
    ke_vir = -0.5 * pow(-1, hd) * dot(f, h[hd]);
    h[hd] = bc->pos_bc(h[hd] + v*dt);
    prv_rnd[hd] = new_rnd;
}

void motor::kill_head(int hd)
{
    state[hd] = motor_state::dead;
}

void motor::deactivate_head(int hd)
{
    state[hd] = motor_state::inactive;
}

// move head in the direction of the other head
// so that the heads are 'mld' apart
void motor::relax_head(int hd)
{
    h[hd] = bc->pos_bc(h[pr(hd)] - pow(-1, hd)*mld*direc);
}

// compute 'disp', 'len', and 'direc' based on positions
void motor::update_angle()
{
    disp = bc->rij_bc(h[1] - h[0]);
    len = abs(disp);
    direc.zero();
    if (len != 0) direc = disp / len;
}

// compute unbinding position
// head must be bound
vec_type motor::generate_off_pos(int hd)
{
    vec_type ldir = filament_network->get_attached_direction(fp_index[hd]);
    double c = dot(ldir, ldir_bind[hd]);
    double s = cross(ldir, ldir_bind[hd]);

    vec_type bind_disp_rot = {bind_disp[hd].x*c - bind_disp[hd].y*s, bind_disp[hd].x*s + bind_disp[hd].y*c};

    return bc->pos_bc(h[hd] - bind_disp_rot);
}

// attempts to detach hd with maximum rate offrate
// the detachment position is determined by generate_off_pos
// returns true if detachment succeeds, and false otherwise
// assumes that the head is bound
bool motor::try_detach(int hd, mc_prob &p)
{
    vec_type hpos_new = generate_off_pos(hd);
    double offrate = filament_network->at_barbed_end(fp_index[hd]) ? kend : koff;
    if (state[pr(hd)] == motor_state::bound)
        offrate = filament_network->at_barbed_end(fp_index[hd]) ? kend2 : koff2;

    boost::optional<double> opt_p = p(offrate);
    if (opt_p) {
        double prob = offrate * metropolis_prob(hd, {-1, -1}, hpos_new);
        if (*opt_p < prob) {
            detach_head(hd, hpos_new);
            return true;
        }
    }
    return false;
}

// stepping and detachment kinetics of a single bound head
void motor::walk(int hd)
{
    if (vs == 0.0) return;
    if (vs > 0.0 && filament_network->at_barbed_end(fp_index[hd])) return;
    if (vs < 0.0 && filament_network->at_pointed_end(fp_index[hd])) return;

    //calculate motor velocity
    double vm = vs;
    if (state[pr(hd)] != motor_state::free) {
        vec_type dir = filament_network->get_attached_direction(fp_index[hd]);
        double f = pow(-1, hd) * dot(force, dir);
        double factor = 1.0 - f / stall_force;
        if (factor < 0.0) factor = 0.0;
        if (factor > 2.0) factor = 2.0;
        vm = factor * vs;
    }

    // update relative position
    filament_network->add_attached_pos(fp_index[hd], vm * dt);
}

// updates head position from filament
// assumes head is bound
void motor::update_position_attached(int hd)
{
    h[hd] = filament_network->get_attached_pos(fp_index[hd]);
}

// add force to filament bound to head
// assumes that head is bound
// Using the lever rule to propagate force as outlined in Nedelec F 2002
void motor::filament_update_hd(int hd, vec_type f)
{
    filament_network->add_attached_force(fp_index[hd], f);
}

// add force to filaments bound to heads
// heads may be in any state
void motor::filament_update()
{
    if (state[0] == motor_state::bound) this->filament_update_hd(0, force + b_force[0]);
    if (state[1] == motor_state::bound) this->filament_update_hd(1, -force + b_force[1]);

    b_force[0].zero();
    b_force[1].zero();
}

// detach head, and move it to the unbinding position
void motor::detach_head(int hd)
{
    vec_type hpos_new = generate_off_pos(hd);
    this->detach_head(hd, hpos_new);
}

// detach head, and move it to the specified position
void motor::detach_head(int hd, vec_type newpos)
{
    detach_head_without_moving(hd);
    h[hd] = newpos;
}

// detach head, and leave it at the same position
void motor::detach_head_without_moving(int hd)
{
    state[hd] = motor_state::free;
    filament_network->del_attached(fp_index[hd]);
    fp_index[hd] = {-1, -1};
}

array<int, 2> motor::get_f_index()
{
    array<int, 2> fl0 = filament_network->get_attached_fl(fp_index[0]);
    array<int, 2> fl1 = filament_network->get_attached_fl(fp_index[1]);
    return {fl0[0], fl1[0]};
}

array<int, 2> motor::get_l_index()
{
    array<int, 2> fl0 = filament_network->get_attached_fl(fp_index[0]);
    array<int, 2> fl1 = filament_network->get_attached_fl(fp_index[1]);
    return {fl0[1], fl1[1]};
}

vec_type motor::get_force()
{
    return force;
}

array<vec_type, 2> motor::get_b_force()
{
    return b_force;
}

double motor::get_stretching_energy()
{
    return abs2(force) / (2.0 * mk);
}

array<double, 2> motor::get_bending_energy()
{
    return b_eng;
}

virial_type motor::get_virial() {
    double k = tension / len;
    return {k * disp.x * disp.x, k * disp.x * disp.y,
            k * disp.x * disp.y, k * disp.y * disp.y};
}

double motor::get_stretching_energy_fene()
{
    double ext = abs(mld - len);

    if (max_ext - ext > eps_ext )
        return -0.5*mk*max_ext*max_ext*log(1-(ext/max_ext)*(ext/max_ext));
    else
        return 0.25*mk*ext*ext*(max_ext/eps_ext);

}

double motor::get_kinetic_energy_vel()
{
    return ke_vel;
}

double motor::get_kinetic_energy_vir()
{
    return ke_vir;
}

string motor::to_string()
{
    array<int, 2> fl0 = filament_network->get_attached_fl(fp_index[0]);
    array<int, 2> fl1 = filament_network->get_attached_fl(fp_index[1]);
    return fmt::format("\n"
            "head 0 position = ({}, {})\t "
            "head 1 position = ({}, {})\n"

            "state = ({}, {})\t "
            "f_index = ({}, {})\t "
            "l_index = ({}, {})\n"

            "viscosity = {}\t "
            "max binding distance = {}\t "
            "stiffness = {}\t "
            "stall force = {}\t "
            "length = {}\n"

            "kon = {}\t "
            "koff = {}\t "
            "kend = {}\t "
            "dt = {}\t "
            "temp = {}\t "
            "damp = {}\n"

            "tension = ({}, {})\n",

        h[0].x, h[0].y, h[1].x, h[1].y,
        static_cast<int>(state[0]),  static_cast<int>(state[1]),
        fl0[0],  fl1[0], fl0[1], fl1[1],
        vs, max_bind_dist, mk, stall_force, mld,
        kon, koff, kend, dt, temperature, damp,
        force.x, force.y);
}

vector<double> motor::output()
{
    array<int, 2> fl0 = filament_network->get_attached_fl(fp_index[0]);
    array<int, 2> fl1 = filament_network->get_attached_fl(fp_index[1]);
    return {h[0].x, h[0].y, disp.x, disp.y,
        double(fl0[0]), double(fl1[0]),
        double(fl0[1]), double(fl1[1])};
}

string motor::write()
{
    array<int, 2> fl0 = filament_network->get_attached_fl(fp_index[0]);
    array<int, 2> fl1 = filament_network->get_attached_fl(fp_index[1]);
    return fmt::format("\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            h[0].x, h[0].y,
            disp.x, disp.y,
            fl0[0], fl1[0],
            fl0[1], fl1[1]);
}

void motor::revive_head(int hd)
{
    state[hd] = motor_state::free;
}

void motor::update_d_strain(double g)
{
    array<double, 2> fov = bc->get_fov();
    h[0] = bc->pos_bc({h[0].x + g * h[0].y / fov[1], h[0].y});
    h[1] = bc->pos_bc({h[1].x + g * h[1].y / fov[1], h[1].y});
}

void motor::set_binding_two(double ron2, double roff2, double rend2){
    kon2  = ron2*dt;
    koff2 = roff2*dt;
    kend2 = rend2*dt;
}
