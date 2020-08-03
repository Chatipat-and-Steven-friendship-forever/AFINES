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

// TODO:
// - decouple damp from mld
// - add rotational force to head so forces aren't immediately applied to filaments
// - force/energy calculation less insane
// - fix virial calculation
// - make dead, inactive, and static motor head work properly
// - add check for state/fp_index consistency

motor::motor(vector<double> mvec,
        double mlen, filament_ensemble * network,
        double delta_t,
        double temp,
        double v0,
        double stiffness,
        double ron, double roff, double rend,
        double fstall, double rcut,
        double vis)
{
    bc = network->get_box();
    filament_network = network;

    // [parameters]
    dt = delta_t;
    temperature = temp;
    damp = 6 * pi * vis * mlen;
    bd_prefactor = sqrt(temperature / (2 * damp * dt));

    // attach/detach
    kon = kon2 = ron * dt;
    koff = koff2 = roff * dt;
    kend = kend2 = rend * dt;
    max_bind_dist = rcut;
    max_bind_dist_sq = rcut * rcut;

    // walk
    vs = v0;
    stall_force = fstall;

    // stretch
    mk = stiffness;
    mld = mlen;

    // bend
    kb = 0.0;  // deactivated
    th0 = 0.0;

    // align
    kalign = 0.0;  // deactivated
    par_flag = 0;

    // external
    ext = nullptr;  // deactivated

    // [state]

    // positions
    h[0] = bc->pos_bc({mvec[0], mvec[1]});
    h[1] = bc->pos_bc({mvec[0] + mvec[2], mvec[1] + mvec[3]});

    // filament and spring indices for each head
    array<int, 2> f_index = {int(mvec[4]), int(mvec[5])};
    array<int, 2> l_index = {int(mvec[6]), int(mvec[7])};

    // for bound heads,
    // set up attachment location with filaments
    // set detachment location as don't move
    // assumes that head positions are on filaments

    // for unbound heads,
    // set in same state as after detachment

    if (f_index[0] == -1 && l_index[0] == -1) {
        state[0] = motor_state::free;
        fp_index[0] = {-1, -1};
    } else {
        this->attach_head(0, h[0], {f_index[0], l_index[0]});
    }
    if (f_index[1] == -1 && l_index[1] == -1) {
        state[1] = motor_state::free;
        fp_index[1] = {-1, -1};
    } else {
        this->attach_head(1, h[1], {f_index[1], l_index[1]});
    }

    // set to N(0, 1) to prevent cooling
    prv_rnd[0] = vec_randn();
    prv_rnd[1] = vec_randn();

    // [derived]
    this->step();

    // [forces] and [thermo]
    // clear everything since setup isn't done yet
    this->clear_forces();
}

// begin [settings]

void motor::set_binding_two(double ron2, double roff2, double rend2){
    kon2  = ron2*dt;
    koff2 = roff2*dt;
    kend2 = rend2*dt;
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
    par_flag = 1;
}

void motor::set_align(double k)
{
    kalign = k;
    par_flag = 0;
}

void motor::set_antipar(double k)
{
    kalign = k;
    par_flag = -1;
}

double motor::get_kalign()
{
    return kalign;
}

int motor::get_align()
{
    return par_flag;
}

void motor::set_external(external *ext_)
{
    ext = ext_;
}

// move head in the direction of the other head
// so that the heads are 'mld' apart
void motor::relax_head(int hd)
{
    h[hd] = bc->pos_bc(h[pr(hd)] - pow(-1, hd)*mld*direc);
    this->step();
}

void motor::kill_head(int hd)
{
    state[hd] = motor_state::dead;
}

void motor::revive_head(int hd)
{
    state[hd] = motor_state::free;
}

void motor::deactivate_head(int hd)
{
    state[hd] = motor_state::inactive;
}

// end [settings]

// begin [state]

// return motor state with a given head number
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

// get total force
array<vec_type, 2> motor::get_force()
{
    return force;
}

// get spring force
array<vec_type, 2> motor::get_s_force()
{
    return s_force;
}

// get bending forces on motor heads
// part of the bending forces is already applied to filaments
array<vec_type, 2> motor::get_b_force()
{
    return b_force;
}

// get external force
array<vec_type, 2> motor::get_ext_force()
{
    return ext_force;
}

// get projected force
array<double, 2> motor::get_force_proj()
{
    return f_proj;
}

// end [state]

// begin [forces]
// assume that derived state is computed

// clear forces/energies/virials
void motor::clear_forces()
{
    force[0].zero();
    force[1].zero();

    tension = 0.0;
    s_force[0].zero();
    s_force[1].zero();
    s_eng = 0.0;
    vir_stretch.zero();

    b_force[0].zero();
    b_force[1].zero();
    b_eng[0] = 0.0;
    b_eng[1] = 0.0;
    vir_bend.zero();

    align_eng = 0.0;
    vir_align.zero();

    ext_force[0].zero();
    ext_force[1].zero();
    ext_eng[0] = 0.0;
    ext_eng[1] = 0.0;
    vir_ext.zero();

    f_proj[0] = 0.0;
    f_proj[1] = 0.0;
}

// update all forces/energies/virials
// adds them to filaments if needed
// (call this ONCE)
void motor::update_force()
{
    // spring forces
    tension = mk * (len - mld);
    vec_type sf = -tension * direc;
    s_force[0] = -sf;
    s_force[1] = sf;
    s_eng = 0.5 * mk * (len - mld) * (len - mld);
    vir_stretch = -0.5 * outer(disp, sf);

    // bending forces
    // partially applied to filaments
    if (kb > 0.0) {

        // bending forces are added (not assigned), so clear them first
        // also, bending forces may not be activated
        b_eng[0] = b_eng[1] = 0.0;
        b_force[0].zero();
        b_force[1].zero();
        vir_bend.zero();

        // bending forces are activated only when both heads are bound
        if (state[0] == motor_state::bound && state[1] == motor_state::bound) {
            this->update_bending(0);
            this->update_bending(1);
        }
    }

    // alignment forces
    // all applied to filaments
    if (kalign != 0.0) {

        // in case alignment isn't activated
        align_eng = 0.0;
        vir_align.zero();

        // alignment is activated only when both heads are bound
        if (state[0] == motor_state::bound && state[1] == motor_state::bound) {
            update_alignment();
        }
    }

    // external forces
    if (ext) {
        this->update_external(0);
        this->update_external(1);
    }

    force[0] = s_force[0] + b_force[0] + ext_force[0];
    force[1] = s_force[1] + b_force[1] + ext_force[1];

    // update projected force for walking
    // computed from other forces, so should be called last
    if (vs != 0.0) {
        this->update_force_proj(0);
        this->update_force_proj(1);
    }

    this->filament_update();
}

// updates bending forces
// applies part of the force to filaments (the other part is in filament_update)
void motor::update_bending(int hd)
{
    vec_type delr1 = filament_network->get_attached_disp(fp_index[hd]);
    vec_type delr2 = pow(-1, hd) * disp;

    bend_result_type result = bend_harmonic(kb, th0, delr1, delr2);

    filament_network->add_attached_disp_force(fp_index[hd], result.force1);

    b_force[pr(hd)] += result.force2;
    b_force[hd] -= result.force2;

    b_eng[hd] += result.energy;

    vir_bend += -0.5 * outer(delr1, result.force1);
    vir_bend += -0.5 * outer(delr2, result.force2);
}

// update forces and energies for alignment
void motor::update_alignment()
{
    vec_type delr0 = filament_network->get_attached_disp(fp_index[0]);
    vec_type delr1 = filament_network->get_attached_disp(fp_index[1]);

    bend_result_type result = bend_angle(delr0, delr1);

    // cosine of angle between vectors
    // if parallel, 1; if antiparallel, -1
    double c = result.energy;

    // penalize NOT being parallel/antiparallel by at most kalign
    double a;
    if (par_flag == 1) {
        a = -kalign;
        align_eng = kalign * (1.0 - c);
    } else if (par_flag == -1) {
        a = kalign;
        align_eng = kalign * (1.0 + c);
    } else {
        if (c > 0.0) {
            a = -kalign;
            align_eng = kalign * (1.0 - c);
        } else {
            a = kalign;
            align_eng = kalign * (1.0 + c);
        }
    }

    vec_type f0 = a * result.force1;
    vec_type f1 = a * result.force2;

    filament_network->add_attached_disp_force(fp_index[0], f0);
    filament_network->add_attached_disp_force(fp_index[1], f1);

    vir_align += -0.5 * outer(delr0, f0);
    vir_align += -0.5 * outer(delr1, f1);
}

// update external forces and energies
void motor::update_external(int hd)
{
    ext_result_type result = ext->compute(h[hd]);
    ext_eng[hd] = result.energy;
    ext_force[hd] = result.force;
    vir_ext = result.virial;
}

// update projected force for walking
void motor::update_force_proj(int hd)
{
    f_proj[hd] = 0.0;  // zero if unbound
    if (state[hd] == motor_state::bound && state[pr(hd)] != motor_state::free) {
        vec_type dir = filament_network->get_attached_direction(fp_index[hd]);
        f_proj[hd] = dot(force[hd], dir);
    }
}

// add force to filaments bound to heads
// heads may be in any state
// Using the lever rule to propagate force as outlined in Nedelec F 2002
void motor::filament_update()
{
    if (state[0] == motor_state::bound)
        filament_network->add_attached_force(fp_index[0], force[0]);
    if (state[1] == motor_state::bound)
        filament_network->add_attached_force(fp_index[1], force[1]);
}

// end [forces]

// begin [dynamics]
// assume that forces are calculated
// does NOT assume that derived state is calculated
// motor::step should be called after all these are done

// apply shear
// called automatically by box
void motor::update_d_strain(double g)
{
    array<double, 2> fov = bc->get_fov();
    h[0] = bc->pos_bc({h[0].x + g * h[0].y / fov[1], h[0].y});
    h[1] = bc->pos_bc({h[1].x + g * h[1].y / fov[1], h[1].y});
}

// Brownian dynamics for unbound particles
// does NOT check that particles are unbound
void motor::brownian_relax(int hd)
{
    vec_type new_rnd = vec_randn();
    vec_type v = force[hd] / damp + bd_prefactor * (new_rnd + prv_rnd[hd]);
    h[hd] = bc->pos_bc(h[hd] + v*dt);
    prv_rnd[hd] = new_rnd;
}

// stepping kinetics of a single bound head
void motor::walk(int hd)
{
    if (vs == 0.0) return;
    if (vs > 0.0 && filament_network->at_barbed_end(fp_index[hd])) return;
    if (vs < 0.0 && filament_network->at_pointed_end(fp_index[hd])) return;

    //calculate motor velocity
    double vm = vs;
    if (state[pr(hd)] != motor_state::free) {
        double factor = 1.0 - f_proj[hd] / stall_force;
        if (factor < 0.0) factor = 0.0;
        if (factor > 2.0) factor = 2.0;
        vm = factor * vs;
    }

    // update relative position
    filament_network->add_attached_pos(fp_index[hd], vm * dt);
}

// update positions from filaments if needed,
// and compute 'disp', 'len', and 'direc' based on positions
void motor::step()
{
    if (state[0] == motor_state::bound) h[0] = filament_network->get_attached_pos(fp_index[0]);
    if (state[1] == motor_state::bound) h[1] = filament_network->get_attached_pos(fp_index[1]);
    disp = bc->rij_bc(h[1] - h[0]);
    len = abs(disp);
    direc.zero();
    if (len != 0) direc = disp / len;
}

// end [dynamics]

// begin [attach/detach]
// assumes that derived state is computed
// if positions are changed, derived states are recomputed

// metropolis algorithm
// compute attachment/detachment probability between current and proposed state
// probability is in range [0, 1]
double motor::metropolis_prob(int hd, array<int, 2> fl_idx, vec_type newpos)
{
    // stretching
    double len_old = bc->dist_bc(h[hd] - h[pr(hd)]) - mld;
    double len_new = bc->dist_bc(newpos - h[pr(hd)]) - mld;
    double dE = 0.5 * mk * (len_new * len_new - len_old * len_old);

    // bending
    if (kb > 0.0 && state[pr(hd)] == motor_state::bound) {
        vec_type delr1, delr2;
        if (state[hd] == motor_state::free) {
            // attach: singly bound -> doubly bound

            // new bending energy of attaching head
            delr1 = filament_network->get_disp(fl_idx[0], fl_idx[1]);
            delr2 = pow(-1, hd) * disp;
            dE += bend_harmonic_energy(kb, th0, delr1, delr2);

            // new bending energy of other head
            delr1 = filament_network->get_attached_disp(fp_index[pr(hd)]);
            delr2 = pow(-1, pr(hd)) * disp;
            dE += bend_harmonic_energy(kb, th0, delr1, delr2);

        }
        if (state[hd] == motor_state::bound) {
            // detach: doubly bound -> singly bound

            // old bending energy of detaching head
            delr1 = filament_network->get_attached_disp(fp_index[hd]);
            delr2 = pow(-1, hd) * disp;
            dE -= bend_harmonic_energy(kb, th0, delr1, delr2);

            // old bending energy of other head
            delr1 = filament_network->get_attached_disp(fp_index[pr(hd)]);
            delr2 = pow(-1, pr(hd)) * disp;
            dE -= bend_harmonic_energy(kb, th0, delr1, delr2);

        }
    }

    // alignment
    if (kalign != 0.0 && state[pr(hd)] == motor_state::bound) {
        if (state[hd] == motor_state::free) {
            // attach: singly bound -> doubly bound

            // spring to be attached
            vec_type delr1 = filament_network->get_disp(fl_idx[0], fl_idx[1]);

            // other attached spring
            vec_type delr2 = filament_network->get_attached_disp(fp_index[pr(hd)]);

            dE += alignment_penalty(delr1, delr2);
        }
        if (state[hd] == motor_state::bound) {
            // detach: doubly bound -> singly bound

            // spring to be detached
            vec_type delr1 = filament_network->get_attached_disp(fp_index[hd]);

            // other attached spring
            vec_type delr2 = filament_network->get_attached_disp(fp_index[pr(hd)]);

            dE -= alignment_penalty(delr1, delr2);
        }
    }

    // external
    if (ext) {
        dE += ext->compute(newpos).energy - ext->compute(h[hd]).energy;
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
    double c = angle(delr1, delr2);

    // penalize NOT being parallel/antiparallel by at most kalign
    if (par_flag == 1) {
        return kalign * (1.0 - c);
    } else if (par_flag == -1) {
        return kalign * (1.0 + c);
    } else {
        if (c > 0.0) {
            return kalign * (1.0 - c);
        } else {
            return kalign * (1.0 + c);
        }
    }
}

// ATTACHMENT

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
        vec_type intpoint = filament_network->intpoint(fl[0], fl[1], h[hd]);

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
    this->step();
}

// DETACHMENT

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

// detach head, and leave it at the same position
void motor::detach_head_without_moving(int hd)
{
    state[hd] = motor_state::free;
    filament_network->del_attached(fp_index[hd]);
    fp_index[hd] = {-1, -1};
    ldir_bind[hd].zero();
    bind_disp[hd].zero();
}

// detach head, and move it to the specified position
void motor::detach_head(int hd, vec_type newpos)
{
    detach_head_without_moving(hd);
    h[hd] = newpos;
    this->step();
}

// end [attach/detach]

// begin [thermo]
// assume that derived state and forces/etc are computed

// energy

double motor::get_stretching_energy()
{
    return s_eng;
}

array<double, 2> motor::get_bending_energy()
{
    return b_eng;
}

double motor::get_alignment_energy()
{
    return align_eng;
}

array<double, 2> motor::get_external_energy()
{
    return ext_eng;
}

// virial

virial_type motor::get_stretching_virial()
{
    return vir_stretch;
}

virial_type motor::get_bending_virial()
{
    return vir_bend;
}

virial_type motor::get_alignment_virial()
{
    return vir_align;
}

virial_type motor::get_external_virial()
{
    return vir_ext;
}

// end [thermo]

// begin [output]

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
    return fmt::format("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            h[0].x, h[0].y,
            disp.x, disp.y,
            fl0[0], fl1[0],
            fl0[1], fl1[1]);
}

// end [output]
