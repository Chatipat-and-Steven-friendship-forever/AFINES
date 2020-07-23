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


    // filament and spring indices for each head
    f_index = {int(mvec[4]), int(mvec[5])};
    array<int,2 > mylindex = {int(mvec[6]), int(mvec[7])};

    state[0] = (f_index[0] == -1 && mylindex[0] == -1) ? motor_state::free : motor_state::bound;
    state[1] = (f_index[1] == -1 && mylindex[1] == -1) ? motor_state::free : motor_state::bound;

    filament_network = network;
    init_l_index(0, mylindex[0]);
    init_l_index(1, mylindex[1]);
    damp        =(6*pi*vis*mld);
    bd_prefactor= sqrt(temperature/(2*damp*dt));

    /********for FENE springs*********/
    max_ext     = max_ext_ratio*mlen;
    eps_ext     = 0.01*max_ext;
    /********************************/

    tension     = 0;
    force.zero(); // force on the spring
    ke_vel = 0; //assume m = 1
    ke_vir = 0;
    pos_a_end = {{0, 0}}; // pos_a_end = distance from pointy end -- by default 0
                        // i.e., if l_index[hd] = j, then pos_a_end[hd] is the distance to the "j+1"th bead

    h[0] = bc->pos_bc({mvec[0], mvec[1]});
    h[1] = bc->pos_bc({mvec[0] + mvec[2], mvec[1] + mvec[3]});
    //force can be non-zero and angle is determined from disp vector
    this->update_angle();
    this->update_force();

    ldir_bind[0].zero();
    ldir_bind[1].zero();
    bind_disp[0].zero();
    bind_disp[1].zero();

    at_barbed_end = {{false, false}};

    if (state[0] == motor_state::bound){
        pos_a_end[0] = bc->dist_bc(filament_network->get_end(f_index[0], l_index[0]) - h[0]);
        ldir_bind[0] = filament_network->get_direction(f_index[0], l_index[0]);

    }
    if (state[1] == motor_state::bound){
        pos_a_end[1] = bc->dist_bc(filament_network->get_end(f_index[1], l_index[1]) - h[1]);
        ldir_bind[1] = filament_network->get_direction(f_index[1], l_index[1]);
    }

    prv_rnd[0].zero();
    prv_rnd[1].zero();
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

//metropolis algorithm with rate constant
//NOTE: while fl_idx doesn't matter for this xlink implementation, it does for "spacers"
double motor::metropolis_prob(int hd, array<int, 2> fl_idx, vec_type newpos, double maxprob)
{
    double prob = maxprob;
    double stretch  = bc->dist_bc(newpos - h[pr(hd)]) - mld;
    double delE = 0.5*mk*stretch*stretch - this->get_stretching_energy();

    if( delE > 0 )
        prob *= exp(-delE/temperature);

    return prob;
}

bool motor::allowed_bind(int hd, array<int, 2> fl_idx){
    return (f_index[pr(hd)] != fl_idx[0] || l_index[pr(hd)] != fl_idx[1]);
}

void motor::attach_head(int hd, vec_type intpoint, array<int, 2> fl)
{
    // update state
    state[hd] = motor_state::bound;
    f_index[hd] = fl[0];
    set_l_index(hd, fl[1]);

    // record displacement of head and orientation of spring for future unbinding move
    ldir_bind[hd] = filament_network->get_direction(f_index[hd], l_index[hd]);
    bind_disp[hd] = bc->rij_bc(intpoint - h[hd]);

    // update head position
    h[hd] = intpoint;

    // update relative head position
    vec_type ref = filament_network->get_end(f_index[hd], l_index[hd]);
    pos_a_end[hd] = bc->dist_bc(ref - intpoint);

    // even if the head is at the barbed end upon binding,
    // it could have negative velocity,
    // so always set this to false until it steps
    at_barbed_end[hd] = false;
}

// does the same thing as motor::attach, but faster
bool motor::attach_opt(int hd)
{
    double mf_rand = rng_u();
    double onrate = (state[pr(hd)] == motor_state::bound) ? kon2 : kon;
    vector<array<int, 2>> *attach_list = filament_network->get_attach_list(h[hd]);
    int count = attach_list->size();
    if (onrate * count > 1.0) throw std::runtime_error("onrate * count > 1");
    int i = floor(mf_rand / onrate);
    if (i < 0) throw std::logic_error("attach list index < 0");
    if (i < count) {
        double remprob = mf_rand - onrate * i;
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

        double prob = metropolis_prob(hd, fl, intpoint, onrate);

        if (remprob < prob) {
            attach_head(hd, intpoint, fl);
            return true;
        }
    }
    return false;
}

// attempt to attach unbound head to a filament
// does NOT check that the head in unbound
//check for attachment of unbound heads given head index (0 for head 1, and 1 for head 2)
bool motor::attach(int hd)
{
    set<tuple<double, array<int, 2>, vec_type>> dist_sq_sorted;
    for (array<int, 2> fl : *filament_network->get_attach_list(h[hd])) {
        filament *f = filament_network->get_filament(fl[0]);
        spring *s = f->get_spring(fl[1]);
        vec_type intpoint = s->intpoint(h[hd]);
        vec_type dr = bc->rij_bc(intpoint - h[hd]);
        double dist_sq = abs2(dr);
        dist_sq_sorted.insert(tuple<double, array<int, 2>, vec_type>(dist_sq, fl, intpoint));
    }

    double not_off_prob = 0.0;
    double mf_rand = rng_u();
    double onrate = (state[pr(hd)] == motor_state::bound) ? kon2 : kon;
    for (auto it : dist_sq_sorted) {
        double dist_sq; array<int, 2> fl; vec_type intpoint;
        std::tie(dist_sq, fl, intpoint) = it;
        if (dist_sq > max_bind_dist_sq) //since it's sorted, all the others will be farther than max_bind_dist too
            break;

        //head can't bind to the same filament spring the other head is bound to
        else if (allowed_bind(hd, fl)) {
            not_off_prob += metropolis_prob(hd, fl, intpoint, onrate);

            if (mf_rand < not_off_prob) {
                attach_head(hd, intpoint, fl);
                return true;
            }
        }
    }
    return false;
}

void motor::update_force()
{
    tension = mk*(len - mld);
    force = tension * direc;
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
    vec_type v = pow(-1, hd) * force / damp + bd_prefactor * (new_rnd + prv_rnd[hd]);
    ke_vel = abs2(v);
    ke_vir = -0.5 * pow(-1, hd) * dot(force, h[hd]);
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
    if (len != 0) direc = disp / len;
    else direc.zero();
}

// compute unbinding position
// head must be bound
vec_type motor::generate_off_pos(int hd)
{
    vec_type ldir = filament_network->get_direction(f_index[hd], l_index[hd]);
    double c = dot(ldir, ldir_bind[hd]);
    double s = cross(ldir, ldir_bind[hd]);

    vec_type bind_disp_rot = {bind_disp[hd].x*c - bind_disp[hd].y*s, bind_disp[hd].x*s + bind_disp[hd].y*c};

    return bc->pos_bc(h[hd] - bind_disp_rot);
}

// stepping and detachment kinetics of a single bound head
void motor::step_onehead(int hd)
{

    vec_type hpos_new = generate_off_pos(hd);
    double offrate = at_barbed_end[hd] ? kend : koff;
    if (state[pr(hd)] == motor_state::bound)
        offrate = at_barbed_end[hd] ? kend2 : koff2;

    double off_prob = metropolis_prob(hd, {{0,0}}, hpos_new, offrate);

    //cout<<"\nDEBUG: at barbed end? : "<<at_barbed_end[hd]<<"; off_prob = "<<off_prob;
    // attempt detachment
    if ( event(off_prob) ) this->detach_head(hd, hpos_new);
    else{

        //calculate motor velocity
        if (vs != 0 && !(at_barbed_end[hd])){
            double vm = vs;
            if (state[pr(hd)] != motor_state::free) {
                vm = my_velocity(vs,
                        pow(-1, hd)*dot(force, filament_network->get_direction(f_index[hd], l_index[hd])),
                        stall_force);
            }
            this->update_pos_a_end(hd, pos_a_end[hd]+dt*vm); // update relative position
        }
        if (state[hd] == motor_state::bound) this->update_position_attached(hd);  // update absolute position
    }
}

// set pos_a_end relative to the start of the spring the head is bound to
// assumes head is bound
void motor::update_pos_a_end(int hd, double pos)
{
    double spring_length = filament_network->get_llength(f_index[hd],l_index[hd]);
    if (pos >= spring_length) { // "passed" the spring
        if (l_index[hd] == 0) { // the barbed end of the filament
            at_barbed_end[hd] = true;
            pos_a_end[hd] = spring_length;
        } else {
            // Move the motor to the next spring on the filament
            // At the projected new position along that filament*/
            this->set_l_index(hd, l_index[hd] - 1);
            pos_a_end[hd] = pos - spring_length;
        }
    } else if (pos < 0) { //this shouldn't be possible if vm > 0
        if (l_index[hd] == (filament_network->get_filament(f_index[hd])->get_nsprings() - 1)) {
            // the pointed end of the filament
            pos_a_end[hd] = 0; //move head to pointed end
        } else {
            // Move the motor to the previous spring on the filament
            // At the projected new position along that filament
            this->set_l_index(hd, l_index[hd] + 1);
            pos_a_end[hd] = pos + filament_network->get_llength(f_index[hd],l_index[hd]);
        }
    } else {
        pos_a_end[hd] = pos;
    }
}

// updates head position from filament
// assumes head is bound
void motor::update_position_attached(int hd)
{
    vec_type pos = filament_network->get_end(f_index[hd], l_index[hd])
        - pos_a_end[hd] * filament_network->get_direction(f_index[hd], l_index[hd]);
    h[hd] = bc->pos_bc(pos);
}

// add force to filament bound to head
// assumes that head is bound
// Using the lever rule to propagate force as outlined in Nedelec F 2002
void motor::filament_update_hd(int hd, vec_type f)
{
    double pos_ratio = pos_a_end[hd]/filament_network->get_llength(f_index[hd], l_index[hd]);
    filament_network->update_forces(f_index[hd], l_index[hd],   f * pos_ratio);
    filament_network->update_forces(f_index[hd], l_index[hd]+1, f * (1 - pos_ratio));
}

// add force to filaments bound to heads
// heads may be in any state
void motor::filament_update()
{
    if (state[0] == motor_state::bound) this->filament_update_hd(0, force);
    if (state[1] == motor_state::bound) this->filament_update_hd(1, -force);
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
    this->set_l_index(hd, -1);
    f_index[hd]=-1;
    pos_a_end[hd]=0;
}

array<int, 2> motor::get_f_index(){
    return f_index;
}


array<int, 2> motor::get_l_index(){
    return l_index;
}


array<double, 2> motor::get_pos_a_end(){
    return pos_a_end;
}

vec_type motor::get_force()
{
    return force;
}

double motor::get_stretching_energy()
{
    return abs2(force) / (2.0 * mk);
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

            "distance from end of spring = ({}, {})\t "
            "tension = ({}, {})\n",

        h[0].x, h[0].y, h[1].x, h[1].y,
        static_cast<int>(state[0]),  static_cast<int>(state[1]),
        f_index[0],  f_index[1], l_index[0], l_index[1],
        vs, max_bind_dist, mk, stall_force, mld,
        kon, koff, kend, dt, temperature, damp,
        pos_a_end[0], pos_a_end[1], force.x, force.y);
}

vector<double> motor::output()
{
    return {h[0].x, h[0].y, disp.x, disp.y,
        double(f_index[0]), double(f_index[1]),
        double(l_index[0]), double(l_index[1])};
}

void motor::init_l_index(int hd, int idx)
{
    l_index[hd] = idx;
//    if (idx != -1)
    if (state[hd] == motor_state::bound)
        this->add_to_spring(hd);
}

void motor::set_f_index(int hd, int idx)
{
    f_index[hd] = idx;
}

void motor::set_l_index(int hd, int idx)
{
    if (l_index[hd] == -1) {
        // initially unbound, then binds (l_index == -1) ==> add_to_spring
        l_index[hd] = idx;
        this->add_to_spring(hd);
    } else if (idx == -1) {
        // initially bound, unbinds (idx = -1) ==> remove_from_spring
        this->remove_from_spring(hd);
        l_index[hd] = idx;
    } else {
        // initially bound, switches (otherwise) ==> both
        this->remove_from_spring(hd);
        l_index[hd] = idx;
        this->add_to_spring(hd);
    }
}

void motor::set_pos_a_end(int hd, double pos)
{
    pos_a_end[hd] = pos;
}

double motor::get_pos_a_end(int hd)
{
    return pos_a_end[hd];
}

void motor::add_to_spring(int hd)
{
    filament_network->get_filament(f_index[hd])->get_spring(l_index[hd])->add_mot(this, hd);
}

void motor::remove_from_spring(int hd)
{
    filament_network->get_filament(f_index[hd])->get_spring(l_index[hd])->remove_mot(this);
}

string motor::write()
{
    return fmt::format("\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            h[0].x, h[0].y,
            disp.x, disp.y,
            f_index[0], f_index[1],
            l_index[0], l_index[1]);
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

void motor::inc_l_index(int hd)
{
    // NOTE: this function DOES NOT add the motor to a different spring;
    // it just increments the l_index of the spring
    // in cases where the spring pointer hasn't changed
    l_index[hd] += 1;
}

void motor::identify()
{
    cout<<"\nI am a motor";
}

void motor::set_binding_two(double ron2, double roff2, double rend2){
    kon2  = ron2*dt;
    koff2 = roff2*dt;
    kend2 = rend2*dt;
}
