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
    force       = {{0,0}}; // force on the spring
    ke_vel = 0; //assume m = 1
    ke_vir = 0;
    pos_a_end = {{0, 0}}; // pos_a_end = distance from pointy end -- by default 0
                        // i.e., if l_index[hd] = j, then pos_a_end[hd] is the distance to the "j+1"th bead

    array<double, 2> posH0 = boundary_check(0, mvec[0], mvec[1]);
    array<double, 2> posH1 = boundary_check(1, mvec[0] + mvec[2], mvec[1] + mvec[3]);
    hx[0] = posH0[0];
    hy[0] = posH0[1];
    hx[1] = posH1[0];
    hy[1] = posH1[1];
    //force can be non-zero and angle is determined from disp vector
    this->update_angle();
    this->update_force();

    ldir_bind[0] = {{0,0}};
    ldir_bind[1] = {{0,0}};
    bind_disp[0] = {{0,0}};
    bind_disp[1] = {{0,0}};

    at_barbed_end = {{false, false}};

    if (state[0] == motor_state::bound){
        pos_a_end[0] = bc->dist_bc({filament_network->get_end(f_index[0], l_index[0])[0] - hx[0],
                                    filament_network->get_end(f_index[0], l_index[0])[1] - hy[0]});
        ldir_bind[0] = filament_network->get_direction(f_index[0], l_index[0]);

    }
    if (state[1] == motor_state::bound){
        pos_a_end[1] = bc->dist_bc({filament_network->get_end(f_index[1], l_index[1])[0] - hx[1],
                                    filament_network->get_end(f_index[1], l_index[1])[1] - hy[1]});
        ldir_bind[1] = filament_network->get_direction(f_index[1], l_index[1]);
    }

    prv_rnd_x = {{0,0}};
    prv_rnd_y = {{0,0}};

}

//return motor state with a given head number
array<motor_state, 2> motor::get_states()
{
    return state;
}

array<double, 2> motor::get_hx()
{
    return hx;
}


array<double, 2> motor::get_hy()
{
    return hy;
}

//metropolis algorithm with rate constant
//NOTE: while fl_idx doesn't matter for this xlink implementation, it does for "spacers"
double motor::metropolis_prob(int hd, array<int, 2> fl_idx, array<double, 2> newpos, double maxprob)
{
    double prob = maxprob;
    double stretch  = bc->dist_bc({newpos[0] - hx[pr(hd)], newpos[1] - hy[pr(hd)]}) - mld;
    double delE = 0.5*mk*stretch*stretch - this->get_stretching_energy();

    if( delE > 0 )
        prob *= exp(-delE/temperature);

    return prob;
}

bool motor::allowed_bind(int hd, array<int, 2> fl_idx){
    return (f_index[pr(hd)] != fl_idx[0] || l_index[pr(hd)] != fl_idx[1]);
}

void motor::attach_head(int hd, array<double, 2> intpoint, array<int, 2> fl)
{
    // update state
    state[hd] = motor_state::bound;
    f_index[hd] = fl[0];
    this->set_l_index(hd, fl[1]);

    // record displacement of head and orientation of spring for future unbinding move
    ldir_bind[hd] = filament_network->get_direction(f_index[hd], l_index[hd]);
    bind_disp[hd] = bc->rij_bc({intpoint[0] - hx[hd], intpoint[1] - hy[hd]});

    // update head position
    hx[hd] = intpoint[0];
    hy[hd] = intpoint[1];

    // update relative head position
    array<double, 2> ref = filament_network->get_end(f_index[hd], l_index[hd]);
    pos_a_end[hd] = bc->dist_bc({ref[0] - hx[hd], ref[1] - hy[hd]});

    // even if the head is at the barbed end upon binding,
    // it could have negative velocity,
    // so always set this to false until it steps
    at_barbed_end[hd] = false;
}

bool motor::attach_opt(int hd)
{
    double mf_rand = rng_u();
    double onrate = (state[pr(hd)] == motor_state::bound) ? kon2 : kon;
    vector<array<int, 2>> *attach_list = filament_network->get_attach_list(hx[hd], hy[hd]);
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
        array<double ,2> intpoint = s->intpoint({hx[hd], hy[hd]});

        // don't bind if binding site is further away than the cutoff
        array<double, 2> dr = bc->rij_bc({intpoint[0] - hx[hd], intpoint[1] - hy[hd]});
        double dist_sq = dr[0] * dr[0] + dr[1] * dr[1];
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

//check for attachment of unbound heads given head index (0 for head 1, and 1 for head 2)
bool motor::attach(int hd)
{
    set<tuple<double, array<int, 2>, array<double, 2>>> dist_sq_sorted;
    for (array<int, 2> fl : *filament_network->get_attach_list(hx[hd], hy[hd])) {
        filament *f = filament_network->get_filament(fl[0]);
        spring *s = f->get_spring(fl[1]);
        array<double, 2> intpoint = s->intpoint({hx[hd], hy[hd]});
        array<double, 2> dr = bc->rij_bc({intpoint[0] - hx[hd], intpoint[1] - hy[hd]});
        double dist_sq = dr[0] * dr[0] + dr[1] * dr[1];
        dist_sq_sorted.insert(tuple<double, array<int, 2>, array<double, 2>>(dist_sq, fl, intpoint));
    }

    double not_off_prob = 0.0;
    double mf_rand = rng_u();
    double onrate = (state[pr(hd)] == motor_state::bound) ? kon2 : kon;
    for (auto it : dist_sq_sorted) {
        double dist_sq; array<int, 2> fl; array<double, 2> intpoint;
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
    force = {{tension*direc[0], tension*direc[1]}};
//    force = {{tension*disp[0]/len, tension*disp[1]/len}};
}

/* Taken from hsieh, jain, larson, jcp 2006; eqn (5)
 * Adapted by placing a cutoff, similar to how it's done in LAMMPS src/bond_fene.cpp*/

void motor::update_force_fraenkel_fene()
{
    double ext = abs(mld - len);
    double scaled_ext, mkp;

    if (max_ext - ext > eps_ext )
        scaled_ext = ext/max_ext;
    else
        scaled_ext = (max_ext - eps_ext)/max_ext;

    mkp = mk/(1-scaled_ext*scaled_ext)*(len-mld);
    force = {{mkp*direc[0], mkp*direc[1]}};

}


void motor::brownian_relax(int hd)
{

    double new_rnd_x= rng_n(), new_rnd_y = rng_n();

    double vx =  pow(-1,hd)*force[0] / damp + bd_prefactor*(new_rnd_x + prv_rnd_x[hd]);
    double vy =  pow(-1,hd)*force[1] / damp + bd_prefactor*(new_rnd_y + prv_rnd_y[hd]);
    ke_vel = vx*vx + vy*vy;
    ke_vir = -(0.5)*(pow(-1,hd))*(force[0]*hx[hd] + force[1]*hy[hd]);
    array<double, 2> pos = boundary_check(hd, hx[hd] + vx*dt, hy[hd] + vy*dt);
    hx[hd] = pos[0];
    hy[hd] = pos[1];

    prv_rnd_x[hd] = new_rnd_x;
    prv_rnd_y[hd] = new_rnd_y;

}

void motor::kill_head(int hd)
{
    state[hd] = motor_state::dead;
}

void motor::deactivate_head(int hd)
{
    state[hd] = motor_state::inactive;
}

void motor::relax_head(int hd)
{
    array<double, 2> newpos = boundary_check(hd, hx[pr(hd)] - pow(-1, hd)*mld*direc[0], hy[pr(hd)] - pow(-1, hd)*mld*direc[1]);
    hx[hd] = newpos[0];
    hy[hd] = newpos[1];
}


void motor::update_angle()
{
    disp  = bc->rij_bc({hx[1]-hx[0], hy[1]-hy[0]});
    len   = hypot(disp[0], disp[1]);
    if ( len != 0 )
        direc = {{disp[0]/len, disp[1]/len}};
    else
        direc = {{0, 0}};
}

array<double, 2> motor::boundary_check(int hd, double x, double y)
{
    return bc->pos_bc({x, y}, {(x - hx[hd])/dt, (y - hy[hd])/dt}, dt);
}

array<double, 2> motor::generate_off_pos(int hd){
    array<double, 2> ldir = filament_network->get_direction(f_index[hd], l_index[hd]);
    double c = dot(  ldir, ldir_bind[hd]);
    double s = cross(ldir, ldir_bind[hd]);

    array<double, 2> bind_disp_rot = {{bind_disp[hd][0]*c - bind_disp[hd][1]*s, bind_disp[hd][0]*s + bind_disp[hd][1]*c}};

    return bc->pos_bc(
            {hx[hd] - bind_disp_rot[0], hy[hd] - bind_disp_rot[1]},
            {-bind_disp_rot[0]/dt, -bind_disp_rot[1]/dt}, dt);
}


//stepping and detachment kinetics of a single bound head
void motor::step_onehead(int hd)
{

    array<double, 2> hpos_new = generate_off_pos(hd);
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


void motor::update_pos_a_end(int hd, double pos)
{
//    cout<<"\nDEBUG: new pos = "<<pos;
    double spring_length = filament_network->get_llength(f_index[hd],l_index[hd]);
    if (pos >= spring_length) { // "passed" the spring
        if (l_index[hd] == 0){ // the barbed end of the filament
            at_barbed_end[hd] = true;
            pos_a_end[hd] = spring_length;
        }
        else{
            /*Move the motor to the next spring on the filament
             *At the projected new position along that filament*/
            this->set_l_index(hd, l_index[hd]-1);
            pos_a_end[hd] = pos - spring_length;
        }
    }
    else if (pos < 0) { //this shouldn't be possible if vm > 0
        if (l_index[hd] == (filament_network->get_filament(f_index[hd])->get_nsprings() - 1)){ // the pointed end of the filament
            pos_a_end[hd]=0; //move head to pointed end
        }
        else{
            /*Move the motor to the previous spring on the filament
             *At the projected new position along that filament*/
            this->set_l_index(hd, l_index[hd] + 1);
            pos_a_end[hd] = pos + filament_network->get_llength(f_index[hd],l_index[hd]);
        }
    }
    else {
        pos_a_end[hd] = pos;
    }

}


void motor::update_position_attached(int hd){

    double posx = filament_network->get_end(f_index[hd],l_index[hd])[0]-pos_a_end[hd]*filament_network->get_direction(f_index[hd],l_index[hd])[0];
    double posy = filament_network->get_end(f_index[hd],l_index[hd])[1]-pos_a_end[hd]*filament_network->get_direction(f_index[hd],l_index[hd])[1];

    array<double, 2> newpos = boundary_check(hd, posx, posy);

    hx[hd] = newpos[0];
    hy[hd] = newpos[1];

}

// Using the lever rule to propagate force as outlined in Nedelec F 2002

void motor::filament_update_hd(int hd, array<double, 2> f)
{
    double pos_ratio = pos_a_end[hd]/filament_network->get_llength(f_index[hd], l_index[hd]);
    filament_network->update_forces(f_index[hd], l_index[hd],   f[0] *    pos_ratio , f[1] *    pos_ratio );
    filament_network->update_forces(f_index[hd], l_index[hd]+1, f[0] * (1-pos_ratio), f[1] * (1-pos_ratio));
}

void motor::filament_update()
{
    if (state[0] == motor_state::bound) this->filament_update_hd(0, force);
    if (state[1] == motor_state::bound) this->filament_update_hd(1, {{-force[0], -force[1]}});
}

void motor::detach_head(int hd)
{
    array<double, 2> hpos_new = generate_off_pos(hd);
    this->detach_head(hd, hpos_new);
}

void motor::detach_head(int hd, array<double, 2> newpos)
{
    detach_head_without_moving(hd);
    hx[hd] = newpos[0];
    hy[hd] = newpos[1];
}

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


array<double, 2> motor::get_force(){
    return force;
}


double motor::get_stretching_energy(){
    return (force[0]*force[0]+force[1]*force[1])/(2*mk);
}

array<array<double, 2>, 2> motor::get_virial() {
    double k = tension / len;
    return {
        array<double, 2>{k * disp[0] * disp[0], k * disp[0] * disp[1]},
        array<double, 2>{k * disp[0] * disp[1], k * disp[1] * disp[1]}
    };
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
    char buffer[1000];
    string out ="";

    sprintf(buffer, "\
            \nhead 0 position = (%f, %f)\t head 1 position=(%f,%f)\
            \nstate = (%d, %d)\t f_index = (%d, %d)\t l_index = (%d, %d)\
            \nviscosity = %f\t max binding distance = %f\t stiffness = %f\t stall force = %f\t length = %f\
            \nkon = %f\t koff = %f\t kend = %f\t dt = %f\t temp = %f\t damp = %f\
            \ndistance from end of spring = (%f, %f)\t tension = (%f, %f)\n",
            hx[0], hy[0], hx[1], hy[1],
            static_cast<int>(state[0]),  static_cast<int>(state[1]),
            f_index[0],  f_index[1], l_index[0], l_index[1],
            vs, max_bind_dist, mk, stall_force, mld,
            kon, koff, kend, dt, temperature, damp,
            pos_a_end[0], pos_a_end[1], force[0], force[1]);
    return buffer;
}

vector<double> motor::output()
{
    return {hx[0], hy[0], disp[0], disp[1],
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
    /* cases:
            initially unbound, then binds (l_index == -1) ==> add_to_spring
            initially bound, unbinds (idx = -1) ==> remove_from_spring
            initially bound, switches (otherwise) ==> both
    */
//    cout<<"\nDEBUG: changing l_index["<<hd<<"] of "<<this<<" from "<<l_index[hd]<<" to "<<idx;

    if (l_index[hd] == -1){
        l_index[hd] = idx;
        this->add_to_spring(hd);
    }
    else if(idx == -1){
        this->remove_from_spring(hd);
        l_index[hd] = idx;
    }
    else{
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
    return "\n" + std::to_string(hx[0]) + "\t" + std::to_string(hy[0])
        +  "\t" + std::to_string(disp[0]) + "\t" + std::to_string(disp[1])
        +  "\t" + std::to_string(f_index[0]) + "\t" + std::to_string(f_index[1])
        +  "\t" + std::to_string(l_index[0]) + "\t" + std::to_string(l_index[1]);
}

void motor::revive_head(int hd)
{
    state[hd] = motor_state::free;
}

void motor::update_d_strain(double g)
{
    array<double, 2> fov = bc->get_fov();
    array<double, 2> pos0 = boundary_check(0, hx[0] + g * hy[0] / fov[1], hy[0]);
    array<double, 2> pos1 = boundary_check(1, hx[1] + g * hy[1] / fov[1], hy[1]);
    hx[0] = pos0[0];
    hx[1] = pos1[0];
    hy[0] = pos0[1];
    hy[1] = pos1[1];
}

void motor::inc_l_index(int hd){
//    this->set_l_index(hd, l_index[hd]+1);
/* NOTE: this function DOES NOT add the motor to a different spring; it just increments the l_index of the spring
 * in cases where the spring pointer hasn't changed*/
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
