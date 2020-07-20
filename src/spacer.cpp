/*
 * spacer.cpp
 *  
 *
 *  Created by Simon Freedman
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

#include "spacer.h"
#include "filament_ensemble.h"
#include "potentials.h"
//spacer class

spacer::spacer(vector<double> mvec,
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

    max_bind_dist = rcut;
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
    array<int, 2> mylindex = {int(mvec[6]), int(mvec[7])};

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
    b_eng       = {{0,0}}; // filament / xlink bending energy

    b_force[0].zero(); //b_force[0] = bending force on head 0 due to h0-h1-spring angle in cartesian coords
    b_force[1].zero(); //b_force[1] = bending force on head 1 due to h1-h0-spring angle in cartesian coords

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

spacer::~spacer(){}

void spacer::set_bending(double force_constant, double ang){
    kb  = force_constant;
    th0 = ang;
}

void spacer::update_force()
{ 
    //cout<<"\nDEBUG: using spacer update_force()";
    if (state[0] == motor_state::bound && state[1] == motor_state::bound){
        update_bending(0);
        update_bending(1);
    }
    else
    {
        b_force[0].zero();
        b_force[1].zero();
    }

    update_angle(); // need to recalculate this, since heads might've moved

    tension = mk*(len - mld);
    force = tension * direc;
}
  //Measure distance to FARTHER END of bead filament that the spacer is bound to
  //
int spacer::get_further_end(int hd, int findex, int lindex)
{
    return (pos_a_end[hd] > 0.5*filament_network->get_llength(findex, lindex));
}

vec_type spacer::disp_from_bead(int hd, int findex, int aindex)
{
    vec_type pos = filament_network->get_filament(findex)->get_bead_position(aindex);
    return bc->rij_bc(pos - h[hd]);
}

void spacer::update_bending(int hd)
{
    int bead_further_end = get_further_end(hd, f_index[hd], l_index[hd]);

    vec_type delr1 = disp_from_bead(hd, f_index[hd], l_index[hd] + bead_further_end);
    vec_type delr2 = pow(-1, hd)*disp;

    bend_result_type result = bend_harmonic(kb, th0, delr1, delr2);

    // apply force to each of 3 atoms
    filament_network->update_forces(f_index[hd], l_index[hd] + bead_further_end, result.force1);
    b_force[hd] -= result.force1;

    b_force[pr(hd)] += result.force2;
    b_force[hd] -= result.force2;

    b_eng[hd] = result.energy;
}

double spacer::get_kb(){
    return kb;
}

double spacer::get_th0(){
    return th0;
}

void spacer::identify(){
    cout<<"I am a spacer";
}

void spacer::brownian_relax(int hd)
{
    vec_type new_rnd {rng_n(), rng_n()};
    vec_type v = pow(-1, hd) * force + b_force[hd] / damp + bd_prefactor * (new_rnd + prv_rnd[hd]);
    ke_vel = abs2(v);
    ke_vir = -0.5 * pow(-1, hd) * dot(force, h[hd]);
    h[hd] = bc->pos_bc(h[hd] + v*dt);
    prv_rnd[hd] = new_rnd;
}

void spacer::filament_update()
{
    if (state[0]==motor_state::bound) this->filament_update_hd(0, force + b_force[0]);
    if (state[1]==motor_state::bound) this->filament_update_hd(1, -force + b_force[1]);

    //reset bending force
    b_force[0].zero();
    b_force[1].zero();
}

array<vec_type, 2> spacer::get_b_force()
{
    return b_force;
}


//metropolis algorithm with rate constant
double spacer::metropolis_prob(int hd, array<int, 2> fl_idx, vec_type newpos, double maxprob)
{
    double prob = maxprob;

    double stretch  = bc->dist_bc(newpos - h[pr(hd)]) - mld;
    double delEs = 0.5 * mk * stretch * stretch - get_stretching_energy();

    double bend_eng = 0.0;
    if (state[hd] == motor_state::free && state[pr(hd)] == motor_state::bound) { //it's trying to attach
        vec_type delr1 = disp_from_bead(hd, fl_idx[0], fl_idx[1] + get_further_end(hd, fl_idx[0], fl_idx[1]));
        vec_type delr2 = pow(-1, hd)*disp;
        bend_eng = bend_harmonic_energy(kb, th0, delr1, delr2);
    }
    double delEb = bend_eng - b_eng[hd];

    double delE = delEs + delEb;
    if (delE > 0.0) prob *= exp(-delE/temperature);

    return prob;
}

bool spacer::allowed_bind(int hd, array<int, 2> fl_idx)
{
//    cout<<"\nDEBUG: using spacer allowed bind";
    return (fl_idx[0] != f_index[pr(hd)]);
}

array<double, 2> spacer::get_bending_energy()
{
    return b_eng;
}
