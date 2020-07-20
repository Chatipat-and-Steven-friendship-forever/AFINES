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
    kb = bending_stiffness;

    ke_vel = 0.0;
    ke_vir = 0.0;
    spring_l0 = spring_length;

    if (beadvec.size() > 0)
    {
        vector<double> &entry = beadvec[0];
        if (entry.size() != 4) throw runtime_error("Wrong number of arguments in beadvec.");
        beads.push_back(new bead(entry[0], entry[1], entry[2], entry[3]));
        prv_rnds.push_back({{0,0}});
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
            prv_rnds.push_back({{0,0}});

        }
    }

    bd_prefactor = sqrt(temperature/(2*dt*damp));

    this->init_ubend();
    fracture_force_sq = fracture_force*fracture_force;
}

filament::~filament(){

    //cout<<"DELETING FILAMENT\n";
    int nr = beads.size(), nl = springs.size();
    for (int i = 0; i < nr; i ++)
    {
        //cout<<"\nDEBUG: deleting pointer "<<beads[i];
        delete beads[i];
    }
    for (int i = 0; i < nl; i ++)
        delete springs[i];

    beads.clear();
    springs.clear();
    prv_rnds.clear();
}

void filament::add_bead(vector<double> a, double spring_length, double stretching_stiffness, double max_ext_ratio){
    if (a.size() != 4) throw runtime_error("Wrong number of arguments in bead.");
    beads.push_back(new bead({a[0], a[1], a[2], a[3]}));
    prv_rnds.push_back({{0,0}});
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
    double vx, vy, fx, fy, fx_brn, fy_brn, x, y;
    array<double, 2> new_rnds;
    array<double, 2> newpos;
    ke_vel = 0.0;
    ke_vir = 0.0;
    int sa = int(beads.size());
    int la = int(springs.size());
    for (int i = 0; i < sa; i++){
        new_rnds = {{rng_n(), rng_n()}};
        fx = beads[i]->get_force()[0];
        fy = beads[i]->get_force()[1];
        fx_brn = bd_prefactor*damp*(new_rnds[0] + prv_rnds[i][0]);
        fy_brn = bd_prefactor*damp*(new_rnds[1] + prv_rnds[i][1]);
        vx  = fx/damp  + fx_brn/damp;
        vy  = fy/damp  + fy_brn/damp;
        //        cout<<"\nDEBUG: Fx("<<i<<") = "<<beads[i]->get_force()[0]<<"; v = ("<<vx<<" , "<<vy<<")";
        x = beads[i]->get_xcm();
        y = beads[i]->get_ycm();
        prv_rnds[i] = new_rnds;
        //cout<<"\nDEBUG: bead force = ("<<beads[i]->get_force()[0]<<" , "<<beads[i]->get_force()[1]<<")";
        ke_vel += vx*vx + vy*vy;
        ke_vir += -(0.5)*((fx*x + fy*y) + (fx_brn*x + fy_brn*y));
        newpos = bc->pos_bc({x + vx*dt, y + vy*dt});
        beads[i]->set_xcm(newpos[0]);
        beads[i]->set_ycm(newpos[1]);
        beads[i]->reset_force();
    }

    for (int i = 0; i < la; i++)
        springs[i]->step();

}

vector<filament *> filament::update_stretching(double t)
{
    vector<filament *> newfilaments;
    array<double, 2> spring_force;
    if(springs.size() == 0)
        return newfilaments;

    for (unsigned int i=0; i < springs.size(); i++) {
        springs[i]->update_force();
        spring_force = springs[i]->get_force();
        if (spring_force[0]*spring_force[0]+spring_force[1]*spring_force[1] > fracture_force_sq){
            newfilaments = this->fracture(i);
            break;
        }
        else
            springs[i]->filament_update();
    }

    return newfilaments;
}

spring *filament::get_spring(int i)
{
    return springs[i];
}

void filament::update_d_strain(double g)
{
    for (size_t i = 0; i < beads.size(); i++) {
        beads[i]->set_xcm(beads[i]->get_xcm() + g * beads[i]->get_ycm() / bc->get_ybox());
    }
}

box *filament::get_box()
{
    return bc;
}

void filament::update_forces(int index, double f1, double f2)
{
    beads[index]->update_force(f1,f2);
}

void filament::pull_on_ends(double f)
{
    if (beads.size() < 2) return;
    int last = beads.size() - 1;
    array<double, 2> dr = bc->rij_bc({beads[last]->get_xcm() - beads[0]->get_xcm(),
                                      beads[last]->get_ycm() - beads[0]->get_ycm()});
    double len = hypot(dr[0], dr[1]);

    beads[ 0  ]->update_force(-0.5*f*dr[0]/len, -0.5*f*dr[1]/len);
    beads[last]->update_force( 0.5*f*dr[0]/len,  0.5*f*dr[1]/len);
}

void filament::affine_pull(double f)
{
    if (beads.size() < 2) return;
    int last = beads.size() - 1;
    array<double, 2> dr = bc->rij_bc({beads[last]->get_xcm() - beads[0]->get_xcm(),
                                      beads[last]->get_ycm() - beads[0]->get_ycm()});
    double len = hypot(dr[0], dr[1]);
    //cout<<"\nDEBUG: angle = "<<ang;
    double frac, fcos = f*dr[0]/len, fsin = f*dr[1]/len;

    for (int i = 0; i <= last; i++){
        frac = (double(i)/double(last)-0.5);
        beads[i]->update_force(frac*fcos, frac*fsin);
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

string filament::write_beads(int fil){
    string all_beads;
    for (unsigned int i =0; i < beads.size(); i++)
    {
        all_beads += beads[i]->write() + "\t" + std::to_string(fil);
    }

    return all_beads;
}

string filament::write_springs(int fil){
    string all_springs;
    for (unsigned int i =0; i < springs.size(); i++)
    {
        all_springs += springs[i]->write() + "\t" + std::to_string(fil);
    }

    return all_springs;
}

string filament::write_thermo(int fil)
{
    return
        "\n" + std::to_string(this->get_kinetic_energy_vel()) +
        "\t" + std::to_string(this->get_kinetic_energy_vir()) +
        "\t" + std::to_string(this->get_potential_energy()) +
        "\t" + std::to_string(this->get_total_energy()) +
        "\t" + std::to_string(fil);
}

vector<vector<double>> filament::get_beads(unsigned int first, unsigned int last)
{
    vector<vector<double>> newbeads;
    for (size_t i = first; i < last; i++) {
        if (i >= beads.size()) {
            break;
        } else {
            bead *b = beads[i];
            newbeads.push_back({
                    b->get_xcm(), b->get_ycm(),
                    b->get_length(), b->get_viscosity()});
        }
    }
    return newbeads;
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

string filament::to_string(){

    // Note: not including springs in to_string, because spring's to_string includes filament's to_string
    char buffer[200];
    string out = "\n";

    for (unsigned int i = 0; i < beads.size(); i++)
        out += beads[i]->to_string();

    sprintf(buffer, "temperature = %f\tdt = %f\tfracture_force=%f\n",
            temperature, dt, fracture_force);

    return out + buffer;

}

inline double filament::angle_between_springs(int i, int j){

    array<double, 2> delr1, delr2;
    double r1,r2, c;

    // 1st bond
    delr1 = springs[i]->get_disp();
    r1    = springs[i]->get_length();

    // 2nd bond
    delr2 = springs[j]->get_disp();
    r2    = springs[j]->get_length();

    // cos angle
    c = (delr1[0]*delr2[0] + delr1[1]*delr2[1]) / (r1*r2);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    return acos(c);
}


void filament::update_bending(double t)
{
    if (springs.size() <= 1 || kb == 0) return;

    virial_clear(bending_virial);
    ubend = 0.0;
    for (int n = 0; n < int(springs.size())-1; n++) {

        array<double, 2> delr1 = springs[n]->get_disp();
        array<double, 2> delr2 = springs[n+1]->get_disp();

        bend_result_type result = bend_harmonic(kb, 0.0, delr1, delr2);

        ubend += result.energy;

        // apply force to each of 3 atoms
        beads[n  ]->update_force(-result.force1[0], -result.force1[1]);
        beads[n+1]->update_force(result.force1[0], result.force1[1]);

        beads[n+1]->update_force(-result.force2[0], -result.force2[1]);
        beads[n+2]->update_force(result.force2[0], result.force2[1]);

        bending_virial[0][0] += result.force1[0] * delr1[0] + result.force2[0] * delr2[0];
        bending_virial[0][1] += result.force1[0] * delr1[1] + result.force2[0] * delr2[1];
        bending_virial[1][0] += result.force1[1] * delr1[0] + result.force2[1] * delr2[0];
        bending_virial[1][1] += result.force1[1] * delr1[1] + result.force2[1] * delr2[1];
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

array<array<double, 2>, 2> filament::get_bending_virial()
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

array<array<double, 2>, 2> filament::get_stretching_virial()
{
    array<array<double, 2>, 2> vir;
    virial_clear(vir);
    for (spring *s : springs) {
        virial_add(vir, s->get_virial());
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

array<double, 2> filament::get_bead_position(int n)
{
    return {{beads[n]->get_xcm(), beads[n]->get_ycm()}};
}

void filament::print_thermo()
{
    cout
        << "\tKE = "<< this->get_kinetic_energy_vel()
        << "\tKEvir = "<< this->get_kinetic_energy_vir()
        << "\tPE = "<< this->get_potential_energy()
        << "\tTE = "<< this->get_total_energy();
}

double filament::get_end2end()
{
    if (beads.size() < 2) {
        return 0;
    } else {
        return bc->dist_bc({beads[beads.size() - 1]->get_xcm() - beads[0]->get_xcm(),
                            beads[beads.size() - 1]->get_ycm() - beads[0]->get_ycm()});
    }
}

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
        double x2, y2, pos;
        array<double, 2> dir = springs[0]->get_direction();
        x2 = beads[1]->get_xcm();
        y2 = beads[1]->get_ycm();
        //add a bead
        array<double, 2> newpos = bc->pos_bc({x2-spring_l0*dir[0], y2-spring_l0*dir[1]});
        beads.insert(beads.begin()+1, new bead(newpos[0], newpos[1], beads[0]->get_length(), beads[0]->get_viscosity()));
        prv_rnds.insert(prv_rnds.begin()+1, {{0, 0}});
        //shift all springs from "1" onward forward
        //move backward; otherwise i'll just keep pushing all the motors to the pointed end, i think
        for (int i = int(springs.size()-1); i > 0; i--){
            springs[i]->inc_aindex();
            springs[i]->step();
            //shift all xsprings on these springs forward
            //lmots = springs[i]->get_mots();
            for (map<motor *, int>::iterator it = springs[i]->get_mots().begin(); it != springs[i]->get_mots().end(); ++it)
            {
                it->first->inc_l_index(it->second);
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
        for (map<motor *, int>::iterator it = mots0.begin(); it != mots0.end(); ++it)
        {
            pos = it->first->get_pos_a_end()[it->second];
            if (pos < spring_l0){
                mots1.push_back(it->first);
            }
            else{
                it->first->set_pos_a_end(it->second, pos - spring_l0);
            }
        }
        int hd;
        for (unsigned int i = 0; i < mots1.size(); i++)
        {
            hd = mots0[mots1[i]];
            mots1[i]->set_l_index(hd, 1);
        }
    }
}

void filament::update_length()
{
    if ( kgrow*lgrow > 0 && this->get_nsprings() + 1 <= nsprings_max && rng(0,1) < kgrow*dt){
        grow(lgrow);
    }
}

void filament::set_kgrow(double k){
    kgrow = k;
}

void filament::set_lgrow(double dl){
    lgrow = dl;
}
