/*------------------------------------------------------------------
 spring.cpp : object describing a hookean spring that connects beads
            in a worm-like chain.

 Copyright (C) 2016
 Created by: Simon Freedman, Shiladitya Banerjee, Glen Hocky, Aaron Dinner
 Contact: dinner@uchicago.edu

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version. See ../LICENSE for details.
-------------------------------------------------------------------*/

#include "spring.h"
#include "globals.h"
#include "filament.h"

spring::spring(){ }

spring::spring(double len, double stretching_stiffness, double max_ext_ratio, filament* f,
        array<int, 2> myaindex)
{
    bc = f->get_box();
    kl      = stretching_stiffness;
    l0      = len;
    fil     = f;
    aindex  = myaindex;

    max_ext = max_ext_ratio * l0;
    eps_ext = 0.01*max_ext;

    h0 = {0, 0};
    h1 = {0, 0};

    force = {0, 0};
    llensq = l0*l0;
    llen = l0;
}

spring::~spring(){
}

vec_type spring::get_h0()
{
    return h0;
}

vec_type spring::get_h1()
{
    return h1;
}

// stepping kinetics

void spring::step()
{
    vec_type h0 = fil->get_bead_position(aindex[0]);
    vec_type h1 = fil->get_bead_position(aindex[1]);

    disp   = bc->rij_bc(h1 - h0);
    llensq = abs2(disp);
    llen   = sqrt(llensq);

    if (llen != 0)
        direc = disp / llen;
    else
        direc = {0, 0};

}

void spring::update_force()
{
    double kf = kl*(llen-l0);
    force = kf * direc;
}

/* Taken from hsieh, jain, larson, jcp 2006; eqn (5)
 * Adapted by placing a cutoff, similar to how it's done in LAMMPS src/bond_fene.cpp*/
void spring::update_force_fraenkel_fene()
{
    double ext = abs(l0 - llen);
    double scaled_ext, klp;
    if (max_ext - ext > eps_ext ){
        scaled_ext = ext/max_ext;
    }
    else{
        scaled_ext = (max_ext - eps_ext)/max_ext;
    }

    klp = kl/(1-scaled_ext*scaled_ext)*(llen-l0);
    force = klp * direc;

}

vec_type spring::get_force()
{
    return force;
}

vec_type spring::get_disp()
{
    return disp;
}

void spring::filament_update()
{
    fil->update_forces(aindex[0],  force);
    fil->update_forces(aindex[1], -force);
}

double spring::get_kl(){
    return kl;
}

void spring::set_l0(double myl0){
    l0 = myl0;
}

double spring::get_l0(){
    return l0;
}

double spring::get_fene_ext(){
    return max_ext/l0;
}

double spring::get_length(){
    return llen;
}

double spring::get_length_sq(){
    return llensq;
}

vector<double> spring::output()
{
    return {h0.x, h0.y, disp.x, disp.y};
}

std::string spring::write()
{
    return "\n" + std::to_string(h0.x) + "\t" + std::to_string(h0.y) + "\t" + std::to_string(disp.x) + "\t"
        + std::to_string(disp.y);
}

std::string spring::to_string(){

    char buffer [100];
    sprintf(buffer, "aindex[0] = %d;\t aindex[1] = %d;\t kl = %f;\t l0 = %f\nfilament : \n",
                        aindex[0], aindex[1], kl, l0);

    return buffer + fil->to_string();

}

bool spring::operator==(const spring& that)
{
    /*Note: you can't compare the filament objects because that will lead to infinite recursion;
     * this function requires the filament poiner to be identical to evaluate to true*/
    return (this->aindex[0] == that.aindex[0] && this->aindex[1] == that.aindex[1] &&
            this->kl == that.kl &&
            this->l0 == that.l0 && this->fil == that.fil);
}

bool spring::is_similar(const spring& that)
{

    /* Same as ==; but doesn't compare the filament pointer*/

    return (this->aindex[0] == that.aindex[0] && this->aindex[1] == that.aindex[1] &&
            this->kl == that.kl &&
            this->l0 == that.l0);
}

//shortest(perpendicular) distance between an arbitrary point and the spring
//SO : 849211
vec_type spring::intpoint(vec_type pos)
{
    double l2 = abs2(disp);
    if (l2 == 0) {
        return h0;
    } else {
        //Consider the line extending the spring, parameterized as h0 + tp ( h1 - h0 )
        //tp = projection of pos onto the line
        double tp = bc->dot_bc(pos - h0, h1 - h0)/l2;
        if (tp < 0.0) {
            return h0;
        } else if (tp > 1.0 ) {
            return h1;
        } else{
            //velocity and dt are 0 since not relevant
            return bc->pos_bc(h0 + tp * disp);
        }
    }
}

bool spring::get_line_intersect(spring *l2)
{
    //Reference to Stack Overflow entry by iMalc on Feb 10, 2013
    //Web Address: https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect

    vec_type disp1 = this->get_disp();
    vec_type disp2 = l2->get_disp();

    vec_type disp12 = bc->rij_bc(h0 - l2->get_h0());

    double denom = disp1.x*disp2.y - disp1.y*disp2.x;
    if (denom == 0) return false;
    bool denomPos = denom > 0;

    double s_num = disp1.x*disp12.y - disp1.y*disp12.x;
    double t_num = disp2.x*disp12.y - disp2.y*disp12.x;

    if ((s_num < 0) == denomPos) return false;
    if ((t_num < 0) == denomPos) return false;

    if (((s_num > denom) == denomPos) || ((t_num > denom) == denomPos)) return false;

    //Else Collision have been detected, the filaments do intersect!
    return true;
}

vec_type spring::get_direction()
{
    return direc;
}

double spring::get_stretching_energy(){
    return abs2(force) / (2.0 * kl);
}

virial_type spring::get_virial() {
    double k = kl*(llen-l0)/llen;
    return outer(disp, k * disp);
}

double spring::get_stretching_energy_fene()
{
    double ext = abs(l0 - llen);

    if (max_ext - ext > eps_ext )
        return -0.5*kl*max_ext*max_ext*log(1-(ext/max_ext)*(ext/max_ext));
    else
        return 0.25*kl*ext*ext*(max_ext/eps_ext);

}

void spring::set_aindex(array<int,2> aidx)
{
    aindex = aidx;
}

double spring::get_max_ext()
{
    return max_ext;
}

array<int, 2> spring::get_aindex(){
    return aindex;
}

// functions for growing
void spring::inc_aindex()
{
    aindex = {{aindex[0]+1, aindex[1]+1}};
}

void spring::add_mot(motor * mot, int hd)
{
    mots[mot] = hd;
}

void spring::remove_mot(motor * mot)
{
    //cout<<"\nDEBUG: trying to remove motor";
///    cout<<"\nDEBUG: removing mot "<<mot<<" from mots on spring "<<this;
    mots.erase(mot);
}

map<motor *, int> & spring::get_mots()
{
    map<motor *, int> &ptr = mots;
    return ptr;
}

/*
int spring::get_n_mots(){
    return int(mots.size());
}

motor * spring::get_mot(int i){
    return mots[i];
}

int spring::get_mot_hd(int i){
    return mot_hds[i];
}*/
