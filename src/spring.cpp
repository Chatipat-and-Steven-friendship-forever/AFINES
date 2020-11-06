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

#include "box.h"

spring::spring(box *bc_, double l0_, double kl_)
{
    bc = bc_;
    kl = kl_;
    l0 = l0_;
    llen = l0;
}

vec_type spring::get_h0()
{
    return h0;
}

vec_type spring::get_h1()
{
    return h1;
}

void spring::step(vec_type h0_, vec_type h1_)
{
    h0 = h0_;
    h1 = h1_;

    disp = bc->rij_bc(h1 - h0);
    llen = abs(disp);

    direc.zero();
    if (llen != 0.0) direc = disp / llen;

    double kf = kl * (llen - l0);
    force = kf * direc;  // -dU/dh0
}

vec_type spring::get_force()
{
    return force;
}

vec_type spring::get_disp()
{
    return disp;
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

double spring::get_length(){
    return llen;
}

// compares parameters only
bool spring::operator==(const spring& that)
{
    return this->kl == that.kl && this->l0 == that.l0;
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

vec_type spring::get_direction()
{
    return direc;
}

double spring::get_stretching_energy()
{
    return 0.5 * kl * (llen - l0) * (llen - l0);
}

virial_type spring::get_virial()
{
    return 0.5 * outer(disp, force);
}
