/*------------------------------------------------------------------
 bead.cpp : object describing a circular bead
 
 Copyright (C) 2016 
 Created by: Simon Freedman, Shiladitya Banerjee, Glen Hocky, Aaron Dinner
 Contact: dinner@uchicago.edu

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version. See ../LICENSE for details. 
-------------------------------------------------------------------*/

#include "bead.h"
#include "globals.h"

bead::bead(double xcm, double ycm, double len, double vis)
{
    pos = {xcm, ycm};
    rad=len; //radius
    visc=vis;
    friction = 6*pi*visc*rad;
    force = {0.0, 0.0};
}

bead::bead(const bead& other)
{
    pos = other.pos;
    rad = other.rad;
    visc = other.visc;
    friction = other.friction;
    force = other.force;
}

vec_type bead::get_force()
{
    return force;
}

double bead::get_length()
{
    return rad;
}

void bead::update_force(vec_type f)
{
    force += f;
}

void bead::reset_force()
{
    force.x = force.y = 0.0;
}

vec_type bead::get_pos()
{
    return pos;
}

void bead::set_pos(vec_type p)
{
    pos = p;
}

bool bead::operator==(const bead& that)
{
    double err = eps;
    return close(pos.x , that.pos.x , err)
        && close(pos.y , that.pos.y , err)
        && close(rad , that.rad , err)
        && close(visc , that.visc , err)
        && close(force.x , that.force.x , err)
        && close(force.y , that.force.y , err)
        ;
}

double bead::get_viscosity()
{
    return visc;
}

double bead::get_friction()
{
    return friction;
}
