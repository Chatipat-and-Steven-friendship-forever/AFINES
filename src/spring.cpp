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
    
    hx = {{0,0}};
    hy = {{0,0}};

    force = {{0,0}};
    intpoint = {{0,0}};
    llensq = l0*l0;
    llen = l0;
}
spring::~spring(){ 
};

array<double,2> spring::get_hx(){
    return hx;
}

array<double,2> spring::get_hy(){
    return hy;
}

// stepping kinetics

void spring::step()
{
    hx[0] = fil->get_bead(aindex[0])->get_xcm();
    hx[1] = fil->get_bead(aindex[1])->get_xcm();
    hy[0] = fil->get_bead(aindex[0])->get_ycm();
    hy[1] = fil->get_bead(aindex[1])->get_ycm();

    disp   = bc->rij_bc({hx[1]-hx[0], hy[1]-hy[0]});
    llensq = disp[0]*disp[0] + disp[1]*disp[1];
    llen   = sqrt(llensq);

    if (llen != 0)
        direc = {{disp[0]/llen, disp[1]/llen}};
    else
        direc = {{0, 0}};

}

void spring::update_force()
{
    double kf = kl*(llen-l0);
    force = {{kf*direc[0], kf*direc[1]}};
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
    force = {{klp*direc[0], klp*direc[1]}};

}

void spring::update_force_marko_siggia(double kToverLp)
{
    double xrat = llen/l0, yrat = llen/l0;
    if (xrat != xrat || xrat == 1) xrat = 0;
    if (yrat != yrat || yrat == 1) yrat = 0;
    force = {{kToverLp*(0.25/((1-xrat)*(1-xrat))-0.25+xrat), kToverLp*(0.25/((1-yrat)*(1-yrat))-0.25+yrat)}};  
}

array<double,2> spring::get_force()
{
    return force;
}

array<double,2> spring::get_disp()
{
    return disp;
}

array<double,2> spring::get_neg_disp()
{
    return {{-disp[0], -disp[1]}};
}

void spring::filament_update()
{
    fil->update_forces(aindex[0],  force[0],  force[1]);
    fil->update_forces(aindex[1], -force[0], -force[1]);

}

double spring::get_kl(){
    return kl;
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
    return {hx[0], hy[0], disp[0], disp[1]};
}

std::string spring::write()
{
    return "\n" + std::to_string(hx[0]) + "\t" + std::to_string(hy[0]) + "\t" + std::to_string(disp[0]) + "\t" 
        + std::to_string(disp[1]);
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
double spring::get_distance_sq(double xp, double yp)
{
    array<double, 2> dr = bc->rij_bc({intpoint[0]-xp, intpoint[1]-yp});
    return dr[0]*dr[0] + dr[1]*dr[1];
}

array<double,2> spring::get_intpoint()
{
    return intpoint;
}

void spring::calc_intpoint(double xp, double yp)
{
    double l2 = disp[0]*disp[0]+disp[1]*disp[1];
    if (l2==0){
        intpoint = {{hx[0], hy[0]}};
    }else{
        //Consider the line extending the spring, parameterized as h0 + tp ( h1 - h0 )
        //tp = projection of (xp, yp) onto the line
        double tp = bc->dot_bc({xp-hx[0], yp-hy[0]}, {hx[1]-hx[0], hy[1]-hy[0]})/l2;
        if (tp<0){ 
            intpoint = {{hx[0], hy[0]}};
        }
        else if(tp>1.0){
            intpoint = {{hx[1], hy[1]}};
        }
        else{
            array<double, 2> proj   = {{hx[0] + tp*disp[0], hy[0] + tp*disp[1]}};
            intpoint                = bc->pos_bc(proj, {0.0, 0.0}, 0.0); //velocity and dt are 0 since not relevant
        }
    }
}

array<double, 2> spring::get_direction()
{
    return direc;
}

double spring::get_stretching_energy(){
    return (force[0]*force[0]+force[1]*force[1])/(2*kl);
}

array<array<double, 2>, 2> spring::get_virial() {
    double k = kl*(llen-l0)/llen;
    return {
        array<double, 2>{k * disp[0] * disp[0], k * disp[0] * disp[1]},
        array<double, 2>{k * disp[0] * disp[1], k * disp[1] * disp[1]}
    };
}

double spring::get_stretching_energy_fene()
{
    double ext = abs(l0 - llen);
    
    if (max_ext - ext > eps_ext )
        return -0.5*kl*max_ext*max_ext*log(1-(ext/max_ext)*(ext/max_ext));
    else
        return 0.25*kl*ext*ext*(max_ext/eps_ext);
    
}
