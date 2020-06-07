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

filament::filament(filament_ensemble *net, vector<bead *> beadvec, double spring_length,
        double stretching_stiffness, double max_ext_ratio, double bending_stiffness,
        double deltat, double temp, double frac_force)
{
    filament_network = net;
    
    bc = net->get_box();
    dt = deltat;
    temperature = temp;
    fracture_force = frac_force;
    kb = bending_stiffness;
    y_thresh = 1;
    kinetic_energy = 0;

    if (beadvec.size() > 0)
    {
        beads.push_back(new bead(*(beadvec[0])));
        prv_rnds.push_back({{0,0}});
        damp = beads[0]->get_friction();
    }
    
    //spring em up
    if (beadvec.size() > 1){
        for (unsigned int j = 1; j < beadvec.size(); j++) {

            beads.push_back(new bead(*(beadvec[j])));
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

void filament::add_bead(bead * a, double spring_length, double stretching_stiffness, double max_ext_ratio){
    
    beads.push_back(new bead(*a));
    prv_rnds.push_back({{0,0}});    
    if (beads.size() > 1){
        int j = (int) beads.size() - 1;
        springs.push_back(new spring(spring_length, stretching_stiffness, max_ext_ratio, this, {j-1,  j}));
        springs[j-1]->step();
    }
    if (damp == infty)
        damp = a->get_friction();
}

void filament::set_y_thresh(double y){
    y_thresh = y;
}

void filament::update_positions()
{
    double vx, vy;
    array<double, 2> new_rnds;
    array<double, 2> newpos;
    kinetic_energy = 0;  
    double top_y = y_thresh*bc->get_ybox()/2.;
    int sa = int(beads.size());
    int la = int(springs.size());
    for (int i = 0; i < sa; i++){

        if (fabs(beads[i]->get_ycm()) > top_y) continue;
     
        new_rnds = {{rng_n(), rng_n()}};
        vx  = (beads[i]->get_force()[0])/damp  + bd_prefactor*(new_rnds[0] + prv_rnds[i][0]);
        vy  = (beads[i]->get_force()[1])/damp  + bd_prefactor*(new_rnds[1] + prv_rnds[i][1]);
//        cout<<"\nDEBUG: Fx("<<i<<") = "<<beads[i]->get_force()[0]<<"; v = ("<<vx<<" , "<<vy<<")";
       
        prv_rnds[i] = new_rnds;
        //cout<<"\nDEBUG: bead force = ("<<beads[i]->get_force()[0]<<" , "<<beads[i]->get_force()[1]<<")";
        kinetic_energy += vx*vx + vy*vy;
        newpos = bc->pos_bc({beads[i]->get_xcm() + vx*dt, beads[i]->get_ycm() + vy*dt}, {vx, vy}, dt);
        beads[i]->set_xcm(newpos[0]);
        beads[i]->set_ycm(newpos[1]);
        beads[i]->reset_force(); 
    }

    for (int i = 0; i < la; i++)
        springs[i]->step();

}

void filament::update_positions_range(int lo, int hi)
{
    double vx, vy;
    array<double, 2> new_rnds;
    array<double, 2> newpos;
    kinetic_energy = 0;  
    double top_y = y_thresh*bc->get_ybox()/2.;

    int low = max(0, lo);
    int high = min(hi, (int)beads.size());

    for (int i = low; i < high; i++){
       
        if (fabs(beads[i]->get_ycm()) > top_y) continue;
     
        new_rnds = {{rng_n(), rng_n()}};
        vx  = (beads[i]->get_force()[0])/damp  + bd_prefactor*(new_rnds[0] + prv_rnds[i][0]);
        vy  = (beads[i]->get_force()[1])/damp  + bd_prefactor*(new_rnds[1] + prv_rnds[i][1]);
//        cout<<"\nDEBUG: Fx("<<i<<") = "<<beads[i]->get_force()[0]<<"; v = ("<<vx<<" , "<<vy<<")";
       
        prv_rnds[i] = new_rnds;
        //cout<<"\nDEBUG: bead force = ("<<beads[i]->get_force()[0]<<" , "<<beads[i]->get_force()[1]<<")";
        kinetic_energy += vx*vx + vy*vy;
        newpos = bc->pos_bc({beads[i]->get_xcm() + vx*dt, beads[i]->get_ycm() + vy*dt}, {vx, vy}, dt);
        beads[i]->set_xcm(newpos[0]);
        beads[i]->set_ycm(newpos[1]);
        beads[i]->reset_force(); 
    }

    for (unsigned int i = 0; i < springs.size(); i++)
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
        //springs[i]->update_force_fraenkel_fene();
        if (spring_force[0]*spring_force[0]+spring_force[1]*spring_force[1] > fracture_force_sq){
//        if ((springs[i]->get_force()[0], springs[i]->get_force()[1]) > fracture_force){
            newfilaments = this->fracture(i);
            break;
        }
        else 
            springs[i]->filament_update();
    }
    
    return newfilaments;
}

bead * filament::get_bead(int i)
{
    try
    {
        return beads[i];
    }
    catch (int e)
    {
        throw std::runtime_error("An exception occurred while returning bead " + std::to_string(i) + ".");
    }
}

spring * filament::get_spring(int i)
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
    return {get_kinetic_energy(), get_potential_energy(), get_total_energy(), double(fil)};
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
    return "\n" + std::to_string(this->get_kinetic_energy()) + \
        "\t" + std::to_string(this->get_potential_energy()) + \
        "\t" + std::to_string(this->get_total_energy()) + "\t" + std::to_string(fil);
}

vector<bead *> filament::get_beads(unsigned int first, unsigned int last)
{
    vector<bead *> newbeads;
    for (unsigned int i = first; i < last; i++)
    {
        if (i >= beads.size())
            break;
        else
            newbeads.push_back(new bead(*(beads[i])));
    }
    return newbeads;
}

vector<filament *> filament::fracture(int node){

    vector<filament *> newfilaments;
    cout<<"\n\tDEBUG: fracturing at node "<<node;

    if(springs.size() == 0)
        return newfilaments;

    vector<bead *> lower_half = this->get_beads(0, node+1);
    vector<bead *> upper_half = this->get_beads(node+1, beads.size());

    if (lower_half.size() > 0)
        newfilaments.push_back(
                new filament(filament_network, lower_half, springs[0]->get_l0(), springs[0]->get_kl(), springs[0]->get_fene_ext(), kb,
                    dt, temperature, fracture_force));
    if (upper_half.size() > 0)
        newfilaments.push_back(
                new filament(filament_network, upper_half, springs[0]->get_l0(), springs[0]->get_kl(), springs[0]->get_fene_ext(), kb,
                    dt, temperature, fracture_force));

    for (int i = 0; i < (int)(lower_half.size()); i++) delete lower_half[i];
    for (int i = 0; i < (int)(upper_half.size()); i++) delete upper_half[i];
    
    lower_half.clear();
    upper_half.clear();
    
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
    string out = "";

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

/* --------------------------------------------------------------------------------- */
/* copied and pasted and edited the following code from LAMMPS src/angle_harmonic.cpp*/
/* --------------------------------------------------------------------------------- */
void filament::lammps_bending_update()
{
    array<double, 2> delr1, delr2;
    double f1[2], f3[2];
    double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;

    virial_clear(bending_virial);

    // 1st bond
    delr1 = springs[0]->get_neg_disp();
    rsq1  = springs[0]->get_length_sq();
    r1    = springs[0]->get_length();

    double theta = 0, totThetaSq = 0;

    for (int n = 0; n < int(springs.size())-1; n++) {

        // 2nd bond
        delr2 = springs[n+1]->get_disp();
        rsq2  = springs[n+1]->get_length_sq();
        r2    = springs[n+1]->get_length();

        // angle (cos and sin)
        c = (delr1[0]*delr2[0] + delr1[1]*delr2[1]) / (r1*r2);

        if (c > 1.0) c = 1.0;
        if (c < -1.0) c = -1.0;

        s = sqrt(1.0 - c*c);
        if (s < maxSmallAngle) s = maxSmallAngle;

        theta = acos(c) - pi;
        totThetaSq += theta*theta;

        // force
        a   = -kb * theta / s; //Note, in this implementation, Lp = kb/kT 
        a11 = a*c / rsq1;
        a12 = -a / (r1*r2);
        a22 = a*c / rsq2;

        f1[0] = a11*delr1[0] + a12*delr2[0];
        f1[1] = a11*delr1[1] + a12*delr2[1];
        f3[0] = a22*delr2[0] + a12*delr1[0];
        f3[1] = a22*delr2[1] + a12*delr1[1];
        //cout<<"\nDEBUG: f1x, f1y, f3x f3y = "<<f1[0]<<" , "<<f1[1]<<" , "<<f3[0]<<" , "<<f3[1];

        // apply force to each of 3 atoms
        beads[n  ]->update_force(f1[0], f1[1]);
        beads[n+1]->update_force(-f1[0] - f3[0], -f1[1] - f3[1]);
        beads[n+2]->update_force(f3[0], f3[1]);

        // 1st bond, next iteration
        delr1 = {{-delr2[0], -delr2[1]}};
        rsq1  = rsq2;
        r1    = r2;

        bending_virial[0][0] += a11 * delr1[0] * delr1[0] + 2.0 * a12 * delr1[0] * delr2[0] + a22 * delr2[0] * delr2[0];
        bending_virial[0][1] += a11 * delr1[0] * delr1[1] + a12 * (delr1[0] * delr2[1] + delr1[1] * delr2[0]) + a22 * delr2[0] * delr2[1];
        bending_virial[1][0] += a11 * delr1[0] * delr1[1] + a12 * (delr1[0] * delr2[1] + delr1[1] * delr2[0]) + a22 * delr2[0] * delr2[1];
        bending_virial[1][1] += a11 * delr1[1] * delr1[1] + 2.0 * a12 * delr1[1] * delr2[1] + a22 * delr2[1] * delr2[1];
    }

    ubend = kb*totThetaSq/2;
}

void filament::update_bending(double t)
{
    if(springs.size() > 1 && kb > 0){
          this->lammps_bending_update();
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

void filament::init_ubend(){
   

    if (springs.size() < 2) 
        ubend = 0;
    else{
        double sum = 0, theta;

        for (unsigned int i = 0; i < springs.size() - 1; i++)
        {
            theta = angle_between_springs(i+1, i);
            sum += theta*theta;
        }
        
        ubend = kb*sum/2;
    }
}

double filament::get_stretching_energy()
{
    
    double u = 0;
    for (unsigned int i = 0; i < springs.size(); i++)
        //u += springs[i]->get_stretching_energy_fene();
        u += springs[i]->get_stretching_energy();
    
    return u;

}

double filament::get_kinetic_energy()
{
    return kinetic_energy;
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
    return this->get_potential_energy() + this->get_kinetic_energy();
}

array<double, 2> filament::get_bead_position(int n)
{
    return {{beads[n]->get_xcm(), beads[n]->get_ycm()}};
}

void filament::print_thermo()
{
    cout<<"\tKE = "<<this->get_kinetic_energy()<<"\tPE = "<<this->get_potential_energy()<<\
        "\tTE = "<<this->get_total_energy();
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
