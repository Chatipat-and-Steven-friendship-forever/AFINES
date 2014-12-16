/*
 * motor.h
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __MOTOR_H_INCLUDED__
#define __MOTOR_H_INCLUDED__

//=====================================
// forward declared dependencies
class actin_ensemble;

//=====================================
//included dependences
#include "map"
#include "string"

//motor class
class motor 
{
    public:

        motor(double mx, double my, double mang, double mlen, actin_ensemble* network, int state0, int state1, int
                aindex0, int aindex1, double fovx, double fovy, double delta_t, double v0, double temp, double stiffness, double ron, double roff,
                double rend, double actin_len, double vis, std::string col);

        ~motor();

        int* get_states();

        double* get_hx();

        double* get_hy();

        std::string get_color();

        double tension();

        void attach(int hd);

        void brownian(double t, double gamma);

        void step_onehead(int hd);

        void step_twoheads();

        void actin_update();

        void update_shape();

    private:

        double hx[2],hy[2], xm[2], ym[2], mphi,mld,mobility,fm[2],vm[2],offrate[2], onrate, stretch,forcex[2],forcey[2],torque[2], force_par[2],force_perp[2],vs, pos_temp;

        int state[2], aindex[2];
        
        double dm,fmax,mk, kon, koff, kend, dt, temperature;
        
        std::map<int, double> dist;
        
        std::string color;
        
        actin_ensemble *actin_network;
        
        double pos_a_end[2], fov[2];

        inline void move_end_detach(int hd, double pos);

        inline void reflect(double t, double gamma, double x1, double x2, double y1, double y2);
};

#endif
