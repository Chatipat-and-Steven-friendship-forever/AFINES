#ifndef AFINES_BOX_H
#define AFINES_BOX_H

#include "globals.h"

enum class bc_type {
    clip,
    infinite,
    lees_edwards,
    periodic,
    reflective,
    xperiodic
};

class box {
    public:
        box(string BC_, double xbox_, double ybox_, double delrx_);

        bc_type get_BC();
        array<double, 2> get_fov();
        double get_xbox();
        double get_ybox();
        double get_delrx();

        void add_callback(function<void(double)> callback);
        void update_d_strain(double d_strain);

        array<double, 2> rij_bc(array<double, 2> disp);
        array<double, 2> pos_bc(array<double, 2> pos, array<double, 2> vel, double dt);
        double dist_bc(array<double, 2> disp);
        double dot_bc(array<double, 2> disp1, array<double, 2> disp2);

    protected:
        bc_type string2bc(string);

        bc_type BC;
        double xbox, ybox, delrx;
        vector<function<void(double)>> callbacks;
};

boost::optional<array<double, 2> > seg_seg_intersection_bc(box *bc, array<double, 2> r1, array<double, 2> r2, array<double, 2> r3, array<double, 2> r4);

#endif
