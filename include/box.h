#ifndef __BOX_H_INCLUDED__
#define __BOX_H_INCLUDED__

#include "globals.h"

class box {
    public:
        box(string BC_, double xbox_, double ybox_, double delrx_);

        string get_BC();
        array<double, 2> get_fov();
        double get_xbox();
        double get_ybox();
        double get_delrx();

        void update_d_strain(double d_strain);

        array<double, 2> rij_bc(array<double, 2> disp);
        array<double, 2> pos_bc(array<double, 2> pos, array<double, 2> vel, double dt);
        double dist_bc(array<double, 2> disp);
        double dot_bc(array<double, 2> disp1, array<double, 2> disp2);

    protected:
        string BC;
        double xbox, ybox, delrx;
};

boost::optional<array<double, 2> > seg_seg_intersection_bc(box *bc, array<double, 2> r1, array<double, 2> r2, array<double, 2> r3, array<double, 2> r4);

#endif
