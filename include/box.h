#ifndef AFINES_BOX_H
#define AFINES_BOX_H

#include "globals.h"

enum class bc_type {
    lees_edwards,
    periodic,
    xperiodic,
    nonperiodic
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

        vec_type rij_bc(vec_type disp);
        vec_type pos_bc(vec_type pos);
        double dist_bc(vec_type disp);
        double dot_bc(vec_type disp1, vec_type disp2);

    protected:
        bc_type string2bc(string);

        bc_type BC;
        double xbox, ybox, delrx;
        vector<function<void(double)>> callbacks;
};

boost::optional<vec_type> seg_seg_intersection_bc(box *bc, vec_type r1, vec_type r2, vec_type r3, vec_type r4);

#endif
