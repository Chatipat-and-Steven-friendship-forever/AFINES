#ifndef AFINES_EXV_H
#define AFINES_EXV_H

#include "filament.h"
#include "quadrants.h"
#include "box.h"

class excluded_volume
{
    public:
        excluded_volume(box *bc_, double rmax_, double kexv_)
        {
            bc = bc_;
            rmax = rmax_;
            kexv = kexv_;
        }

        void update_spring_forces(vector<filament *> &network, int f);
        void update_spring_forces_from_quads(quadrants *quads, vector<filament *> &network);
        void update_force_between_filaments(
                vector<filament *> &network, int n1, int l1, int n2, int l2);
        void update_excluded_volume(vector<filament *> &network, int f);

        double get_pe_exv() { return pe_exv; }

    protected:
        box *bc;
        double rmax, kexv;
        double pe_exv;
};

#endif
