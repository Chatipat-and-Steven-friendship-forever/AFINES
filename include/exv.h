#ifndef AFINES_EXV_H
#define AFINES_EXV_H

#include "globals.h"

class excluded_volume
{
    public:
        excluded_volume(class box *bc_, double rmax_, double kexv_)
        {
            bc = bc_;
            rmax = rmax_;
            kexv = kexv_;
            pe_exv = 0.0;
        }

        void update_spring_forces(class filament_ensemble *net);
        void update_spring_forces_from_quads(
                class quadrants *quads, class filament_ensemble *net);
        void update_force_between_filaments(
                class filament_ensemble *net, array<int, 2> fl1, array<int, 2> fl2);
        void update_excluded_volume(class filament_ensemble *net);

        double get_pe_exv() { return pe_exv; }
        virial_type get_vir_exv() { return vir_exv; }

    protected:
        class box *bc;
        double rmax, kexv;
        double pe_exv;
        virial_type vir_exv;
};

#endif
