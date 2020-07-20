#ifndef AFINES_THERMO_H
#define AFINES_THERMO_H

class thermo {
    public:
        double get_stretching_energy();
        double get_bending_energy();
    protected:
        double ke_vel, ke_vir;
        double pe_stretch, pe_bend;
        virial_type vir_stretch, vir_bend;
};

#endif
