#include "box.h"

box::box(string BC_, double xbox_, double ybox_, double delrx_)
{
    BC = BC_;
    xbox = xbox_;
    ybox = ybox_;
    delrx = delrx_;
}

string box::get_BC()
{
    return BC;
}

array<double, 2> box::get_fov()
{
    return {xbox, ybox};
}

double box::get_xbox()
{
    return xbox;
}

double box::get_ybox()
{
    return ybox;
}

double box::get_delrx()
{
    return delrx;
}

void box::update_d_strain(double d_strain)
{
    delrx += d_strain;
}

array<double, 2> box::rij_bc(array<double, 2> disp)
{
    return ::rij_bc(BC, disp[0], disp[1], xbox, ybox, delrx);
}

array<double, 2> box::pos_bc(array<double, 2> pos, array<double, 2> vel, double dt)
{
    return ::pos_bc(BC, delrx, dt, {xbox, ybox}, vel, pos);
}

double box::dist_bc(array<double, 2> disp)
{
    return ::dist_bc(BC, disp[0], disp[1], xbox, ybox, delrx);
}

double box::dot_bc(array<double, 2> disp1, array<double, 2> disp2)
{
    return ::dot_bc(BC, disp1[0], disp1[1], disp2[0], disp2[1], xbox, ybox, delrx);
}
