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

// modified from www.cplusplus.com/forum/articles/3638/-
double roundhalfup(double value)
{
    return floor(value + 0.5);
}

// using the minimum image convention
// Allen and Tildesley, page 30 (periodic)
// Allen and Tildesley, page 247 (Lees-Edwards)
array<double, 2> box::rij_bc(array<double, 2> disp)
{
    double dx = disp[0];
    double dy = disp[1];
    if (BC == "PERIODIC") {
        dx -= xbox * roundhalfup(dx / xbox);
        dy -= ybox * roundhalfup(dy / ybox);

    } else if (BC == "XPERIODIC") {
        dx -= xbox * roundhalfup(dx / xbox);

    } else if (BC == "LEES-EDWARDS") {
        double cory = roundhalfup(dy / ybox);
        dx -= delrx * cory;
        dx -= xbox * roundhalfup(dx / xbox);
        dy -= ybox * cory;

    }
    return {dx, dy};
}

array<double, 2> box::pos_bc(array<double, 2> pos, array<double, 2> vel, double dt)
{
    double x = pos[0];
    double y = pos[1];

    if (BC == "PERIODIC") {
        x -= xbox * roundhalfup(x / xbox);
        y -= ybox * roundhalfup(y / ybox);

    } else if (BC == "LEES-EDWARDS") {
        double cory = roundhalfup(y / ybox);
        x -= delrx * cory;
        x -= xbox * roundhalfup(x / xbox);
        y -= ybox * cory;

    } else if (BC == "REFLECTIVE") {
        double xlo = -0.5 * xbox + 2.0 * delrx * y / ybox;
        double xhi =  0.5 * xbox + 2.0 * delrx * y / ybox;
        double ylo = -0.5 * ybox;
        double yhi =  0.5 * ybox;
        if (x <= xlo || x >= xhi)
            x -= 2.0 * dt * vel[0];
        if (y <= ylo || y >= yhi)
            y -= 2.0 * dt * vel[1];

    } else if (BC == "INFINITE") {
        double xlo = -0.5 * xbox + 2.0 * delrx * y / ybox;
        double xhi =  0.5 * xbox + 2.0 * delrx * y / ybox;
        double ylo = -0.5 * ybox;
        double yhi =  0.5 * ybox;
        if (x <= xlo)
            x = xlo;
        else if (x >= xhi)
            x = xhi;
        if (y <= ylo)
            y = ylo;
        else if (y >= yhi)
            y = yhi;

    } else if (BC == "XPERIODIC") {
        double xlo = -0.5 * xbox;
        double xhi =  0.5 * xbox;
        double ylo = -0.5 * ybox;
        double yhi =  0.5 * ybox;
        if (x < xlo)
            x += xbox;
        else if (x > xhi)
            x -= xbox;
        if (y < ylo)
            y = ylo;
        else if (y > yhi)
            y = yhi;

    } else if (BC == "CLIP") {
        double xlo = -0.5 * xbox;
        double xhi =  0.5 * xbox;
        double ylo = -0.5 * ybox;
        double yhi =  0.5 * ybox;
        if (x < xlo || x > xhi || y < ylo || y > yhi)
            throw runtime_error("Coordinate outside of box.");

    } else {
        throw runtime_error("Boundary condition not recognized.");

    }

    return {x, y};
}

double box::dist_bc(array<double, 2> disp)
{
    array<double, 2> rij = rij_bc(disp);
    return hypot(rij[0], rij[1]);
}

double box::dot_bc(array<double, 2> disp1, array<double, 2> disp2)
{
    array<double, 2> rij1 = rij_bc(disp1);
    array<double, 2> rij2 = rij_bc(disp2);
    return rij1[0] * rij2[0] + rij1[1] * rij2[1];
}

boost::optional<array<double, 2> > seg_seg_intersection_bc(box *bc, array<double, 2> r1, array<double, 2> r2, array<double, 2> r3, array<double, 2> r4)
{
    array<double, 2> rij12 = bc->rij_bc({r2[0] - r1[0], r2[1] - r1[1]});
    array<double, 2> rij13 = bc->rij_bc({r3[0] - r1[0], r3[1] - r1[1]});
    array<double, 2> rij34 = bc->rij_bc({r4[0] - r3[0], r4[1] - r3[1]});
    array<double, 2> rij14 = {rij13[0] + rij34[0], rij13[1] + rij34[1]};

    boost::optional<array<double, 2>> inter = seg_seg_intersection({0.0, 0.0}, rij12, rij13, rij14);
    if (inter) {
        return bc->pos_bc({inter->at(0) + r1[0], inter->at(1) + r1[1]}, {0.0, 0.0}, 0.0);
    } else {
        return boost::none;
    }
}
