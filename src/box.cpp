#include "box.h"

box::box(string BC_, double xbox_, double ybox_, double delrx_)
{
    BC = string2bc(BC_);
    xbox = xbox_;
    ybox = ybox_;
    delrx = delrx_;
}

bc_type box::get_BC()
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

void box::add_callback(function<void(double)> callback)
{
    callbacks.push_back(callback);
}

void box::update_d_strain(double d_strain)
{
    delrx += d_strain;
    for (auto callback : callbacks) {
        callback(d_strain);
    }
}

// modified from www.cplusplus.com/forum/articles/3638/-
double roundhalfup(double value)
{
    return floor(value + 0.5);
}

// using the minimum image convention
// Allen and Tildesley, page 30 (periodic)
// Allen and Tildesley, page 247 (Lees-Edwards)
vec_type box::rij_bc(vec_type disp)
{
    double dx = disp.x;
    double dy = disp.y;
    if (BC == bc_type::periodic) {
        dx -= xbox * roundhalfup(dx / xbox);
        dy -= ybox * roundhalfup(dy / ybox);

    } else if (BC == bc_type::xperiodic) {
        dx -= xbox * roundhalfup(dx / xbox);

    } else if (BC == bc_type::lees_edwards) {
        double cory = roundhalfup(dy / ybox);
        dx -= delrx * cory;
        dx -= xbox * roundhalfup(dx / xbox);
        dy -= ybox * cory;

    }
    return {dx, dy};
}

vec_type box::pos_bc(vec_type pos)
{
    double x = pos.x;
    double y = pos.y;

    if (BC == bc_type::periodic) {
        x -= xbox * roundhalfup(x / xbox);
        y -= ybox * roundhalfup(y / ybox);

    } else if (BC == bc_type::lees_edwards) {
        double cory = roundhalfup(y / ybox);
        x -= delrx * cory;
        x -= xbox * roundhalfup(x / xbox);
        y -= ybox * cory;

    } else if (BC == bc_type::xperiodic) {
        double xlo = -0.5 * xbox;
        double xhi =  0.5 * xbox;
        double ylo = -0.5 * ybox;
        double yhi =  0.5 * ybox;
        if (x < xlo)
            x += xbox;
        else if (x > xhi)
            x -= xbox;
        if (y < ylo || y > yhi)
            throw runtime_error("Coordinate outside of box.");

    } else if (BC == bc_type::nonperiodic) {
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

double box::dist_bc(vec_type disp)
{
    return abs(rij_bc(disp));
}

double box::dot_bc(vec_type disp1, vec_type disp2)
{
    return dot(rij_bc(disp1), rij_bc(disp2));
}

boost::optional<vec_type> seg_seg_intersection_bc(box *bc, vec_type r1, vec_type r2, vec_type r3, vec_type r4)
{
    vec_type rij12 = bc->rij_bc(r2 - r1);
    vec_type rij13 = bc->rij_bc(r3 - r1);
    vec_type rij34 = bc->rij_bc(r4 - r3);
    vec_type rij14 = rij13 + rij34;

    boost::optional<vec_type> inter = seg_seg_intersection({}, rij12, rij13, rij14);
    if (inter) {
        return bc->pos_bc(*inter + r1);
    } else {
        return boost::none;
    }
}

bc_type box::string2bc(string BC)
{
    if (BC == "LEES-EDWARDS") return bc_type::lees_edwards;
    if (BC == "PERIODIC") return bc_type::periodic;
    if (BC == "XPERIODIC") return bc_type::xperiodic;
    if (BC == "NONPERIODIC") return bc_type::nonperiodic;
    throw std::runtime_error("Boundary condition " + BC + " not recognized.");
}
