#include "quadrants.h"

quadrants::quadrants(box *bc_, array<int, 2> nq_)
{
    bc = bc_;
    nq = nq_;
    quad_flag = true;

    quads = new vector<array<int, 2>> *[nq[0]];
    quads[0] = new vector<array<int, 2>>[nq[0] * nq[1]];
    for (int x = 0; x < nq[0]; x++)
        quads[x] = quads[0] + x * nq[1];
}

quadrants::~quadrants()
{
    delete[] quads[0];
    delete quads;
}

void quadrants::use_quad(bool flag)
{
    quad_flag = flag;
}

void quadrants::add_spring(spring *s, array<int, 2> fl)
{
    if (!quad_flag)
        all_springs.push_back(fl);
    if (bc->get_BC() == bc_type::periodic || bc->get_BC() == bc_type::lees_edwards)
        add_spring_periodic(s, fl);
    else
        add_spring_nonperiodic(s, fl);
}

vector<array<int, 2>> *quadrants::get_attach_list(array<double, 2> pos)
{
    if (!quad_flag) {
        return &all_springs;
    } else {
        double x = pos[0];
        double y = pos[1];
        array<double, 2> fov = bc->get_fov();
        double delrx = bc->get_delrx();
        int iy = round(nq[1] * (y / fov[1] + 0.5));
        while (iy < 0) {
            iy += nq[1];
            x += delrx;
        }
        while (iy >= nq[1]) {
            iy -= nq[1];
            x -= delrx;
        }
        int ix = round(nq[0] * (x / fov[0] + 0.5));
        while (ix < 0) ix += nq[0];
        while (ix >= nq[0]) ix -= nq[0];
        if (!(0 <= ix && ix < nq[0] && 0 <= iy && iy < nq[1])) {
            throw std::logic_error("Invalid quadrant index.");
        }
        return &quads[ix][iy];
    }
}

void quadrants::clear()
{
    all_springs.clear();
    for (int x = 0; x < nq[0]; x++)
        for (int y = 0; y < nq[1]; y++)
            quads[x][y].clear();
}

void quadrants::check_duplicates()
{
    for (int x = 0; x < nq[0]; x++) {
        for (int y = 0; y < nq[1]; y++) {
            vector<array<int, 2>> &quad = quads[x][y];

            set<array<int, 2>> s(quad.begin(), quad.end());
            if (s.size() != quad.size())
                cerr << "Quadrant (" << x << ", " << y << ") contains duplicates!" << endl;

            /*
            set<array<int, 2>> s;
            for (array<int, 2> fl : quad) {
                if (s.find(fl) == s.end())
                    s.insert(fl);
                else
                    cerr << "Quadrant (" << x << ", " << y << ") contains duplicates!" << endl;
            }
            */
        }
    }
}

void quadrants::add_spring_nonperiodic(spring *s, array<int, 2> fl)
{
    array<double, 2> fov = bc->get_fov();
    array<double, 2> hx = s->get_hx();
    array<double, 2> hy = s->get_hy();
    array<double, 2> disp = s->get_disp();

    double xlo = hx[0], xhi = hx[1];
    if (disp[0] < 0) std::swap(xlo, xhi);

    double ylo = hy[0], yhi = hy[1];
    if (disp[1] < 0) std::swap(ylo, yhi);

    int xlower = floor(nq[0] * (xlo / fov[0] + 0.5));
    int xupper = ceil(nq[0] * (xhi / fov[0] + 0.5));
    if (xlower < 0) {
        cout << "Warning: x-index of quadrant < 0." << endl;
        xlower = 0;
    }
    if (xupper > nq[0]) {
        cout << "Warning: x-index of quadrant > nq[0]." << endl;
        xupper = nq[0];
    }
    if (xlower > xupper) throw std::logic_error("xlower > xupper");

    int ylower = floor(nq[1] * (ylo / fov[1] + 0.5));
    int yupper = ceil(nq[1] * (yhi / fov[1] + 0.5));
    if (ylower < 0) {
        cout << "Warning: y-index of quadrant < 0." << endl;
        ylower = 0;
    }
    if (yupper > nq[1]) {
        cout << "Warning: y-index of quadrant > nq[1]." << endl;
        yupper = nq[1];
    }
    if (ylower > yupper) throw std::logic_error("ylower > yupper");

    for (int i = xlower; i <= xupper; i++)
        for (int j = ylower; j <= yupper; j++)
            quads[i][j].push_back(fl);
}

void quadrants::add_spring_periodic(spring *s, array<int, 2> fl)
{
    array<double, 2> fov = bc->get_fov();
    double delrx = bc->get_delrx();
    array<double, 2> hx = s->get_hx();
    array<double, 2> hy = s->get_hy();
    array<double, 2> disp = s->get_disp();

    double xlo, xhi;
    double ylo, yhi;
    if (disp[1] >= 0) {
        ylo = hy[0];
        yhi = hy[0] + disp[1];
        if (disp[0] >= 0) {
            xlo = hx[0];
            xhi = hx[0] + disp[0];
        } else {
            xlo = hx[0] + disp[0];
            xhi = hx[0];
        }
    } else {
        ylo = hy[1];
        yhi = hy[1] - disp[1];
        if (disp[0] >= 0) {
            xlo = hx[1] - disp[0];
            xhi = hx[1];
        } else {
            xlo = hx[1];
            xhi = hx[1] - disp[0];
        }
    }
    if (xlo > xhi) throw std::logic_error("xlo > xhi");
    if (ylo > yhi) throw std::logic_error("ylo > yhi");

    int ylower = floor(nq[1] * (ylo / fov[1] + 0.5));
    int yupper =  ceil(nq[1] * (yhi / fov[1] + 0.5));
    if (ylower > yupper) throw std::logic_error("ylower > yupper");

    for (int jj = ylower; jj <= yupper; jj++) {
        int j = jj;

        double xlo_new = xlo;
        double xhi_new = xhi;
        while (j < 0) {
            j += nq[1];
            xlo_new += delrx;
            xhi_new += delrx;
        }
        while (j >= nq[1]) {
            j -= nq[1];
            xlo_new -= delrx;
            xhi_new -= delrx;
        }
        if (!(0 <= j && j < nq[1])) throw std::logic_error("y quadrant index out of bounds");

        int xlower = floor(nq[0] * (xlo_new / fov[0] + 0.5));
        int xupper =  ceil(nq[0] * (xhi_new / fov[0] + 0.5));
        if (xlower > xupper) throw std::logic_error("xlower > xupper");

        for (int ii = xlower; ii <= xupper; ii++) {
            int i = ii;

            while (i < 0) i += nq[0];
            while (i >= nq[0]) i -= nq[0];
            if (!(0 <= i && i < nq[0])) throw std::logic_error("x quadrant index out of bounds");

            quads[i][j].push_back(fl);
        }
    }
}
