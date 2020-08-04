#include "quadrants.h"

#include "box.h"
#include "filament_ensemble.h"

quadrants::quadrants(box *bc_, array<int, 2> nq_)
{
    bc = bc_;
    nq = nq_;
    all_flag = false;  // build all_springs
    quad_flag = true;  // build quadrants

    // quadrants at the box edges of a periodic boundary are the same
    int nx, ny;
    if (bc->get_BC() == bc_type::nonperiodic) {
        nx = nq[0] + 1;
        ny = nq[1] + 1;
    } else if (bc->get_BC() == bc_type::xperiodic) {
        nx = nq[0];
        ny = nq[1] + 1;
    } else {
        nx = nq[0];
        ny = nq[1];
    }

    quads = new vector<array<int, 2>> *[nx];
    quads[0] = new vector<array<int, 2>>[nx * ny];
    for (int x = 0; x < nx; x++)
        quads[x] = quads[0] + x * ny;
}

quadrants::~quadrants()
{
    delete[] quads[0];
    delete[] quads;
}

void quadrants::use_quad(bool flag)
{
    quad_flag = flag;
}

void quadrants::use_all(bool flag)
{
    all_flag = flag;
}

void quadrants::add_spring(vec_type h0, vec_type disp, array<int, 2> fl)
{
    if (all_flag) {
        all_springs.push_back(fl);
    }
    if (quad_flag) {

        // compute bounding box of spring

        double xlo = h0.x;
        double xhi = h0.x + disp.x;
        if (disp.x < 0) std::swap(xlo, xhi);

        double ylo = h0.y;
        double yhi = h0.y + disp.y;
        if (disp.y < 0) std::swap(ylo, yhi);

        assert(xlo <= xhi);
        assert(ylo <= yhi);

        if (bc->get_BC() == bc_type::nonperiodic) {
            this->add_spring_nonperiodic(xlo, xhi, ylo, yhi, fl);
        } else if (bc->get_BC() == bc_type::xperiodic) {
            this->add_spring_xperiodic(xlo, xhi, ylo, yhi, fl);
        } else if (bc->get_BC() == bc_type::periodic) {
            this->add_spring_periodic(xlo, xhi, ylo, yhi, fl);
        } else if (bc->get_BC() == bc_type::lees_edwards) {
            this->add_spring_lees_edwards(xlo, xhi, ylo, yhi, fl);
        } else {
            throw std::logic_error("Unknown boundary condition.");
        }
    }
}

void quadrants::build_pairs()
{
    // TODO: build a Verlet list on top of this
    pairs.clear();

    if (quad_flag) {
        for (int ix = 0; ix < nq[0]; ix++) {
            for (int iy = 0; iy < nq[1]; iy++) {
                vector<array<int, 2>> &q = quads[ix][iy];
                for (size_t i = 0; i < q.size(); i++) {
                    for (size_t j = i + 1; j < q.size(); j++) {
                        // insertion order is preserved in each quadrant,
                        // so we don't risk double counting
                        // TODO: add check for double counting, just in case
                        pairset.insert({q[i], q[j]});
                    }
                }
            }
        }
        pairs.assign(begin(pairset), end(pairset));
        pairset.clear();

    } else if (all_flag) {
        for (size_t i = 0; i < all_springs.size(); i++) {
            for (size_t j = i + 1; i < all_springs.size(); i++) {
                pairs.push_back({all_springs[i], all_springs[j]});
            }
        }

    } else {
        throw std::runtime_error("Neighbor list disabled.");
    }
}

vector<array<array<int, 2>, 2>> *quadrants::get_pairs()
{
    return &pairs;
}

vector<array<int, 2>> *quadrants::get_attach_list(vec_type pos)
{
    if (quad_flag) {
        double x = pos.x;
        double y = pos.y;
        array<double, 2> fov = bc->get_fov();

        int ix, iy;
        if (bc->get_BC() == bc_type::nonperiodic) {
            ix = round(nq[0] * (x / fov[0] + 0.5));
            iy = round(nq[1] * (y / fov[1] + 0.5));

            assert(0 <= ix && ix <= nq[0]);
            assert(0 <= iy && iy <= nq[1]);

        } else if (bc->get_BC() == bc_type::xperiodic) {
            ix = round(nq[0] * (x / fov[0] + 0.5));
            iy = round(nq[1] * (y / fov[1] + 0.5));
            ix = imod(ix, nq[0]);

            assert(0 <= ix && ix < nq[0]);
            assert(0 <= iy && iy <= nq[1]);

        } else if (bc->get_BC() == bc_type::periodic) {
            ix = round(nq[0] * (x / fov[0] + 0.5));
            iy = round(nq[1] * (y / fov[1] + 0.5));
            ix = imod(ix, nq[0]);
            iy = imod(iy, nq[1]);

            assert(0 <= ix && ix < nq[0]);
            assert(0 <= iy && iy < nq[1]);

        } if (bc->get_BC() == bc_type::lees_edwards) {
            double delrx = bc->get_delrx();
            iy = round(nq[1] * (y / fov[1] + 0.5));
            while (iy < 0) {
                iy += nq[1];
                x += delrx;
            }
            while (iy >= nq[1]) {
                iy -= nq[1];
                x -= delrx;
            }
            ix = round(nq[0] * (x / fov[0] + 0.5));
            ix = imod(ix, nq[0]);

            assert(0 <= ix && ix < nq[0]);
            assert(0 <= iy && iy < nq[1]);

        }

        vector<array<int, 2>> *q = &quads[ix][iy];
        return q;

    } else if (all_flag) {
        return &all_springs;

    } else {
        throw std::runtime_error("Neighbor list disabled.");
    }
}

void quadrants::clear()
{
    all_springs.clear();
    for (int x = 0; x < nq[0]; x++)
        for (int y = 0; y < nq[1]; y++)
            quads[x][y].clear();
}

// maximum cutoff supported by quadrants
double quadrants::get_cut()
{
    double xcut = 0.5 * bc->get_xbox() / nq[0];
    double ycut = 0.5 * bc->get_ybox() / nq[1];
    double cut = std::min(xcut, ycut);
    return 0.999 * cut;  // account for numerical error
}

// maximum pair list cutoff supported by quadrants
double quadrants::get_paircut()
{
    double xcut = 2.0 * bc->get_xbox() / nq[0];
    double ycut = 2.0 * bc->get_ybox() / nq[1];
    double cut = std::min(xcut, ycut);
    return 0.999 * cut;  // account for numerical error
}

// begin [checks]

void quadrants::check_duplicates()
{
    // quadrants don't contain duplicate springs
    for (int x = 0; x < nq[0]; x++) {
        for (int y = 0; y < nq[1]; y++) {
            vector<array<int, 2>> &quad = quads[x][y];

            set<array<int, 2>> s(quad.begin(), quad.end());
            assert(s.size() == quad.size());
        }
    }

    // pair lists don't contain duplicate pairs
    set<array<array<int, 2>, 2>> ps(pairs.begin(), pairs.end());
    assert(ps.size() == pairs.size());
    // if (fl0, fl1) is in pair list, then (fl1, fl0) isn't
    for (array<array<int, 2>, 2> pair : pairs) {
        std::swap(pair[0], pair[1]);
        assert(ps.find(pair) == ps.end());
    }
}

// quadrants should contain all springs within a cutoff
void quadrants::check_quad(filament_ensemble *net, vector<array<int, 2>> *q, vec_type pos)
{
    double cut = this->get_cut();

    // all springs within cutoff
    std::vector<array<int, 2>> ref_fl;
    for (array<int, 2> fl : all_springs) {
        vec_type intpoint = net->intpoint(fl[0], fl[1], pos);
        double dist = bc->dist_bc(intpoint - pos);
        if (dist < cut) ref_fl.push_back(fl);
    }

    // springs in quadrant within cutoff
    std::vector<array<int, 2>> q_fl;
    for (array<int, 2> fl : *q) {
        vec_type intpoint = net->intpoint(fl[0], fl[1], pos);
        double dist = bc->dist_bc(intpoint - pos);
        if (dist < cut) q_fl.push_back(fl);
    }

    // check that springs in both lists are the same
    // we can do this since the spring addition order
    // is preserved in each quadrant
    assert(ref_fl.size() == q_fl.size());
    for (size_t i = 0; i < ref_fl.size(); i++) {
        assert(q_fl[i] == ref_fl[i]);
    }
}

// quadrant pair lists should contain all spring-spring pairs within a cutoff
void quadrants::check_pairs(filament_ensemble *net)
{
    double cut = this->get_paircut();

    // all pairs within cutoff
    std::vector<array<array<int, 2>, 2>> ref;
    for (int i = 0; i < int(all_springs.size()); i++) {
        array<int, 2> fli = all_springs[i];

        for (int j = i + 1; j < int(all_springs.size()); j++) {
            array<int, 2> flj = all_springs[j];

            // minimum distance between springs
            vec_type p1 = net->get_start(flj[0], flj[1]);
            vec_type p2 = net->get_end(flj[0], flj[1]);
            vec_type p3 = net->get_start(fli[0], fli[1]);
            vec_type p4 = net->get_end(fli[0], fli[1]);
            double dist1 = bc->dist_bc(net->intpoint(fli[0], fli[1], p1) - p1);
            double dist2 = bc->dist_bc(net->intpoint(fli[0], fli[1], p2) - p2);
            double dist3 = bc->dist_bc(net->intpoint(flj[0], flj[1], p3) - p3);
            double dist4 = bc->dist_bc(net->intpoint(flj[0], flj[1], p4) - p4);
            double mindist = std::min(std::min(dist1, dist2), std::min(dist3, dist4));

            if (mindist < cut) ref.push_back({fli, flj});
        }
    }

    std::vector<array<array<int, 2>, 2>> p;
    for (auto pair : pairs) {
        array<int, 2> fli = pair[0];
        array<int, 2> flj = pair[1];

        // minimum distance between springs
        vec_type p1 = net->get_start(flj[0], flj[1]);
        vec_type p2 = net->get_end(flj[0], flj[1]);
        vec_type p3 = net->get_start(fli[0], fli[1]);
        vec_type p4 = net->get_end(fli[0], fli[1]);
        double dist1 = bc->dist_bc(net->intpoint(fli[0], fli[1], p1) - p1);
        double dist2 = bc->dist_bc(net->intpoint(fli[0], fli[1], p2) - p2);
        double dist3 = bc->dist_bc(net->intpoint(flj[0], flj[1], p3) - p3);
        double dist4 = bc->dist_bc(net->intpoint(flj[0], flj[1], p4) - p4);
        double mindist = std::min(std::min(dist1, dist2), std::min(dist3, dist4));

        if (mindist < cut) p.push_back({fli, flj});
    }

    // check that spring-spring pairs in both lists are the same
    // we can do this since the spring addition order
    // is preserved in each quadrant
    assert(ref.size() == p.size());
    for (size_t i = 0; i < ref.size(); i++) {
        assert(ref[i] == p[i]);
    }
}

// end [checks]

// begin [add spring bc]
// implementations of adding a spring to quadrants for different boundary conditions

void quadrants::add_spring_nonperiodic(
        double xlo, double xhi, double ylo, double yhi, array<int, 2> fl)
{
    array<double, 2> fov = bc->get_fov();

    // compute quadrants containing bounding box

    int xlower = floor(nq[0] * (xlo / fov[0] + 0.5));
    int xupper = ceil(nq[0] * (xhi / fov[0] + 0.5));

    int ylower = floor(nq[1] * (ylo / fov[1] + 0.5));
    int yupper = ceil(nq[1] * (yhi / fov[1] + 0.5));

    assert(xlower <= xupper);
    assert(ylower <= yupper);

    // add spring index to quadrants

    for (int ii = xlower; ii <= xupper; ii++) {
        int i = clip(ii, 0, nq[0]);
        for (int jj = ylower; jj <= yupper; jj++) {
            int j = clip(jj, 0, nq[1]);
            quads[i][j].push_back(fl);
        }
    }
}

void quadrants::add_spring_xperiodic(
        double xlo, double xhi, double ylo, double yhi, array<int, 2> fl)
{
    array<double, 2> fov = bc->get_fov();

    // compute quadrants containing bounding box

    int xlower = floor(nq[0] * (xlo / fov[0] + 0.5));
    int xupper = ceil(nq[0] * (xhi / fov[0] + 0.5));

    int ylower = floor(nq[1] * (ylo / fov[1] + 0.5));
    int yupper = ceil(nq[1] * (yhi / fov[1] + 0.5));

    assert(xlower <= xupper);
    assert(ylower <= yupper);

    // add spring index to quadrants

    for (int ii = xlower; ii <= xupper; ii++) {
        int i = imod(ii, nq[0]);
        for (int jj = ylower; jj <= yupper; jj++) {
            int j = clip(jj, 0, nq[1]);
            quads[i][j].push_back(fl);
        }
    }
}

void quadrants::add_spring_periodic(
        double xlo, double xhi, double ylo, double yhi, array<int, 2> fl)
{
    array<double, 2> fov = bc->get_fov();

    // compute quadrants containing bounding box

    int xlower = floor(nq[0] * (xlo / fov[0] + 0.5));
    int xupper = ceil(nq[0] * (xhi / fov[0] + 0.5));

    int ylower = floor(nq[1] * (ylo / fov[1] + 0.5));
    int yupper = ceil(nq[1] * (yhi / fov[1] + 0.5));

    assert(xlower <= xupper);
    assert(ylower <= yupper);

    // add spring index to quadrants

    for (int ii = xlower; ii <= xupper; ii++) {
        int i = imod(ii, nq[0]);
        for (int jj = ylower; jj <= yupper; jj++) {
            int j = imod(jj, nq[1]);
            quads[i][j].push_back(fl);
        }
    }
}

void quadrants::add_spring_lees_edwards(
        double xlo, double xhi, double ylo, double yhi, array<int, 2> fl)
{
    array<double, 2> fov = bc->get_fov();
    double delrx = bc->get_delrx();

    // compute y quadrants containing bounding box
    // note that the x quadrants depend on the y quadrant value

    int ylower = floor(nq[1] * (ylo / fov[1] + 0.5));
    int yupper =  ceil(nq[1] * (yhi / fov[1] + 0.5));

    assert(ylower <= yupper);

    for (int jj = ylower; jj <= yupper; jj++) {
        int j = jj;

        // compute x quadrants containing bounding box
        // x quadrants in upper and lower images are shifted by delrx

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
        assert(0 <= j && j < nq[1]);

        int xlower = floor(nq[0] * (xlo_new / fov[0] + 0.5));
        int xupper =  ceil(nq[0] * (xhi_new / fov[0] + 0.5));

        assert(xlower <= xupper);

        for (int ii = xlower; ii <= xupper; ii++) {
            int i = imod(ii, nq[0]);


            quads[i][j].push_back(fl);
        }
    }
}

// end [add spring bc]
