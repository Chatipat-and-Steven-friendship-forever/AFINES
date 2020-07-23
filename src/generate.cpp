#include "filament_ensemble.h"
#include "motor_ensemble.h"

vector<vector<double>> generate_filament_ensemble(
        box *bc, int npolymer, int nbeads_min, int nbeads_extra, double nbeads_extra_prob,
        double dt, double temp, double radius, double l0,
        vector<array<double, 3>> pos_sets, double kb, double seed)
{
    bool is_straight = false;
    if (seed == -1) {
        is_straight = true;
    } else {
        srand(seed);
    }
    binomial_distribution<int> distribution(nbeads_extra, nbeads_extra_prob);
    default_random_engine generator(seed + 2);
    array<double, 2> fov = bc->get_fov();
    double var = (temp != 0.0) ? temp / kb : 0.0;
    array<double, 2> view = {1.0, 1.0};

    vector<vector<double>> beads;
    for (int i = 0; i < npolymer; i++) {
        double x, y, phi;
        if (i < int(pos_sets.size())) {
            x = pos_sets[i][0];
            y = pos_sets[i][1];
            phi = pos_sets[i][2];
        } else {
            x = view[0] * fov[0] * (rng_u() - 0.5);
            y = view[1] * fov[1] * (rng_u() - 0.5);
            phi = 2.0 * pi * rng_u();
        }
        beads.push_back({x, y, radius, double(i)});
        int nbeads = nbeads_min + distribution(generator);
        for (int j = 1; j < nbeads; j++) {
            vec_type pos = bc->pos_bc({x + l0 * cos(phi), y + l0 * sin(phi)});
            x = pos.x;
            y = pos.y;
            beads.push_back({x, y, radius, double(i)});
            if (!is_straight) phi += sqrt(var) * rng_n();
        }
    }
    return beads;
}

vector<vector<double>> generate_filament_ensemble(
        box *bc, double density, double dt, double temp, double radius, int nbeads, double l0,
        vector<array<double, 3>> pos_sets, double kb, double seed)
{
    bool is_straight = false;
    if (seed == -1) {
        is_straight = true;
    } else {
        srand(seed);
    }
    array<double, 2> fov = bc->get_fov();
    int npolymer = ceil(density * fov[0] * fov[1]) / nbeads;
    double var = (temp != 0.0) ? temp / kb : 0.0;
    array<double, 2> view = {1.0, 1.0};

    vector<vector<double>> beads;
    for (int i = 0; i < npolymer; i++) {
        double x, y, phi;
        if (i < int(pos_sets.size())) {
            x = pos_sets[i][0];
            y = pos_sets[i][1];
            phi = pos_sets[i][2];
        } else {
            x = view[0] * fov[0] * (rng_u() - 0.5);
            y = view[1] * fov[1] * (rng_u() - 0.5);
            phi = 2.0 * pi * rng_u();
            // phi = atan2(1+x-y*y, -1-x*x+y); first example in Mathematica's streamplot documentation
        }
        beads.push_back({x, y, radius, double(i)});
        for (int j = 1; j < nbeads; j++) {
            vec_type pos = bc->pos_bc({x + l0 * cos(phi), y + l0 * sin(phi)});
            x = pos.x;
            y = pos.y;
            beads.push_back({x, y, radius, double(i)});
            if (!is_straight) phi += sqrt(var) * rng_n();
        }
    }
    return beads;
}

vector<vector<double>> generate_motor_ensemble(box *bc, double density, double l0, vector<array<double, 3>> positions)
{
    array<double, 2> fov = bc->get_fov();
    int nm = int(ceil(density * fov[0] * fov[1]));
    double alpha = 1;

    vector<vector<double>> motors;
    for (int i = 0; i < nm; i++) {
        double x, y, a;
        if (i < int(positions.size())) {
            x = positions[i][0];
            y = positions[i][1];
            a = positions[i][2];
        } else {
            x = (fov[0] * alpha - l0) * (rng_u() - 0.5);
            y = (fov[1] * alpha - l0) * (rng_u() - 0.5);
            a = 2.0 * pi * rng_u();
        }
        double dx = l0 * cos(a);
        double dy = l0 * sin(a);
        motors.push_back({{x - 0.5 * dx, y - 0.5 * dy, dx, dy, -1, -1, -1, -1}});
    }
    return motors;
}

vector<vector<double>> spring_spring_intersections(box *bc, vector<vector<double>> beads, double len, double prob)
{
    vector<vector<vec_type>> filaments;
    vector<vec_type> current_filament;
    int last_f_id = 0;
    for (size_t i = 0; i < beads.size(); i++) {
        if (beads[i][3] != last_f_id && current_filament.size() > 0) {
            filaments.push_back(current_filament);
            current_filament.clear();
        }
        current_filament.push_back({beads[i][0], beads[i][1]});
        last_f_id = beads[i][3];
    }
    if (current_filament.size() > 0) {
        filaments.push_back(current_filament);
        current_filament.clear();
    }

    vector<vector<double>> motors;
    for (size_t f1 = 0; f1 < filaments.size(); f1++) {
        for (size_t l1 = 0; l1 < filaments[f1].size() - 1; l1++) {
            vec_type r1 = filaments[f1][l1];
            vec_type r2 = filaments[f1][l1 + 1];
            for (size_t f2 = f1 + 1; f2 < filaments.size(); f2++) {
                for (size_t l2 = 0; l2 < filaments[f2].size() - 1; l2++) {
                    if (f1 == f2 && abs(int(l1) - int(l2)) < 2) continue;
                    vec_type s1 = filaments[f2][l2];
                    vec_type s2 = filaments[f2][l2 + 1];
                    boost::optional<vec_type> inter = seg_seg_intersection_bc(bc, r1, r2, s1, s2);
                    if (inter && rng_u() <= prob) {
                        vec_type disp = bc->rij_bc(s2 - s1);
                        double llen = abs(disp);
                        vec_type direc = disp / llen;
                        motors.push_back({inter->x, inter->y, len * direc.x, len * direc.y,
                                double(f1), double(f2), double(l1), double(l2)});
                    }
                }
            }
        }
    }
    return motors;
}
