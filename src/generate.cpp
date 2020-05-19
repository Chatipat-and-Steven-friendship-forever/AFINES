#include "filament_ensemble.h"
#include "motor_ensemble.h"

vector<vector<double>> generate_filament_ensemble(
        box *bc, int npolymer, int nbeads_min, int nbeads_extra, int nbeads_extra_prob,
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
            x = rng(-0.5 * view[0] * fov[0], 0.5 * view[0] * fov[0]);
            y = rng(-0.5 * view[1] * fov[1], 0.5 * view[1] * fov[1]);
            phi = rng(0.0, 2.0 * pi);
        }
        beads.push_back({x, y, radius, i});
        int nbeads = nbeads_min + distribution(generator);
        for (int j = 1; j < nbeads; j++) {
            array<double, 2> pos = bc->pos_bc(
                    {x + l0 * cos(phi), y + l0 * sin(phi)},
                    {l0 * cos(phi) / dt, l0 * sin(phi) / dt}, dt);
            x = pos[0];
            y = pos[1];
            beads.push_back({x, y, radius, i});
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
            x = rng(-0.5 * view[0] * fov[0], 0.5 * view[0] * fov[0]);
            y = rng(-0.5 * view[1] * fov[1], 0.5 * view[1] * fov[1]);
            phi = rng(0.0, 2.0 * pi);
            // phi = atan2(1+x-y*y, -1-x*x+y); first example in Mathematica's streamplot documentation
        }
        beads.push_back({x, y, radius, i});
        for (int j = 1; j < nbeads; j++) {
            array<double, 2> pos = bc->pos_bc(
                    {x + l0 * cos(phi), y + l0 * sin(phi)},
                    {l0 * cos(phi) / dt, l0 * sin(phi) / dt}, dt);
            x = pos[0];
            y = pos[1];
            beads.push_back({x, y, radius, i});
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
            x = rng(-0.5 * (fov[0] * alpha - l0), 0.5 * (fov[0] * alpha - l0));
            y = rng(-0.5 * (fov[1] * alpha - l0), 0.5 * (fov[1] * alpha - l0));
            a = rng(0.0, 2.0 * pi);
        }
        double dx = l0 * cos(a);
        double dy = l0 * sin(a);
        motors.push_back(vector<double>{x - 0.5 * dx, y - 0.5 * dy, dx, dy, -1, -1, -1, -1});
    }
    return motors;
}
