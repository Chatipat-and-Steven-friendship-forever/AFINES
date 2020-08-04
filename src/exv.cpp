#include "exv.h"

#include "box.h"
#include "quadrants.h"
#include "filament_ensemble.h"

void excluded_volume::update_spring_forces_from_quads(quadrants *quads, filament_ensemble *net)
{
    pe_exv = 0.0;
    vir_exv.zero();

    for (array<array<int, 2>, 2> pair : *quads->get_pairs()) {
        // adjacent springs would yield excluded volume interactions between the same bead
        if (pair[0][0] == pair[1][0] && abs(pair[0][1] - pair[1][1]) < 2) continue;

        update_force_between_filaments(net, pair[0], pair[1]);
    }
}

void excluded_volume::update_spring_forces(filament_ensemble *net)
{
    // loop over all unique pairs of springs
    // which are not on the same filament

    // for every spring in filament f
    for (int f = 0; f < net->get_nfilaments(); f++) {
        for (int i = 0; i < net->get_nsprings(f); i++) {

            // for every spring in filaments g > f
            for (int g = f + 1; g < net->get_nfilaments(); g++) {
                for (int j = 0; j < net->get_nsprings(g); j++) {

                    update_force_between_filaments(net, {f, i}, {g, j});
                }
            }
        }
    }
}

void excluded_volume::update_force_between_filaments(
        filament_ensemble *net, array<int, 2> fl1, array<int, 2> fl2)
{
    // This function calculates the forces applied to the actin beads of a pair of filaments under certain limits.
    // Here, we use distance of closest approach to describe the direction and magnitude of the forces.

    double b = 1/rmax;

    int n1 = fl1[0];
    int l1 = fl1[1];

    int n2 = fl2[0];
    int l2 = fl2[1];

    vec_type h0_1 = net->get_start(n1, l1);
    vec_type h1_1 = net->get_end(n1, l1);

    vec_type h0_2 = net->get_start(n2, l2);
    vec_type h1_2 = net->get_end(n2, l2);

    array<double, 2> len;
    len[0] = net->get_llength(n1, l1);
    len[1] = net->get_llength(n2, l2);

    // compute the nearest point on the other spring
    vec_type p1 = net->intpoint(n1, l1, h0_2);
    vec_type p2 = net->intpoint(n1, l1, h1_2);
    vec_type p3 = net->intpoint(n2, l2, h0_1);
    vec_type p4 = net->intpoint(n2, l2, h1_1);

    // compute the distance to the nearest point
    array<double, 4> r_c;
    r_c[0] = bc->dist_bc(h0_2 - p1);
    r_c[1] = bc->dist_bc(h1_2 - p2);
    r_c[2] = bc->dist_bc(h0_1 - p3);
    r_c[3] = bc->dist_bc(h1_1 - p4);

    // find the minimum distance between the filaments
    // assuming that they don't intersect
    double r = r_c[0];
    int index = 0;
    for(int k = 1; k < 4; k++){
        if (r_c[k] < r) {
            r = r_c[k];
            index = k;
        }
    }

    bool intersect = net->intersect({n1, l1}, {n2, l2});
    if (intersect) throw runtime_error("Intersecting filaments with excluded volume!");

    if (r < rmax) {

        double length = 0.0;
        double len1 = 0.0;
        vec_type dist;
        if (index == 0) {
            r = r_c[0];
            len1 = bc->dist_bc(h0_1 - p1);
            length = len[0];
            dist = bc->rij_bc(p1 - h0_2);
        } else if (index == 1) {
            r = r_c[1];
            len1 = bc->dist_bc(h0_1 - p2);
            length = len[0];
            dist = bc->rij_bc(p2 - h1_2);
        } else if (index == 2) {
            r = r_c[2];
            len1 = bc->dist_bc(h0_2 - p3);
            length = len[1];
            dist = bc->rij_bc(p3 - h0_1);
        } else if (index == 3) {
            r = r_c[3];
            len1 = bc->dist_bc(h0_2 - p4);
            length = len[1];
            dist = bc->rij_bc(p4 - h1_1);
        }

        double len2 = length - len1;
        double r_1 = len2/length;
        double r_2 = len1/length;

        vec_type F = 2*kexv*dist*b*((1/r) - b);

        pe_exv += kexv*pow((1-r*b),2);
        vir_exv += -0.5 * outer(dist, F);

        if (index == 0) {
            net->update_forces(n1, l1, F*r_1);
            net->update_forces(n1, l1+1, F*r_2);
            net->update_forces(n2, l2, -F);
        } else if (index == 1) {
            net->update_forces(n1, l1, F*r_1);
            net->update_forces(n1, l1+1, F*r_2);
            net->update_forces(n2, l2+1, -F);
        } else if (index == 2) {
            net->update_forces(n2, l2, F*r_1);
            net->update_forces(n2, l2+1, F*r_2);
            net->update_forces(n1, l1, -F);
        } else if (index == 3) {
            net->update_forces(n2, l2, F*r_1);
            net->update_forces(n2, l2+1, F*r_2);
            net->update_forces(n1, l1+1, -F);
        }

    }
}

void excluded_volume::update_excluded_volume(filament_ensemble *net)
{
    // 10^6 included to account for m to micron conversion
    double a = 0.004;
    double b = 1/rmax;

    // loop over unique pairs of beads not on the same filament

    // for each bead i on filament f
    for (int f = 0; f < net->get_nfilaments(); f++) {
        for (int i = 0; i < net->get_nbeads(f); i++) {
            vec_type h1 = net->get_pos(f, i);

            // for each bead j on filament g > f
            for (int g = f + 1; g < net->get_nfilaments(); g++) {
                for (int j = 0; j < net->get_nbeads(g); j++) {
                    vec_type h2 = net->get_pos(g, j);

                    // compute forces for the potential
                    // U(r) = a (b r - 1)^2 if r < 1 / b
                    vec_type del = h1 - h2;
                    double r = bc->dist_bc(del);
                    if (r > 0.0 && r <= rmax) {
                        vec_type F = 2*del*a*b*((1/r)-b);
                        net->update_forces(f, i, F);
                        net->update_forces(g, j, -F);
                    }
                }
            }
        }
    }
}
