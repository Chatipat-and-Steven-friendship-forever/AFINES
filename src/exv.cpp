#include "exv.h"

void excluded_volume::update_spring_forces_from_quads(quadrants *quads, vector<filament *> &network)
{
    for (array<array<int, 2>, 2> pair : *quads->get_pairs()) {
        int f1 = pair[0][0];
        int l1 = pair[0][1];
        int f2 = pair[1][0];
        int l2 = pair[1][1];

        // adjacent springs would yield excluded volume interactions between the same bead
        if (f1 == f2 && abs(l1 - l2) < 2) continue;

        update_force_between_filaments(network, f1, l2, f2, l2);
    }
}

void excluded_volume::update_spring_forces(vector<filament *> &network, int f)
{
    //This function loops through every filament and spring in the network and applies the force calulation under certain limits
    int net_sz = network.size();

    // for every spring in filament f
    int lks_sz = network[f]->get_nsprings();
    for (int i = 0; i < lks_sz; i++) {

        // for every spring in filaments g > f
        for (int g = f+1; g < net_sz; g++) {
            int oth_lks_sz = network[g]->get_nsprings();
            for (int j = 0; j < oth_lks_sz; j++) {

                update_force_between_filaments(network, f, i, g, j);
            }
        }
    }
}

void excluded_volume::update_force_between_filaments(
        vector<filament *> &network, int n1, int l1, int n2, int l2)
{
    //This function calculates the forces applied to the actin beads of a pair of filaments under certain limits.
    //Here, we use distance of closest approach to describe the direction and magnitude of the forces.

    double b = 1/rmax;

    spring *s1 = network[n1]->get_spring(l1);
    spring *s2 = network[n2]->get_spring(l2);

    vec_type h0_1 = s1->get_h0();
    vec_type h1_1 = s1->get_h1();

    vec_type h0_2 = s2->get_h0();
    vec_type h1_2 = s2->get_h1();

    array<double, 2> len;
    len[0] = s1->get_length();
    len[1] = s2->get_length();

    // compute the nearest point on the other spring
    vec_type p1 = s1->intpoint(h0_2);
    vec_type p2 = s1->intpoint(h1_2);
    vec_type p3 = s2->intpoint(h0_1);
    vec_type p4 = s2->intpoint(h1_1);

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
        if(r_c[k] < r){
            r = r_c[k];
            index = k;
        }
    }

    bool intersect = s1->get_line_intersect(s2);

    // assume spring length l0 < 2 rmax ?
    if (r < rmax) {

        if (!intersect) {
            // doesn't intersect

            double length=0, len1=0;
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

            if (index == 0) {
                network[n1]->update_forces(l1, F*r_1);
                network[n1]->update_forces(l1+1, F*r_2);
                network[n2]->update_forces(l2, -F);
            } else if (index == 1) {
                network[n1]->update_forces(l1, F*r_1);
                network[n1]->update_forces(l1+1, F*r_2);
                network[n2]->update_forces(l2+1, -F);
            } else if (index == 2) {
                network[n2]->update_forces(l2, F*r_1);
                network[n2]->update_forces(l2+1, F*r_2);
                network[n1]->update_forces(l1, -F);
            } else if (index == 3) {
                network[n2]->update_forces(l2, F*r_1);
                network[n2]->update_forces(l2+1, F*r_2);
                network[n1]->update_forces(l1+1, -F);
            }

        } else {
            // intersects

            // apply constant force?
            // this doesn't look right
            vec_type F = {2*kexv/(rmax*sqrt(2)), 2*kexv/(rmax*sqrt(2))};

            pe_exv += kexv*pow((1-r*b),2);

            network[n1]->update_forces(l1, F);
            network[n1]->update_forces(l1+1, F);
            network[n2]->update_forces(l2, -F);
            network[n2]->update_forces(l2+1, -F);

        }
    }
}

void excluded_volume::update_excluded_volume(vector<filament *> &network, int f)
{
    //For every filament bead on f, for every bead not on f, calculate the force between the two bead using the Jones potential, and update them ( maybe divide by half due to overcaluclations).

    int net_sz = network.size();
    //10^6 included to account for m to microm conversion
    double a = 0.004;
    double b = 1/rmax;

    // for every bead in filament f
    int act_sz = network[f]->get_nbeads();
    for (int i = 0; i < act_sz; i++) {
        vec_type h1 = network[f]->get_bead_position(i);

        // for every bead in filament g > f
        for (int g = f+1; g < net_sz; g++) {
            int act_sz_other = network[g]->get_nbeads();
            for (int j = 0; j < act_sz_other; j++) {
                vec_type h2 = network[g]->get_bead_position(j);

                // compute forces for the potential
                // U(r) = a (b r - 1)^2 if r < 1 / b
                vec_type del = h1 - h2;
                double r = bc->dist_bc(del);
                if (r > 0.0 && r <= rmax) {
                    vec_type F = 2*del*a*b*((1/r)-b);
                    network[f]->update_forces(i,F);
                    network[g]->update_forces(j,-F);
                }
            }
        }
    }
}
