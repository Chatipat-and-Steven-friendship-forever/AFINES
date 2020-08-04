#ifndef AFINES_QUADRANTS_H
#define AFINES_QUADRANTS_H

#include "globals.h"

class quadrants {
    public:
        quadrants(class box *bc, array<int, 2> nq);
        ~quadrants();
        void use_quad(bool flag);
        void use_all(bool flag);

        void add_spring(vec_type h0, vec_type disp, array<int, 2> fl);
        vector<array<int, 2>> *get_attach_list(vec_type pos);
        void build_pairs();
        vector<array<array<int, 2>, 2>> *get_pairs();
        void clear();

        array<int, 2> get_nq() { return nq; }
        vector<array<int, 2>> *get_quad(array<int, 2> q) { return &quads[q[0]][q[1]]; }

        double get_cut();
        double get_paircut();

        // check
        void check_duplicates();
        void check_quad(class filament_ensemble *, vector<array<int, 2>> *, vec_type);
        void check_pairs(class filament_ensemble *);

    protected:
        void add_spring_nonperiodic(double xlo, double xhi, double ylo, double yhi, array<int, 2>);
        void add_spring_xperiodic(double xlo, double xhi, double ylo, double yhi, array<int, 2>);
        void add_spring_periodic(double xlo, double xhi, double ylo, double yhi, array<int, 2>);
        void add_spring_lees_edwards(double xlo, double xhi, double ylo, double yhi, array<int, 2>);

        class box *bc;
        array<int, 2> nq;
        bool all_flag, quad_flag;
        vector<array<int, 2>> all_springs;
        vector<array<int, 2>> **quads;
        vector<array<array<int, 2>, 2>> pairs;
        unordered_set<array<array<int, 2>, 2>, boost::hash<array<array<int, 2>, 2>>> pairset;
};

#endif
