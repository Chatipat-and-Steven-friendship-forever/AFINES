#ifndef AFINES_QUADRANTS_H
#define AFINES_QUADRANTS_H

#include "globals.h"
#include "box.h"
#include "spring.h"

class quadrants {
    public:
        quadrants(box *bc, array<int, 2> nq);
        ~quadrants();
        void use_quad(bool flag);

        void add_spring(spring *s, array<int, 2> fl);
        vector<array<int, 2>> *get_attach_list(vec_type pos);
        void clear();

        array<int, 2> get_nq() { return nq; }
        vector<array<int, 2>> *get_quad(array<int, 2> q) { return &quads[q[0]][q[1]]; }

        void check_duplicates();

    protected:
        void add_spring_nonperiodic(spring *, array<int, 2>);
        void add_spring_periodic(spring *, array<int, 2>);

        box *bc;
        array<int, 2> nq;
        bool quad_flag;
        vector<array<int, 2>> all_springs;
        vector<array<int, 2>> **quads;
};

#endif
