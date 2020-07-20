#ifndef AFINES_VEC_H
#define AFINES_VEC_H

#include "globals.h"

using vec_type = array<double, 2>;

inline vec_type operator-(vec_type a)
{
    return {-a[0], -a[1]};
}

inline virial_type outer(vec_type a, vec_type b)
{
    return {{
        {{a[0] * b[0], a[0] * b[1]}},
        {{a[1] * b[0], a[1] * b[1]}}
    }};
}

inline vec_type operator-(vec_type a, vec_type b)
{
    return {a[0] - b[0], a[1] - b[1]};
}

inline vec_type operator+(vec_type a, vec_type b)
{
    return {a[0] + b[0], a[1] + b[1]};
}

inline vec_type operator*(double a, vec_type b)
{
    return {a * b[0], a * b[1]};
}

inline vec_type operator*(vec_type a, double b)
{
    return {a[0] * b,  a[1] * b};
}

inline vec_type operator/(vec_type a, double b)
{
    return {a[0] / b, a[1] / b};
}

#endif
