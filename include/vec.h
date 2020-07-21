#ifndef AFINES_VEC_H
#define AFINES_VEC_H

struct vec_type {
    double x, y;

    vec_type()
    {
        x = y = 0.0;
    }

    vec_type(double x_, double y_)
    {
        x = x_; y = y_;
    }

    void zero()
    {
        x = y = 0.0;
    }

    vec_type operator-()
    {
        return {-x, -y};
    }

    vec_type &operator+=(vec_type other)
    {
        x += other.x;
        y += other.y;
        return *this;
    }

    vec_type &operator-=(vec_type other)
    {
        x -= other.x;
        y -= other.y;
        return *this;
    }

};

struct virial_type {
    double xx, xy;
    double yx, yy;

    virial_type()
    {
        xx = xy = yx = yy = 0.0;
    }

    virial_type(double xx_, double xy_, double yx_, double yy_)
    {
        xx = xx_; xy = xy_;
        yx = yx_; yy = yy_;
    }

    void zero()
    {
        xx = xy = yx = yy = 0.0;
    }

    virial_type &operator+=(virial_type other)
    {
        xx += other.xx; xy += other.xy;
        yx += other.yx; yy += other.yy;
        return *this;
    }

};

inline virial_type outer(vec_type a, vec_type b)
{
    return {a.x * b.x, a.x * b.y,
            a.y * b.x, a.y * b.y};
}

inline vec_type operator-(vec_type a, vec_type b)
{
    return {a.x - b.x, a.y - b.y};
}

inline vec_type operator+(vec_type a, vec_type b)
{
    return {a.x + b.x, a.y + b.y};
}

inline vec_type operator*(double a, vec_type b)
{
    return {a * b.x, a * b.y};
}

inline vec_type operator*(vec_type a, double b)
{
    return {a.x * b,  a.y * b};
}

inline vec_type operator/(vec_type a, double b)
{
    return {a.x / b, a.y / b};
}

inline double dot(vec_type a, vec_type b)
{
    return a.x * b.x + a.y * b.y;
}

inline double cross(vec_type a, vec_type b)
{
    return a.x * b.y - a.y * b.x;
}

inline double abs2(vec_type a)
{
    return a.x * a.x + a.y * a.y;
}

inline double abs(vec_type a)
{
    return sqrt(abs2(a));
}

inline bool operator<(vec_type a, vec_type b)
{
    return a.x < b.x || (a.x == b.x && a.y < b.y);
}

#endif
