/*
 *  globals.h
 *  
 *  Created by Simon Freedman on 9/12/2014
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef AFINES_GLOBALS_H
#define AFINES_GLOBALS_H

//=====================================
// forward declared dependencies

//=====================================
//included dependences
#include <boost/algorithm/string.hpp>
#include <boost/functional/hash.hpp>
#include <boost/filesystem.hpp>

#include <algorithm> //std::for_each
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream> //std::cout
#include <limits>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include "vec.h"
#include "pcg_random.hpp"

using namespace std;
namespace fs = boost::filesystem;

/* distances in microns, time in seconds, forces in pN * 
 * --> Temp in pN-um                                   */
//const double pi, eps, dt, temperature;
const double pi = 3.14159265358979323;
const double maxSmallAngle = 0.001; //pi/12.0; //Small angles DEFINED as such that sin(t) = t to 2 SigFigs
const double eps = 1e-4;
const double infty = 1e10;
const double actin_mass_density = 2.6e-14; //miligram / micron

// [random]

double rng_u();  // uniform distribution [0, 1)
double rng_n();  // normal distribution mean=0 std=1

pcg64 *get_rng();

// these methods handle rng state only
string gen_seed();  // generate rng state from random device
string get_seed();  // serialize rng state
void set_seed(string);  // seed from serialized rng state

// note that distributions also have state
// these methods handle both rng and distribution states
void rng_save(ostream &); // save entire state
void rng_load(istream &);  // load entire state

inline vec_type vec_randn()
{
    // separate statements to keep call order correct
    double x = rng_n();
    double y = rng_n();
    return {x, y};
}

// end [random]

/*generic functions to be used below*/

int pr(int num);

vector<int> range_bc(string bc, double delrx, int botq, int topq, int low, int high);
vector<int> range_bc(string bc, double delrx, int botq, int topq, int low, int high, int di);
vector<int> int_range(int lo, int hi);
vector<int> int_range(int lo, int hi, int di);

double mean_periodic(const vector<double>& nums, double bnd);
double mean(const vector<double>& nums);
array<double, 2> cm_bc(string bc, const vector<double>& xi, const vector<double>& yi, double xbox, double ybox, double shear_dist);

array<int, 2> coord2quad(const array<double, 2>& fov, const array<int, 2>& nq, const array<double, 2>& pos);
int coord2quad_floor(double fov, int nq, double pos);
int coord2quad_ceil(double fov, int nq, double pos);
int coord2quad(double fov, int nq, double pos);

double angBC(double ang);
double angBC(double ang, double max);

double var(const vector<double>& vals);
double mode_var(const vector<double>& vals, double m);
bool close(double e, double a, double r);
bool are_same(double a, double b);
vector<double> sum_vecs(const vector<double>& v1, const vector<double>& v2);
vector<double *> vec2ptrvec(const vector<double>&, int dim);
vector<double *> str2ptrvec(string, string, string);
vector<array<double,3> > str2arrvec(string, string, string);
vector<vector<double> > file2vecvec(string path, string delim);
vector<vector<double> > traj2vecvec(string path, string delim, double tf);
double get_restart_strain(string path, int tf);
double last_full_timestep(string dirpath);
void write_first_nlines(string path, int n);
void write_first_ntsteps(string path, int n);
void write_first_tsteps(string path, double tstop);

template <typename T> int sgn(T val);
int mysgn(double);

static inline int imod(int a, int b)
{
    while (a < 0) a += b;
    while (a >= b) a -= b;
    assert(0 <= a && a < b);
    return a;
}

static inline int clip(int x, int lo, int hi)
{
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    assert(lo <= x && x <= hi);
    return x;
}

pair<double, array<int, 2> > flip_pair(const pair<array<int, 2>, double> &p);
multimap<double, array<int, 2> > flip_map(const std::unordered_map<array<int, 2>, double, boost::hash<array<int,2>>> &p);
string print_pair(string name, const array<double, 2>& p);

//template <typename A, typename B> pair<B,A> flip_pair(const pair<A,B> &p);
//template <typename A, typename B> multimap<B,A> flip_map(const map<A,B> &src);

map<array<int, 2>, double> transpose(map<array<int, 2>, double> mat);
map<array<int, 2>, double> invert_block_diagonal(map<array<int, 2>, double> mat);
void intarray_printer(array<int,2> a);

std::optional<vec_type> seg_seg_intersection(vec_type, vec_type, vec_type, vec_type);
std::string quads_error_message(std::string, vector<array<int, 2> >, vector<array<int, 2> > );

struct mc_prob
{
    mc_prob()
    {
        prob = rng_u();
        used = 0.0;
    }

    double prob;
    double used;

    std::optional<double> operator()(double needprob)
    {
        std::optional<double> result;
        if (used <= prob && prob < used + needprob) {
            result = {prob - used};
        }
        used += needprob;
        if (used > 1) {
            throw runtime_error("Need too much probability.");
        }
        return result;
    }

};

struct fp_index_type
{
    int f_index;
    int p_index;
};

#endif
