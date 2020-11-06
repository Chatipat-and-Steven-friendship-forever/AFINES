/*------------------------------------------------------------------
 globals.cpp : functions and constants that are used throughout AFiNeS

 Copyright (C) 2016
 Created by: Simon Freedman, Shiladitya Banerjee, Glen Hocky, Aaron Dinner
 Contact: dinner@uchicago.edu

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version. See ../LICENSE for details.
-------------------------------------------------------------------*/

#include "globals.h"

// begin [random]

pcg64 rng;
std::uniform_real_distribution<double> uniform_dist;
std::normal_distribution<double> normal_dist;

double rng_u()
{
    return uniform_dist(rng);
}

double rng_n()
{
    return normal_dist(rng);
}

pcg64 *get_rng()
{
    return &rng;
}

string gen_seed()
{
    pcg_extras::seed_seq_from<std::random_device> seed_source;
    return fmt::format("{}", pcg64(seed_source));
}

void set_seed(string seed)
{
    std::istringstream(seed + '\n') >> rng;
}

string get_seed()
{
    return fmt::format("{}", rng);
}

void rng_load(istream &inp)
{
    inp >> rng >> uniform_dist >> normal_dist;
}

void rng_save(ostream &out)
{
    out << rng << '\n' << uniform_dist << '\n' << normal_dist << '\n';
}

// end [random]

int pr(int num)
{
    if (num==0) {
        return 1;
    }
    else {
        return 0;
    }
}

array<double, 2> cm_bc(string bc, const vector<double>& xi, const vector<double>& yi, double xbox, double ybox, double delrx)
{
    if (bc == "PERIODIC" || bc == "LEES-EDWARDS")
        return {{mean_periodic(xi, xbox) , mean_periodic(yi, ybox)}};
    else
        return {{mean(xi), mean(yi)}};
}

double mean(const vector<double>& nums)
{
    double tot = 0;
    for (double n : nums) tot += n;
    return tot/((double) nums.size());
}

// Source https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
double mean_periodic(const vector<double>& nums, double bnd)
{
    double theta, xitot=0, zetatot=0;
    for (double n : nums){
        theta = n*2*pi/bnd;
        xitot += cos(theta);
        zetatot += sin(theta);
    }
    double thetabar = atan2(zetatot/((double) nums.size()), xitot/((double) nums.size())) + pi;
    return bnd*thetabar/(2*pi);
}

double var(const vector<double>& vals)
{
    double m = mean(vals), sum = 0;
    for (unsigned int i = 0; i < vals.size(); i++){
        sum += (vals[i] - m)*(vals[i] - m);
    }
    return sum / vals.size();
}

double mode_var(const vector<double>& vals, double m)
{
    double sum = 0;
    for (unsigned int i = 0; i < vals.size(); i++){
        sum += (vals[i] - m)*(vals[i] - m);
    }
    return sum / vals.size();
}

vector<double> sum_vecs(const vector<double>& v1, const vector<double>& v2)
{
    vector<double> s;
    if (v1.empty())
        s = v2;
    else if( v2.empty())
        s = v1;
    else if (v1.size() != v2.size())
        return s;
    else{
        for (unsigned int i = 0; i < v1.size(); i++){
            s.push_back(v1[i] + v2[i]);
        }
    }
    return s;
}

//SO:17333 (more info there)
bool are_same(double a, double b)
{
    return fabs(a-b) < std::numeric_limits<double>::epsilon();
}

bool close(double actual, double expected, double err)
{
    if (expected == 0){
        return fabs(expected-actual) < err;
    }
    else{
        return fabs(expected-actual)/expected < err;
    }
}

/* Takes a vector formatted
 * [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
 * And converts it to a vector formatted
 * [{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}] if dim = 4
 * Or
 * [{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}] if dim = 3
 */

vector<double *> vec2ptrvec(const vector<double>& v, int dim)
{
    vector<double *> out;
    double * pos;
    for (unsigned int i = 0; i < v.size(); i+=dim)
    {
        pos = new double[dim];
        for (int j = 0; j < dim; j++)
        {
           pos[j] = v[i+j];
        }
        out.push_back(pos);
    }
    return out;
}

/* Takes a string formatted
 * 1,2,3;4,5,6;7,8,9;10,11,12
 * and converts it into a vector of pointers:
 * [{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}]
 */

vector<double *> str2ptrvec(string pos_str, string pos_dlm, string coord_dlm)
{
    vector<string> posn, posns, coords;
    vector<double *> out;
    double * pos;

    boost::split(posns, pos_str, boost::is_any_of(pos_dlm));

    for(unsigned int i=0; i < posns.size(); i++){

        boost::split(coords, posns[i], boost::is_any_of(coord_dlm));
        pos = new double[coords.size()];

        for(unsigned int j=0; j < coords.size(); j++){
            pos[j] = (double) atof(coords[j].data());
        }
        out.push_back(pos);
    }

    return out;
}

vector<array<double,3> > str2arrvec(string pos_str, string pos_dlm, string coord_dlm)
{
    vector<string> posn, posns, coords;
    vector<array<double,3> > out;
    array<double,3> pos;

    boost::split(posns, pos_str, boost::is_any_of(pos_dlm));

    for(unsigned int i=0; i < posns.size(); i++){

        boost::split(coords, posns[i], boost::is_any_of(coord_dlm));

        for(unsigned int j=0; j < 3; j++) pos[j] = (double) atof(coords[j].data());

        out.push_back(pos);
    }

    return out;
}

vector<vector<double> > file2vecvec(string path, string delim)
{
    vector<vector<double> > out;
    string pos_str = "";
    vector<string> coords;
    vector<double> pos;

    ifstream pos_file;
    pos_file.open(path);

    while(getline(pos_file, pos_str))
    {
        boost::trim_right(pos_str);
        boost::split(coords, pos_str, boost::is_any_of(delim));

        for(unsigned int j=0; j < coords.size(); j++)
            pos.push_back( (double) atof(coords[j].data()) );

        out.push_back(pos);

        pos.clear();
    }

    pos_file.close();

    return out;
}

void intarray_printer(array<int, 2> a)
{
    cout<<"\n{ " <<a[0]<<" , "<<a[1]<<" }";
}

template <typename T> int sgn(T val){
    return (T(0) < val) - (val < T(0));
}

int mysgn(double d){
    return d < 0 ? -1 : 1;
}

vector<int> int_range(int lo, int hi)
{
    vector<int> out;
    for (int i = lo; i<= hi; i++) out.push_back(i);
    return out;
}

vector<int> int_range(int lo, int hi, int di)
{
    vector<int> out;
    for (int i = lo; i != hi + di; i+=di) out.push_back(i);
    return out;
}

vector<int> range_bc(string bc, double delrx, int botq, int topq, int lo, int hi)
{
    vector<int> out;
    if (hi > topq) hi = hi - (topq-botq);
    if (lo < botq) lo = lo + (topq-botq);

    if (lo <= hi)
        out = int_range(lo, hi);
    else if (bc == "PERIODIC" || bc == "LEES-EDWARDS"){
        vector<int> A = int_range(lo, topq-1), B = int_range(botq, hi);
        out.reserve(A.size() + B.size());
        out.insert(out.end(), A.begin(), A.end());
        out.insert(out.end(), B.begin(), B.end());
    }
    else
        out = vector<int>();

    return out;
}

vector<int> range_bc(string bc, double delrx, int botq, int topq, int lo, int hi, int di)
{
    vector<int> out;
    if (hi > topq) hi = hi - (topq-botq);
    if (lo < botq) lo = lo + (topq-botq);

    if ((lo <= hi && di > 0) || (lo > hi && di < 0))
        out = int_range(lo, hi, di);
    else if (bc == "PERIODIC" || bc == "LEES-EDWARDS"){
        vector<int> A, B;
        if ( di > 0 ){
            A = int_range(lo, topq, di);
            B = int_range(botq, hi, di);
        }else{
            A = int_range(lo, botq, di);
            B = int_range(topq - 1, hi, di);
        }

        out.reserve(A.size() + B.size());
        out.insert(out.end(), A.begin(), A.end());
        out.insert(out.end(), B.begin(), B.end());
    }
    else
        out = vector<int>();

    return out;
}

// Method to sort a map by value; source, for more general formulation:
// http://stackoverflow.com/questions/5056645/sorting-stdmap-using-value/5056797#5056797
pair<double, array<int, 2> > flip_pair(const pair<array<int, 2>, double> &p)
{
        return std::pair<double,array<int,2> >(p.second, p.first);
}

multimap<double, array<int, 2> > flip_map(const unordered_map<array<int, 2>, double, boost::hash<array<int,2>>> &src)
{
    multimap<double,array<int,2> > dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()),
            flip_pair);
    return dst;
}

std::optional<vec_type> seg_seg_intersection(
        vec_type p1, vec_type q1,
        vec_type p2, vec_type q2)
{
    // adapted from Graphics Gems, page 304

    // intersection is at
    // p1 + t v1 = p2 + s v2

    vec_type v1 = q1 - p1;
    vec_type v2 = q2 - p2;

    double det = cross(v1, v2);

    if (det == 0.0) return {};  // springs are parallel
    // ignore the collinear case, which should be very rare

    vec_type p12 = p2 - p1;
    double t = cross(p12, v2) / det;
    double s = cross(p12, v1) / det;

    // check if intersection point is on both segments
    if (0.0 <= t && t <= 1.0 && 0.0 <= s && s <= 1.0) {
        return {p1 + t * v1};
    } else {
        return {};
    }
}

string print_pair(string name, const array<double, 2>& p)
{
    return name + ": ("+std::to_string(p[0])+","+std::to_string(p[1])+")";
}

int coord2quad_floor(double fov, int nq, double coord)
{
    int q = floor(nq*(coord/fov+0.5));
    if (q == nq)
        return 0;
    else
        return q;
}

int coord2quad_ceil(double fov, int nq, double coord)
{
    int q = ceil(nq*(coord/fov + 0.5));
    if (q == nq)
        return 0;
    else
        return q;
}

int coord2quad(double fov, int nq, double coord)
{
    int q = round(nq*(coord/fov + 0.5));
    if (q >= nq)
        return 0;
    else
        return q;
}

double angBC(double ang)
{
    return ang - 2*pi*floor(ang / (2*pi) + 0.5);
}

double angBC(double ang, double max)
{
    return ang - max*floor(ang / max + 0.5);
}

std::string quads_error_message(std::string title, vector<array<int, 2> > equads, vector< array<int, 2> > aquads)
{

    cout<<"\nTEST "<< title<< ": Expected Quadrants : don't equal spring Quadrants : \n";
    cout<<"\nActual Quadrants:";
    for_each(aquads.begin(), aquads.end(), intarray_printer);
    cout<<"\nExpected Quadrants:";
    for_each(equads.begin(), equads.end(), intarray_printer);
    return "";
}

vector<vector<double> > traj2vecvec(string path, string delim, double tf)
{
    vector<vector<double> > out;
    string pos_str = "";
    vector<string> coords;
    vector<double> pos;

    ifstream pos_file;
    pos_file.open(path);

    double t = 0;

    while(getline(pos_file, pos_str))
    {
        if (pos_str[0]=='t'){
            boost::split(coords, pos_str, boost::is_any_of(delim));
            t = (double) atof(coords[2].data());
            continue;
        }
        if (t < tf) continue;
        else if (t > tf) break;
        else{
            boost::trim_right(pos_str);
            boost::split(coords, pos_str, boost::is_any_of(delim));

            for(unsigned int j=0; j < coords.size(); j++)
                pos.push_back( (double) atof(coords[j].data()) );

            out.push_back(pos);

            pos.clear();
        }
    }

    pos_file.close();

    return out;
}

double get_restart_strain(string path, int tf)
{

  double restart_strain = 0;
  string str = "";
  vector<string> coords;

  ifstream pe_file;
  pe_file.open(path);

  for (int i = 1; i <= tf; i++)
    getline(pe_file, str);

  boost::trim_right(str);
  boost::split(coords, str, boost::is_any_of("\t "));

  restart_strain = (double) atof(coords[3].data());

  pe_file.close();

  return restart_strain;
}

double last_full_timestep(string path)
{
    string pos_str = "";
    vector<string> coords;

    ifstream pos_file;
    pos_file.open(path);

    int n = 0, pcount = 0;
    double t = 0, tprev = 0;

    while(getline(pos_file, pos_str))
    {
        if (pos_str[0]=='t'){
            boost::trim_right(pos_str);
            boost::split(coords, pos_str, boost::is_any_of("\t "));
            tprev = t;

            t = (double) atof(coords[2].data());
            n = (int) atoi(coords[coords.size()-1].data());

            pcount = 0;
        }
        else pcount++;
    }

    pos_file.close();
    if (pcount == n)
        return t;
    else
        return tprev;

}

void write_first_nlines(string src, int nlines)
{

    string tmp = src + ".tmp";
    fs::path src_path(src), tmp_path(tmp);

    fs::copy_file(src_path, tmp_path, fs::copy_option::overwrite_if_exists);

    ifstream read_file;
    read_file.open(tmp);

    ofstream write_file;
    write_file.open(src);

    int n = 0;

    string pos_str;
    while(getline(read_file, pos_str))
    {
        if (n >= nlines) break;
        write_file << pos_str <<endl;
        n++;
    }

    read_file.close();
    write_file.close();

}

void write_first_ntsteps(string src, int ntsteps)
{

    string tmp = src + ".tmp";
    fs::path src_path(src), tmp_path(tmp);

    fs::copy_file(src_path, tmp_path, fs::copy_option::overwrite_if_exists);

    ifstream read_file;
    read_file.open(tmp);

    ofstream write_file;
    write_file.open(src);

    string pos_str;
    int nt = 0;
    while(getline(read_file, pos_str))
    {
        if (pos_str[0]=='t'){
            nt++;
            if (nt > ntsteps)
                return;
        }
        write_file << pos_str << endl;
    }
}

void write_first_tsteps(string src, double tstop)
{

    string tmp = src + ".tmp";
    fs::path src_path(src), tmp_path(tmp);

    fs::copy_file(src_path, tmp_path, fs::copy_option::overwrite_if_exists);

    vector<string> coords;

    ifstream read_file;
    read_file.open(tmp);

    ofstream write_file;
    write_file.open(src);

    string pos_str;
    while(getline(read_file, pos_str))
    {
        if (pos_str[0]=='t'){

            boost::split(coords, pos_str, boost::is_any_of("\t "));
            if ( atof(coords[2].data()) >= tstop )
                return;

        }
        write_file << pos_str << endl;
    }
}

template int sgn<int>(int);
template int sgn<double>(double);
template int sgn<float>(float);
