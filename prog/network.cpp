#include "filament_ensemble.h"
#include "motor_ensemble.h"
#include "globals.h"
#include "generate.h"

#include <iostream>
#include <fstream>
#include <iterator>
#include <array>
#include <boost/program_options.hpp>
#include <boost/any.hpp>
#include <typeinfo>

namespace po = boost::program_options;

template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
    return os;
}

ostream &operator<<(ostream &os, const virial_type &a)
{
    fmt::print(os, "{}\t{}\t{}\t{}", a.xx, a.xy, a.yx, a.yy);
    return os;
}

int main(int argc, char **argv)
{
    // BEGIN PROGRAM OPTIONS

    // options allowed only on the command line
    string config_file;
    po::options_description generic("Command Line Only Options");
    generic.add_options()
        ("version,v", "print version string")
        ("help,h", "produce help message")
        ("config,c", po::value<string>(&config_file)->default_value("config/network.cfg"), "name of a configuration file")
        ;

    // environment

    string bnd_cnd;
    double xrange, yrange;

    double dt, tinit, tfinal;
    int nframes, nmsgs;

    double viscosity, temperature;

    string dir;
    int myseed;

    bool restart;
    double restart_time;
    double restart_strain;
    bool steven_continuation_flag;
    int continuation_fr;

    double grid_factor;
    bool quad_off_flag;
    int quad_update_period;

    bool circle_flag; double circle_radius, circle_spring_constant;

    po::options_description config_environment("Environment Options");
    config_environment.add_options()
        ("bnd_cnd,bc", po::value<string>(&bnd_cnd)->default_value("PERIODIC"), "boundary conditions")
        ("xrange", po::value<double>(&xrange)->default_value(10), "size of cell in horizontal direction (um)")
        ("yrange", po::value<double>(&yrange)->default_value(10), "size of cell in vertical direction (um)")

        ("dt", po::value<double>(&dt)->default_value(0.0001), "length of individual timestep in seconds")
        ("tinit", po::value<double>(&tinit)->default_value(0), "time that recording of simulation starts")
        ("tfinal", po::value<double>(&tfinal)->default_value(0.01), "length of simulation in seconds")
        ("nframes", po::value<int>(&nframes)->default_value(1000), "number of times between actin/link/motor positions to are printed to file")
        ("nmsgs", po::value<int>(&nmsgs)->default_value(10000), "number of times simulation progress is printed to stdout")

        ("viscosity", po::value<double>(&viscosity)->default_value(0.001), "Dynamic viscosity to determine friction [mg / (um*s)]. At 20 C, is 0.001 for water")
        ("temperature,temp", po::value<double>(&temperature)->default_value(0.004), "Temp in kT [pN-um] that effects magnituded of Brownian component of simulation")

        ("dir", po::value<string>(&dir)->default_value("."), "output directory")
        ("myseed", po::value<int>(&myseed)->default_value(time(NULL)), "Random number generator myseed")

        // restarts
        ("restart", po::value<bool>(&restart)->default_value(false), "if true, will restart simulation from last timestep recorded")
        ("restart_time", po::value<double>(&restart_time)->default_value(-1), "time to restart simulation from")
        ("restart_strain", po::value<double>(&restart_strain)->default_value(0),"the starting strain for restarting simulation")
        ("steven_continuation_flag", po::value<bool>(&steven_continuation_flag)->default_value(false), "flag to continue from last strain")
        ("continuation_fr", po::value<int>(&continuation_fr)->default_value(0),"the last saved frame from the simulation to be continued")

        // quadrants
        ("grid_factor", po::value<double>(&grid_factor)->default_value(2), "number of grid boxes per um^2")
        ("quad_off_flag", po::value<bool>(&quad_off_flag)->default_value(false), "flag to turn off neighbor list updating")
        ("quad_update_period", po::value<int>(&quad_update_period)->default_value(1), "number of timesteps between actin/link/motor position updates to update quadrants")

        // circular confinement
        ("circle_flag", po::value<bool>(&circle_flag)->default_value(false), "flag to add a circular wall")
        ("circle_radius", po::value<double>(&circle_radius)->default_value(INFINITY), "radius of circular wall")
        ("circle_spring_constant", po::value<double>(&circle_spring_constant)->default_value(0.0), "spring constant of circular wall")
        ;

    // filaments

    string actin_in;
    int npolymer, nmonomer, nmonomer_extra; double extra_bead_prob;

    double actin_length;
    string actin_pos_str;

    double link_length, polymer_bending_modulus, link_stretching_stiffness, fracture_force;
    double rmax, kexv;
    double kgrow, lgrow, l0min, l0max; int nlink_max;

    double occ;

    bool freeze_filaments;

    po::options_description config_actin("Filament Options");
    config_actin.add_options()
        ("actin_in", po::value<string>(&actin_in)->default_value(""), "input actin positions file")
        ("npolymer", po::value<int>(&npolymer)->default_value(3), "number of polymers in the network")
        ("nmonomer", po::value<int>(&nmonomer)->default_value(11), "number of monomers per filament (if extra_bead_prob != 0, then minimum #)")
        ("nmonomer_extra", po::value<int>(&nmonomer_extra)->default_value(0), "max # of monomers per filament")
        ("extra_bead_prob", po::value<double>(&extra_bead_prob)->default_value(0.5), "probability of adding an extra bead from nmonomer...nmonomer_extra")

        ("actin_length", po::value<double>(&actin_length)->default_value(0.5), "Length of a single actin monomer")
        ("actin_pos_str", po::value<string> (&actin_pos_str)->default_value(""), "Starting positions of actin polymers, commas delimit coordinates; semicolons delimit positions")

        ("link_length", po::value<double>(&link_length)->default_value(1), "Length of links connecting monomers")
        ("polymer_bending_modulus", po::value<double>(&polymer_bending_modulus)->default_value(0.068), "Bending modulus of a filament")
        ("fracture_force", po::value<double>(&fracture_force)->default_value(100000000), "pN-- filament breaking point")
        ("link_stretching_stiffness,ks", po::value<double>(&link_stretching_stiffness)->default_value(1), "stiffness of link, pN/um")

        // excluded volume
        ("rmax", po::value<double>(&rmax)->default_value(0.25), "cutoff distance for interactions between actins beads and filaments")
        ("kexv", po::value<double>(&kexv)->default_value(1.0), "parameter of exv force calculation")

        // filament growth
        ("kgrow", po::value<double>(&kgrow)->default_value(0), "rate of filament growth")
        ("lgrow", po::value<double>(&lgrow)->default_value(0), "additional length of filament upon growth")
        ("l0min", po::value<double>(&l0min)->default_value(0), "minimum length a link can shrink to before disappearing")
        ("l0max", po::value<double>(&l0max)->default_value(2), "maximum length a link can grow to before breaking into two links")
        ("nlink_max", po::value<int>(&nlink_max)->default_value(25), "maximum number of links allowed on filament")

        // steric occlusion for binding
        ("occ", po::value<double>(&occ)->default_value(0), "closest distance crosslinkers and motors can bind on a filament")

        ("freeze_filaments", po::value<bool>(&freeze_filaments)->default_value(false), "freeze filaments in place")
        ;

    // motors

    string a_motor_in;
    double a_motor_density;
    string a_motor_pos_str;

    bool motor_intersect_flag;
    double a_linkage_prob;

    double a_m_kon, a_m_kend, a_m_koff, a_m_cut;
    double a_m_kon2, a_m_kend2, a_m_koff2;
    double a_motor_length, a_motor_stiffness;
    double a_motor_v, a_m_stall;
    double a_motor_v2, a_m_stall2;
    double a_m_bend, a_m_ang, a_m_kalign;
    int a_m_align;

    bool dead_head_flag; int dead_head;

    po::options_description config_motors("Motor Options");
    config_motors.add_options()
        ("a_motor_in", po::value<string>(&a_motor_in)->default_value(""), "input motor positions file")
        ("a_motor_density", po::value<double>(&a_motor_density)->default_value(0.05), "number of active motors / um^2")
        ("a_motor_pos_str", po::value<string> (&a_motor_pos_str)->default_value(""), "Starting positions of motors, commas delimit coordinates; semicolons delimit positions")

        ("motor_intersect_flag", po::value<bool>(&motor_intersect_flag)->default_value(false), "flag to put a motor at all filament intersections")
        ("a_linkage_prob", po::value<double>(&a_linkage_prob)->default_value(1), "If motor_intersect_flag, probability that two filaments that intersect will be motor-d")

        // binding/unbinding
        ("a_m_kon", po::value<double>(&a_m_kon)->default_value(1),"active motor on rate")
        ("a_m_koff", po::value<double>(&a_m_koff)->default_value(0.1),"active motor off rate")
        ("a_m_kend", po::value<double>(&a_m_kend)->default_value(0.1),"active motor off rate at filament end")
        ("a_m_cut", po::value<double>(&a_m_cut)->default_value(0.063),"cutoff distance for binding (um)")
        ("a_m_kon2", po::value<double>(&a_m_kon2)->default_value(-1), "active motor on rate when both heads bound")
        ("a_m_koff2", po::value<double>(&a_m_koff2)->default_value(-1), "active motor off rate when both heads bound")
        ("a_m_kend2", po::value<double>(&a_m_kend2)->default_value(-1), "active motor off rate at filament end when both heads bound")

        ("a_motor_length", po::value<double>(&a_motor_length)->default_value(0.4),"active motor rest length (um)")
        ("a_motor_stiffness", po::value<double>(&a_motor_stiffness)->default_value(1),"active motor spring stiffness (pN/um)")
        ("a_m_bend", po::value<double>(&a_m_bend)->default_value(0.04),"bending force constant of active motors (pN) (10kT/um by default)")
        ("a_m_ang", po::value<double>(&a_m_ang)->default_value(pi/2),"equilibrium angle of active motor-filament system")
        ("a_m_align", po::value<int>(&a_m_align)->default_value(0), "active motor alignment type")
        ("a_m_kalign", po::value<double>(&a_m_kalign)->default_value(0.0), "active motor alignment penalty")

        // walking
        ("a_motor_v", po::value<double>(&a_motor_v)->default_value(1),"active motor velocity (um/s)")
        ("a_m_stall", po::value<double>(&a_m_stall)->default_value(0.5),"force beyond which motors don't walk (pN)")
        ("a_motor_v2", po::value<double>(&a_motor_v2)->default_value(NAN),"active motor velocity (um/s) for the second head")
        ("a_m_stall2", po::value<double>(&a_m_stall2)->default_value(NAN),"force beyond which motors don't walk (pN) for the second head")

        ("dead_head_flag", po::value<bool>(&dead_head_flag)->default_value(false), "flag to kill head <dead_head> of all motors")
        ("dead_head", po::value<int>(&dead_head)->default_value(0), "index of head to kill")
        ;

    // crosslinkers

    string p_motor_in;
    double p_motor_density;
    string p_motor_pos_str;

    bool link_intersect_flag; double p_linkage_prob;

    double p_m_kon, p_m_kend, p_m_koff, p_m_cut;
    double p_m_kon2, p_m_kend2, p_m_koff2;
    double p_motor_length, p_motor_stiffness;
    double p_motor_v, p_m_stall;
    double p_motor_v2, p_m_stall2;
    double p_m_bend, p_m_ang, p_m_kalign;
    int p_m_align;

    bool p_dead_head_flag; int p_dead_head;
    bool static_cl_flag;

    po::options_description config_crosslinks("Crosslinker Options");
    config_crosslinks.add_options()
        ("p_motor_in", po::value<string>(&p_motor_in)->default_value(""), "input crosslinker positions file")
        ("p_motor_density", po::value<double>(&p_motor_density)->default_value(0.05), "number of passive motors / um^2")
        ("p_motor_pos_str", po::value<string> (&p_motor_pos_str)->default_value(""), "Starting positions of crosslinks, commas delimit coordinates; semicolons delimit positions")

        ("link_intersect_flag", po::value<bool>(&link_intersect_flag)->default_value(false), "flag to put a cross link at all filament intersections")
        ("p_linkage_prob", po::value<double>(&p_linkage_prob)->default_value(1), "If link_intersect_flag, probability that two filaments that intersect will be linked")

        // binding/unbinding
        ("p_m_kon", po::value<double>(&p_m_kon)->default_value(1),"passive motor on rate")
        ("p_m_koff", po::value<double>(&p_m_koff)->default_value(0.1),"passive motor off rate")
        ("p_m_kend", po::value<double>(&p_m_kend)->default_value(0.1),"passive motor off rate at filament end")
        ("p_m_cut", po::value<double>(&p_m_cut)->default_value(0.063),"cutoff distance for binding (um)")
        ("p_m_kon2", po::value<double>(&p_m_kon2)->default_value(-1), "passive motor on rate when both heads bound")
        ("p_m_koff2", po::value<double>(&p_m_koff2)->default_value(-1), "passive motor off rate when both heads bound")
        ("p_m_kend2", po::value<double>(&p_m_kend2)->default_value(-1), "passive motor off rate at filament end when both heads bound")

        ("p_motor_length", po::value<double>(&p_motor_length)->default_value(0.150),"passive motor rest length (um) (default: filamin)")
        ("p_motor_stiffness", po::value<double>(&p_motor_stiffness)->default_value(1),"passive motor spring stiffness (pN/um)")
        ("p_m_bend", po::value<double>(&p_m_bend)->default_value(0.04),"bending force constant of passive motors (pN) (10kT/um by default)")
        ("p_m_ang", po::value<double>(&p_m_ang)->default_value(pi/2),"equilibrium angle of passive motor-filament system")
        ("p_m_align", po::value<int>(&p_m_align)->default_value(0), "passive motor alignment type")
        ("p_m_kalign", po::value<double>(&p_m_kalign)->default_value(0.0), "passive motor alignment penalty")

        // walking
        ("p_motor_v", po::value<double>(&p_motor_v)->default_value(0),"passive motor velocity (um/s)")
        ("p_m_stall", po::value<double>(&p_m_stall)->default_value(0),"force beyond which xlinks don't walk (pN)")
        ("p_motor_v2", po::value<double>(&p_motor_v2)->default_value(NAN),"passive motor velocity (um/s) for the second head")
        ("p_m_stall2", po::value<double>(&p_m_stall2)->default_value(NAN),"force beyond which xlinks don't walk (pN) for the second head")

        ("p_dead_head_flag", po::value<bool>(&p_dead_head_flag)->default_value(false), "flag to kill head <dead_head> of all crosslinks")
        ("p_dead_head", po::value<int>(&p_dead_head)->default_value(0), "index of head to kill")

        ("static_cl_flag", po::value<bool>(&static_cl_flag)->default_value(false), "flag to indicate compeletely static xlinks; i.e, no walking, no detachment")
        ;

    int n_bw_shear;

    double d_strain_freq, d_strain_pct;
    double time_of_dstrain, time_of_dstrain2;

    bool diff_strain_flag, osc_strain_flag;
    bool stress_flag; double stress1, stress_rate1, stress2, stress_rate2;

    bool shear_motor_flag;

    po::options_description config_shear("Shear Options");
    config_shear.add_options()
        ("n_bw_shear", po::value<int>(&n_bw_shear)->default_value(1000000000), "number of timesteps between subsequent shears")

        ("d_strain_freq", po::value<double>(&d_strain_freq)->default_value(1), "differential strain frequency")
        ("d_strain_pct", po::value<double>(&d_strain_pct)->default_value(0), "differential strain amplitude")
        ("time_of_dstrain", po::value<double>(&time_of_dstrain)->default_value(10000), "time when differential strain starts")
        ("time_of_dstrain2", po::value<double>(&time_of_dstrain2)->default_value(10000), "time when second differential strain starts")

        ("diff_strain_flag", po::value<bool>(&diff_strain_flag)->default_value(false), "flag to turn on linear differential strain")
        ("osc_strain_flag", po::value<bool>(&osc_strain_flag)->default_value(false), "flag to turn on oscillatory differential strain")

        // stress
        ("stress_flag", po::value<bool>(&stress_flag)->default_value(false), "flag to turn on constant stress")
        ("stress", po::value<double>(&stress1)->default_value(0.0), "value of constant stress")
        ("stress_rate", po::value<double>(&stress_rate1)->default_value(0.0), "decay rate to specified value of stress, in weird units")
        ("stress2", po::value<double>(&stress2)->default_value(0.0), "second value of constant stress")
        ("stress_rate2", po::value<double>(&stress_rate2)->default_value(0.0), "second decay rate to specified value of stress, in weird units")

        ("shear_motor_flag", po::value<bool>(&shear_motor_flag)->default_value(false), "flag to turn on shearing for motors")
        ;

    po::options_description config;
    config.add(config_environment).add(config_actin).add(config_motors).add(config_crosslinks).add(config_shear);

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << generic << "\n";
        std::cout << config << "\n";
        return 1;
    }

    ifstream ifs(config_file);
    if (!ifs){
        cout<<"can not open config file: "<<config_file<<"\n";
        return 0;
    } else {
        po::store(po::parse_config_file(ifs, config), vm);
        po::notify(vm);
    }

    string tdir   = dir  + "/txt_stack";
    string ddir   = dir  + "/data";

    string afile  = tdir + "/actins.txt";
    string lfile  = tdir + "/links.txt";
    string amfile = tdir + "/amotors.txt";
    string pmfile = tdir + "/pmotors.txt";
    string thfile = ddir + "/filament_e.txt";
    string pefile = ddir + "/pe.txt";

    if (fs::create_directory(fs::path(dir))) cerr << "Directory Created: " << dir << endl;
    if (fs::create_directory(fs::path(tdir))) cerr << "Directory Created: " << tdir << endl;
    if (fs::create_directory(fs::path(ddir))) cerr << "Directory Created: " << ddir << endl;

    // Write the full configuration file
    {
        ofstream o_file(dir + "/data/config_full.cfg");
        for (auto const &it : vm) {
            if (it.first == "config") continue;
            boost::any val = it.second.value();

            if(typeid(bool) == val.type())
                o_file << it.first <<"="<< boost::any_cast<bool>(val) <<endl;
            else if(typeid(int) == val.type())
                o_file << it.first <<"="<< boost::any_cast<int>(val) <<endl;
            else if(typeid(double) == val.type())
                o_file << it.first <<"="<< boost::any_cast<double>(val) <<endl;
            else if(typeid(string) == val.type())
                o_file << it.first <<"="<< boost::any_cast<string>(val) <<endl;
        }
    }

    if (a_m_kon2 == -1) a_m_kon2 = a_m_kon;
    if (a_m_koff2 == -1) a_m_koff2 = a_m_koff;
    if (a_m_kend2 == -1) a_m_kend2 = a_m_kend;

    if (p_m_kon2 == -1) p_m_kon2 = p_m_kon;
    if (p_m_koff2 == -1) p_m_koff2 = p_m_koff;
    if (p_m_kend2 == -1) p_m_kend2 = p_m_kend;

    // END PROGRAM OPTIONS

    int n_bw_stdout = max(int((tfinal)/(dt*double(nmsgs))),1);
    int n_bw_print  = max(int((tfinal)/(dt*double(nframes))),1);
    int unprinted_count = int(double(tinit)/dt);

    vector<vector<double>> actin_pos_vec;
    vector<vector<double>> a_motor_pos_vec, p_motor_pos_vec;

    // BEGIN READ CONFIGURATIONS

    // To get positions from input files:
    if (actin_in.size() > 0)
        actin_pos_vec   = file2vecvec(actin_in, "\t");
    if (a_motor_in.size() > 0)
        a_motor_pos_vec = file2vecvec(a_motor_in, "\t");
    if (p_motor_in.size() > 0)
        p_motor_pos_vec = file2vecvec(p_motor_in, "\t");

    // To restart a whole trajectory from it's last full timestep :
    if (restart){

        double tf_prev  = min(last_full_timestep(afile), last_full_timestep(lfile));
        if (a_motor_density > 0)
            tf_prev = min(tf_prev, last_full_timestep(amfile));
        if (p_motor_density > 0)
            tf_prev = min(tf_prev, last_full_timestep(pmfile));

        if (restart_time == -1 || restart_time > tf_prev)
            restart_time = tf_prev;

        cout<<"\nRestarting from t = "<<restart_time<<endl;

        double nprinted = restart_time / (dt*n_bw_print);

        actin_pos_vec   = traj2vecvec(afile, "\t ", restart_time);
        a_motor_pos_vec = traj2vecvec(amfile, "\t ", restart_time);
        p_motor_pos_vec = traj2vecvec(pmfile, "\t ", restart_time);

        // for actins, links, amotors, pmotors:
        // do:
        //   copy whole file into temp
        //   while hasn't reached tf in temp file:
        //      write from copy into afile
        write_first_tsteps(afile,  restart_time);
        write_first_tsteps(lfile,  restart_time);
        write_first_tsteps(amfile, restart_time);
        write_first_tsteps(pmfile, restart_time);

        write_first_tsteps(thfile, restart_time);
        write_first_nlines(pefile, (int) nprinted);

        // additional outputs below

        tinit = restart_time;

        if (steven_continuation_flag) {
            restart_strain = get_restart_strain(pefile, continuation_fr);
        }

    }

    // END READ CONFIGURATIONS

    // compute derived quantities

    if (polymer_bending_modulus < 0){ //This is a flag for using the temperature for the bending modulus
        polymer_bending_modulus = 10*temperature; // 10um * kT
    }

    double actin_density = double(npolymer*nmonomer)/(xrange*yrange);//0.65;
    cout<<"\nDEBUG: actin_density = "<<actin_density;
    double link_bending_stiffness    = polymer_bending_modulus / link_length;

    // set number of quadrants to 1 if there are no crosslinkers/motors
    int xgrid, ygrid;
    if(a_motor_density == 0 && a_motor_pos_vec.size() == 0 &&
            p_motor_density==0 && p_motor_pos_vec.size() == 0 &&
            !link_intersect_flag && !motor_intersect_flag){
        xgrid = 1;
        ygrid = 1;
    }
    else{
        xgrid  = (int) round(grid_factor*xrange);
        ygrid  = (int) round(grid_factor*yrange);
    }
    double d_strain_amp = d_strain_pct * xrange;

    box *bc = new box(bnd_cnd, xrange, yrange, restart_strain);

    set_seed(myseed);

    // BEGIN GENERATE CONFIGURATIONS

    if (actin_pos_vec.size() == 0 && actin_in.size() == 0) {
        vector<array<double, 3>> actin_position_arrs;
        if (actin_pos_str.size() > 0) {
            // read positions from input strings in config file
            actin_position_arrs = str2arrvec(actin_pos_str, ":", ",");
        }
        actin_pos_vec = generate_filament_ensemble(
                bc, npolymer, nmonomer, nmonomer_extra, extra_bead_prob,
                dt, temperature, actin_length, link_length,
                actin_position_arrs, link_bending_stiffness, myseed);
    }

    if (link_intersect_flag)
        p_motor_pos_vec = spring_spring_intersections(bc, actin_pos_vec, p_motor_length, p_linkage_prob);

    if (motor_intersect_flag)
        a_motor_pos_vec = spring_spring_intersections(bc, actin_pos_vec, a_motor_length, a_linkage_prob);

    if (a_motor_pos_vec.size() == 0 && a_motor_in.size() == 0) {
        vector<array<double, 3>> a_motor_position_arrs;
        if (a_motor_pos_str.size() > 0) {
            a_motor_position_arrs = str2arrvec(a_motor_pos_str, ":", ",");
        }
        a_motor_pos_vec = generate_motor_ensemble(bc, a_motor_density, a_motor_length, a_motor_position_arrs);
    }

    if (p_motor_pos_vec.size() == 0 && p_motor_in.size() == 0) {
        vector<array<double, 3>> p_motor_position_arrs;
        if (p_motor_pos_str.size() > 0) {
            p_motor_position_arrs = str2arrvec(p_motor_pos_str, ":", ",");
        }
        p_motor_pos_vec = generate_motor_ensemble(bc, p_motor_density, p_motor_length, p_motor_position_arrs);
    }

    // END GENERATE CONFIGURATIONS

    // BEGIN CREATE NETWORK OBJECTS

    cout<<"\nCreating actin network..";
    filament_ensemble *net = new filament_ensemble(
            bc, actin_pos_vec, {xgrid, ygrid}, dt,
            temperature, viscosity, link_length,
            link_stretching_stiffness, link_bending_stiffness,
            fracture_force, rmax, kexv);

    // additional options
    net->set_growing(kgrow, lgrow, l0min, l0max, nlink_max);
    if (quad_off_flag) net->get_quads()->use_quad(false);

    cout<<"\nAdding active motors...";
    motor_ensemble *myosins = new motor_ensemble(
            a_motor_pos_vec, dt, temperature,
            a_motor_length, net, a_motor_v, a_motor_stiffness, a_m_kon, a_m_koff,
            a_m_kend, a_m_stall, a_m_cut, viscosity);

    myosins->set_binding_two(a_m_kon2, a_m_koff2, a_m_kend2);
    myosins->set_bending(a_m_bend, a_m_ang);
    if (dead_head_flag) myosins->kill_heads(dead_head);
    if (shear_motor_flag) myosins->use_shear(true);
    if (a_m_kalign != 0.0) {
        if (a_m_align == 1) myosins->set_par(a_m_kalign);
        if (a_m_align == 0) myosins->set_align(a_m_kalign);
        if (a_m_align == -1) myosins->set_antipar(-a_m_kalign);
    }
    if (!std::isnan(a_motor_v2)) myosins->set_velocity(a_motor_v, a_motor_v2);
    if (!std::isnan(a_m_stall2)) myosins->set_stall_force(a_m_stall, a_m_stall2);

    cout<<"Adding passive motors (crosslinkers) ...\n";
    motor_ensemble *crosslks = new motor_ensemble(
            p_motor_pos_vec, dt, temperature,
            p_motor_length, net, p_motor_v, p_motor_stiffness, p_m_kon, p_m_koff,
            p_m_kend, p_m_stall, p_m_cut, viscosity);

    crosslks->set_binding_two(p_m_kon2, p_m_koff2, p_m_kend2);
    crosslks->set_bending(p_m_bend, p_m_ang);
    if (p_dead_head_flag) crosslks->kill_heads(p_dead_head);
    if (shear_motor_flag) crosslks->use_shear(true);
    if (static_cl_flag) crosslks->use_static(true);
    if (p_m_kalign != 0.0) {
        if (p_m_align == 1) crosslks->set_par(p_m_kalign);
        if (p_m_align == 0) crosslks->set_align(p_m_kalign);
        if (p_m_align == -1) crosslks->set_antipar(-p_m_kalign);
    }
    if (!std::isnan(p_motor_v2)) crosslks->set_velocity(p_motor_v, p_motor_v2);
    if (!std::isnan(p_m_stall2)) crosslks->set_stall_force(p_m_stall, p_m_stall2);

    if (circle_flag) {
        net->set_external(new ext_circle(circle_spring_constant, circle_radius));
        myosins->set_external(new ext_circle(circle_spring_constant, circle_radius));
        crosslks->set_external(new ext_circle(circle_spring_constant, circle_radius));
    }

    if (occ > 0.0) {
        myosins->set_occ(occ);
        crosslks->set_occ(occ);
    }

    // compute forces and energies
    net->compute_forces();
    crosslks->compute_forces();
    myosins->compute_forces();

    // END CREATE NETWORK OBJECTS

    // run simulation

    cout<<"\nUpdating motors, filaments and crosslinks in the network..";
    string time_str;

    // open output files
    // for restarts, append instead of writing from the start
    // files are automatically closed by RAII
    ios_base::openmode write_mode = (restart) ? ios_base::app : ios_base::out;
    ofstream file_a(afile, write_mode);
    ofstream file_l(lfile, write_mode);
    ofstream file_am(amfile, write_mode);
    ofstream file_pm(pmfile, write_mode);
    ofstream file_th(thfile, write_mode);
    ofstream file_pe(pefile, write_mode);

    // set up occ
    size_t n_myosins = myosins->get_nmotors();
    size_t n_crosslks = crosslks->get_nmotors();
    vector<size_t> motor_ix(n_myosins + n_crosslks);
    for (size_t i = 0; i < n_myosins + n_crosslks; i++) {
        motor_ix[i] = i;
    }

    int count; double t;
    for (count = 0, t = tinit; t <= tfinal; count++, t += dt) {

        // output to file
        if (t+dt/100 >= tinit && (count-unprinted_count)%n_bw_print==0) {

            if (t>tinit) time_str ="\n";
            time_str += "t = "+to_string(t);

            fmt::print(file_a, "{}\tN = {}", time_str, net->get_nbeads());
            net->write_beads(file_a);

            fmt::print(file_l, "{}\tN = {}", time_str, net->get_nsprings());
            net->write_springs(file_l);

            fmt::print(file_am, "{}\tN = {}", time_str, myosins->get_nmotors());
            myosins->motor_write(file_am);

            fmt::print(file_pm, "{}\tN = {}", time_str, crosslks->get_nmotors());
            crosslks->motor_write(file_pm);

            fmt::print(file_th, "{}\tN = {}", time_str, net->get_nfilaments());
            net->write_thermo(file_th);

            fmt::print(file_pe,
                    "{}\t{}\t{}\t{}\t"

                    "{}\t{}\t{}\t{}\t"
                    "{}\t{}\t{}\t{}\t"
                    "{}\t{}\t{}\t{}\t"

                    "{}\t{}\t{}\t{}\t"
                    "{}\t{}\t{}\t{}\t"
                    "{}\t{}\t{}\t{}\n",

                    t, bc->get_xbox(), bc->get_ybox(), bc->get_delrx(),

                    net->get_stretching_energy(),
                    net->get_bending_energy(),
                    net->get_excluded_energy(),
                    net->get_external_energy(),

                    myosins->get_stretching_energy(),
                    myosins->get_bending_energy(),
                    myosins->get_alignment_energy(),
                    myosins->get_external_energy(),

                    crosslks->get_stretching_energy(),
                    crosslks->get_bending_energy(),
                    crosslks->get_alignment_energy(),
                    crosslks->get_external_energy(),

                    net->get_stretching_virial(),
                    net->get_bending_virial(),
                    net->get_excluded_virial(),
                    net->get_external_virial(),

                    myosins->get_stretching_virial(),
                    myosins->get_bending_virial(),
                    myosins->get_alignment_virial(),
                    myosins->get_external_virial(),

                    crosslks->get_stretching_virial(),
                    crosslks->get_bending_virial(),
                    crosslks->get_alignment_virial(),
                    crosslks->get_external_virial());

            file_a<<std::flush;
            file_l<<std::flush;
            file_am<<std::flush;
            file_pm<<std::flush;
            file_th<<std::flush;
            file_pe<<std::flush;
        }

        // print to stdout
        if (count%n_bw_stdout==0) {
            fmt::print("\nCount: {}\tTime: {} s\tShear: {} um", count, t, bc->get_delrx());
            //net->print_filament_thermo();
            net->print_network_thermo();
            crosslks->print_ensemble_thermo();
            myosins->print_ensemble_thermo();
        }

        // shear
        if (t >= time_of_dstrain && count % n_bw_shear == 0) {
            double d_strain = 0.0;
            if (stress_flag) {
                double stress, stress_rate;
                if (t < time_of_dstrain2) {
                    stress = stress1;
                    stress_rate = stress_rate1;
                } else {
                    stress = stress2;
                    stress_rate = stress_rate2;
                }

                virial_type virial
                    = net->get_potential_virial()
                    + myosins->get_potential_virial()
                    + crosslks->get_potential_virial();

                double xbox = bc->get_xbox();
                double ybox = bc->get_ybox();
                double delrx = bc->get_delrx();
                double f_delrx = -2.0 * virial.yx / ybox;
                // if strain = delrx/ybox, then
                // strain += stress_rate (-dU/dstrain / area + stress) dt
                // units:
                // - stress: energy / area
                // - stress_rate: area / (energy time)
                d_strain = delrx + stress_rate * (stress + f_delrx / xbox) * ybox * dt;
            } else if (osc_strain_flag) {
                //d_strain += d_strain_amp * sin(2*pi*d_strain_freq * ( t - time_of_dstrain) );
                d_strain = restart_strain + d_strain_amp*4*d_strain_freq*((t-time_of_dstrain) - 1/(d_strain_freq*2)*floor(2*(t-time_of_dstrain)*d_strain_freq + 0.5))*pow(-1,floor(2*(t-time_of_dstrain)*d_strain_freq + 0.5));
            } else if (diff_strain_flag) {
                d_strain = restart_strain + d_strain_amp*d_strain_freq*(t - time_of_dstrain);
            }
            bc->update_d_strain(d_strain - bc->get_delrx());
        }

        // Brownian dynamics and motor walking
        if (!freeze_filaments)
            net->integrate();
        crosslks->integrate();
        myosins->integrate();

        // filament growth and fracturing
        // also unbinds motors
        if (!freeze_filaments)
            net->montecarlo();

        if (quad_off_flag) {
            // we want results that are correct regardless of other settings when quadrants are off
            // this just builds a list of all springs, which are then handed to attachment/etc
            net->quad_update_serial();

        } else if (count % quad_update_period == 0) {
            // when quadrants are on, this actually builds quadrants
            net->quad_update_serial();

        }

        // motor attachment/detachment
        if (occ > 0.0) {
            std::shuffle(motor_ix.begin(), motor_ix.end(), get_rng());
            for (size_t i : motor_ix) {
                if (i < n_myosins) {
                    myosins->try_attach_detach(i);
                } else {
                    crosslks->try_attach_detach(i - n_myosins);
                }
            }
        } else {
            crosslks->try_attach_detach();
            myosins->try_attach_detach();
        }

        // compute forces and energies
        if (!freeze_filaments)
            net->compute_forces();
        crosslks->compute_forces();
        myosins->compute_forces();
    }

    file_a << "\n";
    file_l << "\n";
    file_am << "\n";
    file_pm << "\n";
    file_th << "\n";

    //Delete all objects created
    cout<<"\nHere's where I think I delete things\n";

    delete myosins;
    delete crosslks;
    delete net;
    delete bc;

    cout<<"\nTime counts: "<<count;
    cout<<"\nExecuted";
    cout<<"\n Done\n";

    return 0;
}
