#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "box.h"
#include "quadrants.h"
#include "generate.h"
#include "potentials.h"
#include "ext.h"
#include "filament.h"
#include "filament_ensemble.h"
#include "motor.h"
#include "motor_ensemble.h"

namespace py = pybind11;

void vec_bindings(py::module &m)
{
    py::class_<vec_type>(m, "Vector")
        .def(py::init<double, double>())
        .def_readwrite("x", &vec_type::x)
        .def_readwrite("y", &vec_type::y)
        ;

    py::class_<virial_type>(m, "Virial")
        .def(py::init<double, double, double, double>())
        .def_readwrite("xx", &virial_type::xx)
        .def_readwrite("xy", &virial_type::xy)
        .def_readwrite("yx", &virial_type::yx)
        .def_readwrite("yy", &virial_type::yy)
        ;
}

void globals_bindings(py::module &m)
{
    // random
    m.def("rng_u", &rng_u);
    m.def("rng_n", &rng_n);
    m.def("gen_seed", &gen_seed);
    m.def("get_seed", &get_seed);
    m.def("set_seed", &set_seed);

    // misc
    py::class_<mc_prob>(m, "MCProb")
        .def(py::init<>())
        .def("__call__", &mc_prob::operator())
        .def_readwrite("prob", &mc_prob::prob)
        .def_readwrite("used", &mc_prob::used)
        ;
    py::class_<fp_index_type>(m, "FPIndex")
        .def(py::init<>())
        .def_readwrite("f_index", &fp_index_type::f_index)
        .def_readwrite("p_index", &fp_index_type::p_index)
        ;
}

void box_bindings(py::module &m)
{
    py::enum_<bc_type>(m, "BCType")
        .value("lees_edwards", bc_type::lees_edwards)
        .value("periodic", bc_type::periodic)
        .value("xperiodic", bc_type::xperiodic)
        .value("nonperiodic", bc_type::nonperiodic)
        ;

    py::class_<box>(m, "Box")
        // constructors
        .def(py::init<string, double, double, double>())
        // get/set state
        .def_property_readonly("bc", &box::get_BC)
        .def_property_readonly("xbox", &box::get_xbox)
        .def_property_readonly("ybox", &box::get_ybox)
        .def_property_readonly("delrx", &box::get_delrx)
        .def("add_delrx", &box::update_d_strain)
        // methods
        .def("rij", &box::rij_bc)
        .def("pos", &box::pos_bc)
        .def("dist", &box::dist_bc)
        .def("dot", &box::dot_bc)
        ;

    m.def("seg_seg_intersection", &seg_seg_intersection_bc);
}

void quadrants_bindings(py::module &m)
{
    py::class_<quadrants>(m, "Quadrants")
        .def(py::init<box *, array<int, 2>>())
        .def("use_quad", &quadrants::use_quad)
        .def("use_all", &quadrants::use_all)
        .def("add_spring", &quadrants::add_spring)
        .def("attach_list", &quadrants::get_attach_list)
        .def("build_pairs", &quadrants::build_pairs)
        .def("pairs", &quadrants::get_pairs)
        .def("clear", &quadrants::clear)
        .def_property_readonly("nq", &quadrants::get_nq)
        .def("quad", &quadrants::get_quad)
        .def_property_readonly("cut", &quadrants::get_cut)
        .def_property_readonly("paircut", &quadrants::get_paircut)
        .def("check_duplicates", &quadrants::check_duplicates)
        .def("check_quad", &quadrants::check_quad)
        .def("check_pairs", &quadrants::check_pairs)
        ;
}

void generate_bindings(py::module &m)
{
    m.def("generate_filaments",
            static_cast<vector<vector<double>>(*)(
                box *,
                int, int, int, double,
                double, double, double,
                vector<array<double, 3>>,
                double)>(&generate_filament_ensemble));
    m.def("generate_filaments",
            static_cast<vector<vector<double>>(*)(
                box *,
                double, double, double, int, double,
                vector<array<double, 3>>,
                double)>(&generate_filament_ensemble));
    m.def("generate_motors", &generate_motor_ensemble);
    m.def("spring_spring_intersections", &spring_spring_intersections);
}

void potentials_bindings(py::module &m)
{
    py::class_<stretch_fene>(m, "StretchFene")
        .def(py::init<double, double, double>())
        .def("tension", &stretch_fene::tension)
        .def("energy", &stretch_fene::energy)
        .def_readwrite("kl", &stretch_fene::kl)
        .def_readwrite("l0", &stretch_fene::l0)
        .def_readwrite("max_ext", &stretch_fene::max_ext)
        .def_readwrite("eps_ext", &stretch_fene::eps_ext)
        ;

    py::class_<bend_result_type>(m, "BendResult")
        .def_readwrite("energy", &bend_result_type::energy)
        .def_readwrite("force1", &bend_result_type::force1)
        .def_readwrite("force2", &bend_result_type::force2)
        ;

    m.def("bend_cos_angle", &bend_angle);
    m.def("cos_angle", &angle);

    m.def("bend_harmonic", &bend_harmonic);
    m.def("bend_harmonic_energy", &bend_harmonic_energy);
}

void ext_bindings(py::module &m)
{
    py::class_<ext_result_type>(m, "ExternalResult")
        .def(py::init<>())
        .def_readwrite("energy", &ext_result_type::energy)
        .def_readwrite("force", &ext_result_type::force)
        .def_readwrite("virial", &ext_result_type::virial)
        ;

    py::class_<external>(m, "External")
        ;

    py::class_<ext_circle, external>(m, "ExternalCircle")
        .def(py::init<double, double>())
        .def("compute", &ext_circle::compute)
        ;

    py::class_<ext_circle_fene, external>(m, "ExternalCircleFene")
        .def(py::init<double, double, double>())
        .def("compute", &ext_circle_fene::compute)
        ;

    py::class_<ext_rectangle, external>(m, "ExternalRectangle")
        .def(py::init<double, double, double>())
        .def("compute", &ext_rectangle::compute)
        ;

    py::class_<ext_xperiodic, external>(m, "ExternalXPeriodic")
        .def(py::init<double, double>())
        .def("compute", &ext_xperiodic::compute)
        ;
}

void filament_bindings(py::module &m)
{
    py::class_<filament>(m, "Filament")
        .def(py::init<
                box *,
                vector<vec_type>, vector<vec_type>,
                double, double,
                double, double, double,
                double, double, double>())
        .def("add_bead", &filament::add_bead)
        // get/set bead state
        .def_property_readonly("num_beads", &filament::get_nbeads)
        .def("position", &filament::get_bead_position)
        .def("force", &filament::get_force)
        // spring methods
        .def_property_readonly("num_springs", &filament::get_nsprings)
        .def("start", &filament::get_start)
        .def("end", &filament::get_end)
        .def("length", &filament::get_llength)
        .def("displacement", &filament::get_disp)
        .def("direction", &filament::get_direction)
        .def("tension", &filament::get_tension)
        .def("stretching_energy",
                [](filament *f, int i) { return f->get_stretching_energy(i); })
        .def("intpoint", &filament::intpoint)
        .def("kl", &filament::get_kl)
        .def("l0", &filament::get_l0)
        // get thermo
        .def_property_readonly("pe_stretch",
                [](filament *f) { return f->get_stretching_energy(); })
        .def_property_readonly("pe_bend", &filament::get_bending_energy)
        .def_property_readonly("vir_stretch", &filament::get_stretching_virial)
        .def_property_readonly("vir_bend", &filament::get_bending_virial)
        // output state
        .def("output_beads", &filament::output_beads)
        .def("output_springs", &filament::output_springs)
        // output thermo
        .def("output_thermo", &filament::output_thermo)
        // methods
        .def("add_delrx", &filament::update_d_strain)
        .def("integrate", &filament::update_positions)
        .def("add_stretch_forces", &filament::update_stretching)
        .def("add_bend_forces", &filament::update_bending)
        .def("add_forces", &filament::update_forces)
        // attachment
        .def("new_attached", &filament::new_attached)
        .def("del_attached", &filament::del_attached)
        .def("get_attached_l", &filament::get_attached_l)
        .def("get_attached_position", &filament::get_attached_pos)
        .def("add_attached_force", &filament::add_attached_force)
        .def("add_attached_position", &filament::add_attached_pos)
        .def("at_barbed_end", &filament::at_barbed_end)
        .def("at_pointed_end", &filament::at_pointed_end)
        // fracturing
        .def("try_fracture", &filament::try_fracture)
        .def("fracture", &filament::fracture)
        .def("detach_all_motors", &filament::detach_all_motors)
        // growing
        .def("set_max_l0", &filament::set_l0_max)
        .def("set_max_springs", &filament::set_nsprings_max)
        .def("set_min_l0", &filament::set_l0_min)
        .def("set_grow_rate", &filament::set_kgrow)
        .def("set_grow_length", &filament::set_lgrow)
        .def("try_grow", &filament::update_length)
        .def("grow", &filament::grow)
        // pulling
        .def_property_readonly("end_to_end", &filament::get_end2end)
        .def("pull_on_ends", &filament::pull_on_ends)
        .def("affine_pull", &filament::affine_pull)
        ;
}

void filament_ensemble_bindings(py::module &m)
{
    py::class_<filament_ensemble>(m, "FilamentEnsemble")
        .def(py::init<
                box *, vector<vector<double>>, array<int, 2>,
                double, double, double,
                double, double,
                double, double,
                double, double, double>())
        .def("set_external", &filament_ensemble::set_external)
        // quadrants
        .def_property_readonly("quads", &filament_ensemble::get_quads)
        .def("update_quads", &filament_ensemble::quad_update_serial)
        .def("attach_list", &filament_ensemble::get_attach_list)
        // state
        .def_property_readonly("num_filaments", &filament_ensemble::get_nfilaments)
        .def_property_readonly("box", &filament_ensemble::get_box)
        .def("intersect", &filament_ensemble::intersect)
        // bead methods
        .def_property_readonly("num_beads",
                [](filament_ensemble *f) { return f->get_nbeads(); })
        .def("num_filament_beads",
                [](filament_ensemble *f, int i) { return f->get_nbeads(i); })
        .def("position", &filament_ensemble::get_pos)
        .def("force", &filament_ensemble::get_force)
        .def("is_polymer_start", &filament_ensemble::is_polymer_start)
        // spring methods
        .def("num_springs",
                [](filament_ensemble *f) { return f->get_nsprings(); })
        .def("num_filament_springs",
                [](filament_ensemble *f, int i) { return f->get_nsprings(i); })
        .def("displacement", &filament_ensemble::get_disp)
        .def("direction", &filament_ensemble::get_direction)
        .def("start", &filament_ensemble::get_start)
        .def("end", &filament_ensemble::get_end)
        .def("intpoint", &filament_ensemble::intpoint)
        .def("length", &filament_ensemble::get_llength)
        // attached methods
        .def("new_attached", &filament_ensemble::new_attached)
        .def("del_attached", &filament_ensemble::del_attached)
        .def("get_attached_fl", &filament_ensemble::get_attached_fl)
        .def("get_attached_position", &filament_ensemble::get_attached_pos)
        .def("add_attached_force", &filament_ensemble::add_attached_force)
        .def("add_attached_position", &filament_ensemble::add_attached_pos)
        .def("get_attached_direction", &filament_ensemble::get_attached_direction)
        .def("get_attached_displacement", &filament_ensemble::get_attached_disp)
        .def("add_attached_displacement_force", &filament_ensemble::add_attached_disp_force)
        .def("at_pointed_end", &filament_ensemble::at_pointed_end)
        .def("at_barbed_end", &filament_ensemble::at_barbed_end)
        // get thermo
        .def_property_readonly("pe", &filament_ensemble::get_potential_energy)
        .def_property_readonly("pe_stretch", &filament_ensemble::get_stretching_energy)
        .def_property_readonly("pe_bend", &filament_ensemble::get_bending_energy)
        .def_property_readonly("pe_excluded", &filament_ensemble::get_excluded_energy)
        .def_property_readonly("pe_external", &filament_ensemble::get_external_energy)
        .def_property_readonly("vir", &filament_ensemble::get_potential_virial)
        .def_property_readonly("vir_stretch", &filament_ensemble::get_stretching_virial)
        .def_property_readonly("vir_bend", &filament_ensemble::get_bending_virial)
        .def_property_readonly("vir_excluded", &filament_ensemble::get_excluded_virial)
        .def_property_readonly("vir_external", &filament_ensemble::get_external_virial)
        // dynamics
        .def("integrate", &filament_ensemble::integrate)
        .def("add_delrx", &filament_ensemble::update_d_strain)
        .def("set_growing", &filament_ensemble::set_growing)
        .def("try_grow", &filament_ensemble::try_grow)
        .def("try_fracture", &filament_ensemble::try_fracture)
        // update forces/energies
        .def("compute_forces", &filament_ensemble::compute_forces)
        .def("add_stretch_forces", &filament_ensemble::update_stretching)
        .def("add_bend_forces", &filament_ensemble::update_bending)
        .def("add_excluded_forces", &filament_ensemble::update_excluded_volume)
        .def("add_external_forces", &filament_ensemble::update_external)
        .def("add_forces", &filament_ensemble::update_forces)
        .def("update_energies", &filament_ensemble::update_energies)
        // output
        .def("output_beads", &filament_ensemble::output_beads)
        .def("output_springs", &filament_ensemble::output_springs)
        .def("output_thermo", &filament_ensemble::output_thermo)
        ;
}

void motor_bindings(py::module &m)
{
    py::enum_<motor_state>(m, "MotorState")
        .value("free", motor_state::free)
        .value("bound", motor_state::bound)
        .value("dead", motor_state::dead)
        .value("inactive", motor_state::inactive)
        ;

    py::class_<motor>(m, "Motor")
        // constructors
        .def(py::init<
                vector<double>,
                double, filament_ensemble *,
                double, double, double, double,
                double, double, double,
                double, double,
                double>())
        // settings
        .def("set_binding_two", &motor::set_binding_two)
        .def("set_bending", &motor::set_bending)
        .def_property_readonly("kb", &motor::get_kb)
        .def_property_readonly("ang0", &motor::get_th0)
        .def("set_parallel_alignment", &motor::set_par)
        .def("set_alignment", &motor::set_align)
        .def("set_antiparallel_alignment", &motor::set_antipar)
        .def_property_readonly("alignment_type", &motor::get_align)
        .def_property_readonly("kalign", &motor::get_kalign)
        .def("set_external", &motor::set_external)
        // manipulate state
        .def("relax_head", &motor::relax_head)
        .def("kill_head", &motor::kill_head)
        .def("revive_head", &motor::revive_head)
        .def("deactivate_head", &motor::deactivate_head)
        // get state
        .def_property_readonly("state", &motor::get_states)
        .def_property_readonly("start", &motor::get_h0)
        .def_property_readonly("end", &motor::get_h1)
        .def_property_readonly("f_index", &motor::get_f_index)
        .def_property_readonly("l_index", &motor::get_l_index)
        // get forces
        .def_property_readonly("forces", &motor::get_force)
        .def_property_readonly("stretch_forces", &motor::get_s_force)
        .def_property_readonly("bend_forces", &motor::get_b_force)
        .def_property_readonly("external_forces", &motor::get_ext_force)
        .def_property_readonly("projected_forces", &motor::get_force_proj)
        // update forces
        .def("clear_forces", &motor::clear_forces)
        .def("update_forces", &motor::update_force)
        .def("add_bend_forces", &motor::update_bending)
        .def("add_align_forces", &motor::update_alignment)
        .def("add_external_forces", &motor::update_external)
        .def("update_projected_forces", &motor::update_force_proj)
        .def("apply_forces_to_filaments", &motor::filament_update)
        // dynamics
        .def("add_delrx", &motor::update_d_strain)
        .def("brownian_relax", &motor::brownian_relax)
        .def("walk", &motor::walk)
        .def("step", &motor::step)
        // monte carlo
        .def("metropolis_prob", &motor::metropolis_prob)
        .def("alignment_penalty", &motor::alignment_penalty)
        // attach
        .def("try_attach", &motor::try_attach)
        .def("allowed_bind", &motor::allowed_bind)
        .def("attach_head", &motor::attach_head)
        // detach
        .def("try_detach", &motor::try_detach)
        .def("generate_off_pos", &motor::generate_off_pos)
        .def("detach_head_without_moving", &motor::detach_head_without_moving)
        .def("detach_head", &motor::detach_head)
        // get thermo
        .def_property_readonly("pe_stretch", &motor::get_stretching_energy)
        .def_property_readonly("pe_bend", &motor::get_bending_energy)
        .def_property_readonly("pe_align", &motor::get_alignment_energy)
        .def_property_readonly("pe_external", &motor::get_external_energy)
        .def_property_readonly("vir_stretch", &motor::get_stretching_virial)
        .def_property_readonly("vir_bend", &motor::get_bending_virial)
        .def_property_readonly("vir_align", &motor::get_alignment_virial)
        .def_property_readonly("vir_external", &motor::get_external_virial)
        // output state
        .def("output", &motor::output)
        ;
}

void motor_ensemble_bindings(py::module &m)
{
    py::class_<motor_ensemble>(m, "MotorEnsemble")
        .def(py::init<
                vector<vector<double>>,double, double,
                double, filament_ensemble *, double, double,
                double, double, double,
                double, double, double>())
        .def("add_motor", &motor_ensemble::add_motor)
        .def_property_readonly("num_motors", &motor_ensemble::get_nmotors)
        .def("motor", &motor_ensemble::get_motor)
        // settings
        .def("set_binding_two", &motor_ensemble::set_binding_two)
        .def("set_bending", &motor_ensemble::set_bending)
        .def("use_shear", &motor_ensemble::use_shear)
        .def("use_static", &motor_ensemble::use_static)
        .def("set_parallel_alignment", &motor_ensemble::set_par)
        .def("set_antiparallel_alignment", &motor_ensemble::set_antipar)
        .def("set_alignment", &motor_ensemble::set_align)
        .def("kill_heads", &motor_ensemble::kill_heads)
        .def("unbind_all_heads", &motor_ensemble::unbind_all_heads)
        .def("revive_heads", &motor_ensemble::revive_heads)
        .def("set_external", &motor_ensemble::set_external)
        // dynamics
        .def("montecarlo", &motor_ensemble::montecarlo)
        .def("integrate", &motor_ensemble::integrate)
        .def("add_delrx", &motor_ensemble::update_d_strain)
        .def("compute_forces", &motor_ensemble::compute_forces)
        .def("update_energies", &motor_ensemble::update_energies)
        // get thermo
        .def_property_readonly("pe", &motor_ensemble::get_potential_energy)
        .def_property_readonly("pe_stretch", &motor_ensemble::get_stretching_energy)
        .def_property_readonly("pe_bend", &motor_ensemble::get_bending_energy)
        .def_property_readonly("pe_align", &motor_ensemble::get_alignment_energy)
        .def_property_readonly("pe_external", &motor_ensemble::get_external_energy)
        .def_property_readonly("vir", &motor_ensemble::get_potential_virial)
        .def_property_readonly("vir_stretch", &motor_ensemble::get_stretching_virial)
        .def_property_readonly("vir_bend", &motor_ensemble::get_bending_virial)
        .def_property_readonly("vir_align", &motor_ensemble::get_alignment_virial)
        .def_property_readonly("vir_external", &motor_ensemble::get_external_virial)
        // output state
        .def("output", &motor_ensemble::output)
        ;
}

PYBIND11_MODULE(afines, m)
{
    vec_bindings(m);
    globals_bindings(m);
    box_bindings(m);
    quadrants_bindings(m);
    generate_bindings(m);
    potentials_bindings(m);
    ext_bindings(m);
    filament_bindings(m);
    filament_ensemble_bindings(m);
    motor_bindings(m);
    motor_ensemble_bindings(m);
}
