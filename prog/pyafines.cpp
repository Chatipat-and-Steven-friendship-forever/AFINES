#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "box.h"
#include "bead.h"
#include "spring.h"
#include "filament.h"
#include "filament_ensemble.h"
#include "motor.h"
#include "motor_ensemble.h"
#include "generate.h"
#include "globals.h"

namespace py = pybind11;

PYBIND11_MODULE(pyafines, m) {

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

    py::class_<bead>(m, "Bead")
        // constructors
        .def(py::init<double, double, double, double>())
        .def(py::init<bead const &>())
        // get parameters
        .def_property_readonly("radius", &bead::get_length)
        .def_property_readonly("friction", &bead::get_friction)
        .def_property_readonly("viscosity", &bead::get_viscosity)
        // get/set state
        .def_property("x", &bead::get_xcm, &bead::set_xcm)
        .def_property("y", &bead::get_ycm, &bead::set_ycm)
        .def_property("fx", &bead::get_fx, &bead::set_fx)
        .def_property("fy", &bead::get_fy, &bead::set_fy)
        .def("zero_forces", &bead::reset_force)
        // output state
        .def("write", &bead::write)
        ;

    py::class_<spring>(m, "Spring")
        // constructors
        .def(py::init<double, double, double, filament *, array<int, 2>>())
        // get/set state
        .def_property_readonly("x", &spring::get_hx)
        .def_property_readonly("y", &spring::get_hy)
        .def_property_readonly("length", &spring::get_length)
        .def_property_readonly("force", &spring::get_force)
        .def_property_readonly("direction", &spring::get_direction)
        .def_property_readonly("displacement", &spring::get_disp)
        // get parameters
        .def_property_readonly("kl", &spring::get_kl)
        .def_property_readonly("l0", &spring::get_l0)
        // get thermo
        .def_property_readonly("pe_stretch", &spring::get_stretching_energy)
        .def_property_readonly("vir_stretch", &spring::get_virial)
        // output state
        .def("write", &spring::write)
        // methods
        .def("get_positions_from_filament", &spring::step)
        .def("update_forces", &spring::update_force)
        .def("apply_forces_to_filament", &spring::filament_update)
        // intpoint
        .def_property_readonly("intpoint", &spring::get_intpoint)
        .def("update_intpoint", &spring::calc_intpoint)
        ;

    py::class_<filament>(m, "Filament")
        .def(py::init<filament_ensemble *, vector<bead *>, double, double, double, double, double, double, double, double>())
        // get/set state
        .def_property_readonly("num_beads", &filament::get_nbeads)
        .def_property_readonly("num_springs", &filament::get_nsprings)
        .def("bead", &filament::get_bead)
        .def("spring", &filament::get_spring)
        .def("add_bead", &filament::add_bead)
        // get thermo
        .def_property_readonly("ke", &filament::get_kinetic_energy)
        .def_property_readonly("pe", &filament::get_potential_energy)
        .def_property_readonly("te", &filament::get_total_energy)
        .def_property_readonly("pe_stretch", &filament::get_stretching_energy)
        .def_property_readonly("pe_bend", &filament::get_bending_energy)
        .def_property_readonly("vir_stretch", &filament::get_stretching_virial)
        .def_property_readonly("vir_bend", &filament::get_bending_virial)
        // output state
        .def("write_beads", &filament::write_beads)
        .def("write_springs", &filament::write_springs)
        // output thermo
        .def("write_thermo", &filament::write_thermo)
        // methods
        .def("add_delrx", &filament::update_d_strain)
        .def("update_pos", &filament::update_positions)
        .def("add_stretch_forces", &filament::update_stretching)
        .def("add_bend_forces", &filament::update_bending)
        .def("add_forces", &filament::update_forces)
        .def("fracture", &filament::fracture)
        // pulling
        .def_property_readonly("end_to_end", &filament::get_end2end)
        .def("pull_on_ends", &filament::pull_on_ends)
        .def("affine_pull", &filament::affine_pull)
        ;

    py::class_<filament_ensemble>(m, "FilamentEnsemble")
        .def(py::init<box *, vector<vector<double>>, array<int, 2>, double, double, double, double, double, double, double, double, bool>())
        // quadrants
        .def("disable_quads", &filament_ensemble::turn_quads_off)
        .def("update_quads", &filament_ensemble::quad_update_serial)
        .def("attach_map", &filament_ensemble::get_dist)
        .def("attach_map_all", &filament_ensemble::get_dist_all)
        .def("attach_list", &filament_ensemble::get_attach_list)
        .def_property_readonly("nq", &filament_ensemble::get_nq)
        // get/set state
        .def_property_readonly("friction", &filament_ensemble::get_bead_friction)
        .def_property_readonly("box", &filament_ensemble::get_box)
        .def_property_readonly("num_filaments", &filament_ensemble::get_nfilaments)
        .def_property_readonly("num_beads", &filament_ensemble::get_nbeads)
        .def_property_readonly("num_springs", &filament_ensemble::get_nsprings)
        .def("filament", &filament_ensemble::get_filament)
        // get thermo
        .def_property_readonly("pe_stretch", &filament_ensemble::get_stretching_energy)
        .def_property_readonly("pe_bend", &filament_ensemble::get_bending_energy)
        .def_property_readonly("vir_stretch", &filament_ensemble::get_stretching_virial)
        .def_property_readonly("vir_bend", &filament_ensemble::get_bending_virial)
        // output state
        .def("write_beads", &filament_ensemble::write_beads)
        .def("write_springs", &filament_ensemble::write_springs)
        // output thermo
        .def("write_thermo", &filament_ensemble::write_thermo)
        // methods
        .def("add_delrx", &filament_ensemble::update_delrx)
        .def("update_pos", &filament_ensemble::update_positions)
        .def("add_stretch_forces", &filament_ensemble::update_stretching)
        .def("add_bend_forces", &filament_ensemble::update_bending)
        .def("add_forces", &filament_ensemble::update_forces)
        .def("update_energies", &filament_ensemble::update_energies)
        .def("step", &filament_ensemble::update)
        // misc
        .def("add_circle_wall", &filament_ensemble::set_circle_wall)
        .def("spring_spring_intersections", &filament_ensemble::spring_spring_intersections)
        .def_property_readonly("broken", &filament_ensemble::get_broken)
        .def("clear_broken", &filament_ensemble::clear_broken)
        ;

    py::class_<motor>(m, "Motor")
        // constructors
        .def(py::init<array<double, 4>, double, filament_ensemble *, array<int, 2>, array<int, 2>, array<int, 2>, double, double, double, double, double, double, double, double, double, double, double>())
        // get/set state
        .def_property_readonly("f_index", &motor::get_f_index)
        .def_property_readonly("l_index", &motor::get_l_index)
        .def_property_readonly("pos_a_end", &motor::get_pos_a_end)
        .def_property_readonly("force", &motor::get_force)
        .def_property_readonly("state", &motor::get_states)
        .def_property_readonly("hx", &motor::get_hx)
        .def_property_readonly("hy", &motor::get_hy)
        // get thermo
        .def_property_readonly("ke", &motor::get_kinetic_energy)
        .def_property_readonly("pe_stretch", &motor::get_stretching_energy)
        .def_property_readonly("vir_stretch", &motor::get_virial)
        // output state
        .def("write", &motor::write)
        // change state
        .def("relax_head", &motor::relax_head)
        .def("kill_head", &motor::kill_head)
        .def("detach_head", (void (motor::*)(int)) &motor::detach_head)
        .def("detach_head", (void (motor::*)(int, array<double, 2>)) &motor::detach_head)
        .def("detach_head_without_moving", &motor::detach_head_without_moving)
        .def("revive_head", &motor::revive_head)
        .def("deactivate_head", &motor::deactivate_head)
        // methods
        .def("add_delrx", &motor::update_d_strain)
        .def("get_positions_from_filament", &motor::update_position_attached)
        .def("brownian_relax", &motor::brownian_relax)
        .def("update_angle", &motor::update_angle)
        .def("update_force", &motor::update_force)
        .def("attach", &motor::attach)
        .def("attach_opt", &motor::attach_opt)
        .def("step_head", &motor::step_onehead)
        .def("apply_forces_to_filament", &motor::filament_update_hd)
        // misc
        .def("allowed_bind", &motor::allowed_bind)
        .def("generate_off_pos", &motor::generate_off_pos)
        .def("metropolis_prob", &motor::metropolis_prob)
        .def("update_pos_a_end", &motor::update_pos_a_end)
        ;

    py::class_<motor_ensemble>(m, "MotorEnsemble")
        // constructors
        .def(py::init<vector<vector<double>>, double, double, double, filament_ensemble *, double, double, double, double, double, double, double, double, double, bool>())
        .def_property_readonly("num_motors", &motor_ensemble::get_nmotors)
        // get thermo
        .def("pe_stretch", &motor_ensemble::get_potential_energy)
        .def("vir_stretch", &motor_ensemble::get_virial)
        // output state
        .def("write", &motor_ensemble::motor_write)
        // change state
        .def("kill_heads", &motor_ensemble::kill_heads)
        .def("unbind_all_heads", &motor_ensemble::unbind_all_heads)
        .def("revive_heads", &motor_ensemble::revive_heads)
        // methods
        .def("add_motor", &motor_ensemble::add_motor)
        .def("add_delrx", &motor_ensemble::update_d_strain)
        .def("walk", &motor_ensemble::motor_walk)
        .def("step", &motor_ensemble::motor_update)
        .def("detach_from_broken", &motor_ensemble::check_broken_filaments)
        .def("update_energies", &motor_ensemble::update_energies)
        ;

    m.def("generate_filaments", (vector<vector<double>> (*)(box *, int, int, int, int, double, double, double, double, vector<array<double, 3>>, double, double)) &generate_filament_ensemble);
    m.def("generate_filaments", (vector<vector<double>> (*)(box *, double, double, double, double, int, double, vector<array<double, 3>>, double, double)) &generate_filament_ensemble);
    m.def("generate_motors", generate_motor_ensemble);

}
