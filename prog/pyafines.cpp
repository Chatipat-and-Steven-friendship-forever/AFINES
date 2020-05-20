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

    py::class_<filament_ensemble>(m, "FilamentEnsemble");
    py::class_<motor>(m, "Motor");
    py::class_<motor_ensemble>(m, "MotorEnsemble");

    m.def("generate_filaments", (vector<vector<double>>(*)(box *, int, int, int, int, double, double, double, double, vector<array<double, 3>>, double, double)) generate_filament_ensemble);
    m.def("generate_filaments", (vector<vector<double>>(*)(box *, double, double, double, double, int, double, vector<array<double, 3>>, double, double)) generate_filament_ensemble);
    m.def("generate_motors", generate_motor_ensemble);

}
