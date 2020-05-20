#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

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

    py::class_<spring>(m, "Spring");
    py::class_<filament>(m, "Filament");
    py::class_<filament_ensemble>(m, "FilamentEnsemble");
    py::class_<motor>(m, "Motor");
    py::class_<motor_ensemble>(m, "MotorEnsemble");

    m.def("generate_filaments", (vector<vector<double>>(*)(box *, int, int, int, int, double, double, double, double, vector<array<double, 3>>, double, double)) generate_filament_ensemble);
    m.def("generate_filaments", (vector<vector<double>>(*)(box *, double, double, double, double, int, double, vector<array<double, 3>>, double, double)) generate_filament_ensemble);
    m.def("generate_motors", generate_motor_ensemble);

}
