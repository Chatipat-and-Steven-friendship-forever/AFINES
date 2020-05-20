#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "bead.h"
#include "spring.h"
#include "filament.h"
#include "filament_ensemble.h"
#include "motor.h"
#include "motor_ensemble.h"
#include "globals.h"

namespace py = pybind11;

PYBIND11_MODULE(pyafines, m) {
    py::class_<bead>(m, "Bead");
    py::class_<spring>(m, "Spring");
    py::class_<filament>(m, "Filament");
    py::class_<filament_ensemble>(m, "FilamentEnsemble");
    py::class_<motor>(m, "Motor");
    py::class_<motor_ensemble>(m, "MotorEnsemble");
}
