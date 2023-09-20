#include <pybind11/pybind11.h>
#include "ibd_call.cpp"

namespace py = pybind11;

PYBIND11_MODULE(ibd_call, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("call_ibds", &call_ibds, "A function that calls IBDs(matches) in sites of specified length");
}
