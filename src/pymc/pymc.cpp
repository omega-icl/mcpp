#include <pybind11/pybind11.h>

namespace py = pybind11;

void mc_mcfunc    ( py::module_ & );
void mc_interval  ( py::module_ & );
void mc_mccormick ( py::module_ & );
void mc_supmodel  ( py::module_ & );
void mc_ffunc     ( py::module_ & );
void mc_fflin     ( py::module_ & );
void mc_ffvect    ( py::module_ & );
void mc_ffcustom  ( py::module_ & );

PYBIND11_MODULE( pymc, m )
{

  m.doc() = "Python interface of library MC++";

  mc_mcfunc    ( m );
  mc_interval  ( m );
  mc_mccormick ( m );
  mc_supmodel  ( m );
  mc_ffunc     ( m );
  mc_fflin     ( m );
  mc_ffvect    ( m );
  mc_ffcustom  ( m );
  
  m.attr("__version__") = "4";
}

