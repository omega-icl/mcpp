#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ffunc.hpp"
#include "smon.hpp"
#include "spoly.hpp"

namespace py = pybind11;

void
mc_ffmon(py::module_& m)
{
  typedef mc::FFVar const* KEY;
  typedef mc::lt_FFVar COMP;
  typedef mc::SMon<KEY, COMP> FFMon;
  typedef mc::SPoly<KEY, COMP> FFPoly;

  py::class_<FFMon> pyFFMon(m, "FFMon");

  // --- Exceptions Nested Class ---
  py::class_<FFMon::Exceptions> pyFFMonExceptions(pyFFMon, "Exceptions");
  py::enum_<FFMon::Exceptions::TYPE>(pyFFMonExceptions, "TYPE")
      .value("SUB", FFMon::Exceptions::SUB,
             "Subtraction of a monomial that is not a proper subset")
      .value("DIV", FFMon::Exceptions::DIV,
             "Division by a factor greater than the greatest common exponent")
      .value("CONV", FFMon::Exceptions::CONV,
             "Conversion with different variable indexing failed")
      .export_values();

  pyFFMonExceptions.def("ierr", &FFMon::Exceptions::ierr, "Return error flag")
      .def("what", &FFMon::Exceptions::what, "Return error description");

  // --- FFMon Main Class ---
  pyFFMon
      // Constructors
      .def(py::init<>(), "Default constructor (constant 1)")
      .def(py::init<KEY const&, unsigned const>(), py::arg("var"),
           py::arg("ord") = 1, "Constructor from variable and exponent")
      .def(py::init<std::map<KEY, unsigned, COMP> const&>(), py::arg("expr"),
           "Constructor from dictionary of variables and exponents")
      .def(py::init<FFMon const&>(), "Copy constructor")
      // Accessors
      .def_readonly("tord", &FFMon::tord, "Monomial total order")
      .def_property_readonly(
          "expr",
          //&FFMon::expr,
          [](FFMon& self)
          {
            auto const& expr = self.expr;
            return std::vector<std::pair<KEY, unsigned> >(expr.cbegin(),
                                                          expr.cend());
          },
          "List of variables and exponents tuples")
      // Methods
      .def("exp", &FFMon::exp, py::arg("var"),
           "Exponent of variable 'var' in monomial")
      .def("gcexp", &FFMon::gcexp, "Greatest common exponent of terms")
      .def("lexp", &FFMon::lexp, "Least exponent among all terms")
      .def("gexp", &FFMon::gexp, "Greatest exponent among all terms")
      .def("hull", &FFMon::hull, py::arg("mon"), "Union with other monomial")
      .def("inter", &FFMon::inter, py::arg("mon"), "Test for intersection")
      .def("subset", &FFMon::subset, py::arg("mon"), "Test for proper subset")
      .def("subseteq", &FFMon::subseteq, py::arg("mon"), "Test for subset")
      .def(
          "display", &FFMon::display, py::arg("basis") = 0,
          "String representation with basis option - 0: Monomial; 1: Chebyshev")
      // Operators
      .def(
          "__hash__",
          [](FFMon const& self)
          {
            // Hash is computed based on 'tord' and the 'expr' map content.
            // We convert the map to a tuple of items (which is sorted by key in
            // C++ map) to ensure consistent hashing in Python.
            py::list expr_items;
            for (auto const& [var, ord] : self.expr)
            {
              expr_items.append(py::make_tuple(var, ord));
            }
            // Return hash((tord, tuple(items)))
            return py::hash(py::make_tuple(self.tord, py::tuple(expr_items)));
          },
          "Compute hash of the monomial")
      .def(
          "__getitem__", [](FFMon const& self, KEY const& var)
          { return self[var]; }, "Extract sub-monomial for variable")
      .def(py::self += py::self)
      .def(py::self -= py::self)
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self *= unsigned())
      .def(py::self /= unsigned())
      .def(py::self * unsigned())
      .def(py::self / unsigned())

      // Comparison operators
      .def(
          "__eq__", [](FFMon const& Mon1, FFMon const& Mon2)
          { return Mon1 == Mon2; }, py::is_operator())
      .def(
          "__ne__", [](FFMon const& Mon1, FFMon const& Mon2)
          { return !(Mon1 == Mon2); }, py::is_operator())

      // String representation
      .def("__str__", [](FFMon const& self)
           { return self.display(FFPoly::options.BASIS); })
      .def("__repr__", [](FFMon const& self)
           { return self.display(FFPoly::options.BASIS); });
}
