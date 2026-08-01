#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace pybind11
{
namespace detail
{

// Partial specialization for std::multiset
template <typename Type, typename Alloc>
struct type_caster<std::multiset<Type, std::less<Type>, Alloc>>
{
  using value_type = std::multiset<Type, std::less<Type>, Alloc>;

  // Python -> C++ (Convert list to multiset)
  bool
  load(handle src, bool convert)
  {
    if (!isinstance<pybind11::list>(src)) return false;
    auto l = reinterpret_borrow<pybind11::list>(src);
    for (auto item : l)
    {
      value.insert(item.cast<Type>());
    }
    return true;
  }

  // C++ -> Python (Convert multiset to list)
  static handle
  cast(const value_type& src, return_value_policy policy, handle parent)
  {
    pybind11::list l;
    for (const auto& item : src)
    {
      // convert each C++ element to a Python object and append
      l.append(pybind11::cast(item, policy, parent));
    }
    return l.release();
  }

  PYBIND11_TYPE_CASTER(value_type,
                       _("List[") + make_caster<Type>::name + _("]"));
};

}  // namespace detail
}  // namespace pybind11

namespace py = pybind11;

void mc_mcfunc(py::module_&);
void mc_interval(py::module_&);
void mc_mccormick(py::module_&);
void mc_specbnd(py::module_&);
void mc_tmodel(py::module_&);
void mc_cmodel(py::module_&);
void mc_scmodel(py::module_&);
void mc_sicmodel(py::module_&);
void mc_supmodel(py::module_&);
void mc_polimage(py::module_&);
void mc_ellimage(py::module_&);
void mc_ffdep(py::module_&);
void mc_ffinv(py::module_&);
void mc_ffunc(py::module_&);
void mc_ffmon(py::module_&);
void mc_ffpoly(py::module_&);
void mc_smon(py::module_&);
void mc_spoly(py::module_&);
void mc_slift(py::module_&);
void mc_fflin(py::module_&);
void mc_ffmlp(py::module_&);
void mc_ffdagext(py::module_&);
void mc_ffcustom(py::module_&);
void mc_ffvect(py::module_&);

PYBIND11_MODULE(pymcpp, m)
{
  m.doc() = R"doc(
Python bindings of MC++: guaranteed bounds and convex/concave
relaxations of factorable functions, for global and robust
optimization.

A factorable function can be represented as a directed acyclic graph
(``FFGraph``/``FFVar``) and evaluated in different arithmetics, or the
arithmetics can be used directly via operator overloading. The main
types are:

- ``Interval``: interval bounds on the range of a function;
- ``McCormick``: McCormick convex/concave relaxations with subgradient
  propagation;
- ``TModel``/``CModel``/``SCModel``/``SICModel``: Taylor and Chebyshev
  models; ``Specbnd``: spectral bounds; ``SupModel``: superposition
  relaxations; ``PolImg``: polyhedral relaxations; ``EllImg``:
  ellipsoidal calculus;
- ``FFGraph``/``FFVar``: expression DAGs supporting evaluation in all
  of the above, symbolic differentiation and Taylor expansion.

Example: bound f(x,y) = x*(exp(x)-y)**2 on [-2,1] x [-1,2]:

>>> import pymcpp
>>> x, y = pymcpp.Interval(-2, 1), pymcpp.Interval(-1, 2)
>>> f = x * (pymcpp.exp(x) - y)**2
>>> print(f)
[ -2.76512e+01 :  1.38256e+01 ]

See the Jupyter notebooks in the ``notebook/`` directory of the MC++
distribution (interval.ipynb, mccormick.ipynb, ffunc.ipynb, ...) for
worked examples of each arithmetic.
)doc";

  mc_mcfunc(m);
  mc_interval(m);
  mc_mccormick(m);
  mc_specbnd(m);
  mc_tmodel(m);
  mc_cmodel(m);
  mc_scmodel(m);
  mc_sicmodel(m);
  mc_supmodel(m);
  mc_polimage(m);
  mc_ellimage(m);
  mc_smon(m);
  mc_spoly(m);
  mc_ffdep(m);
  mc_ffinv(m);
  mc_ffunc(m);
  mc_ffmon(m);
  mc_ffpoly(m);
  mc_slift(m);
  mc_fflin(m);
  mc_ffmlp(m);
  mc_ffdagext(m);
  mc_ffcustom(m);
  mc_ffvect(m);

  m.attr("__version__") = "5.0.2";
}
