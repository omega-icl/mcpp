// Copyright (C) 2012, 2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODESTRUCT_HPP
#define MC__ODESTRUCT_HPP

namespace mc
{

//! @brief C++ class defining the structure of parametric ODE problems
class ODESTRUCT
{
protected:
  //! @brief Number of independent variables (parameters)
  unsigned int _np;
  //! @brief Number of dependent variables (differential states)
  unsigned int _nx;
  //! @brief Number of time-invariant functions (parameters)
  unsigned int _ni;
  
public:
  //! @brief Constructor
  ODESTRUCT
    ( const unsigned int np, const unsigned int nx, const unsigned int ni=0 )
    : _np(np), _nx(nx), _ni(ni)
    {}

  //! @brief Returns number of independent variables (parameters)
  unsigned int np() const
    { return _np; }
  //! @brief Returns number of dependent variables (differential states)
  unsigned int nx() const
    { return _nx; }
  //! @brief Returns number of time-invariant functions (parameters)
  unsigned int ni() const
    { return _ni; }

  template <typename TX, typename TP, typename TT>
  TX INV
    ( const unsigned int ii, const TP*p, const TX*x, const TT&t,
      const unsigned int is )
    {
      throw std::runtime_error("invariant undefined");
    }

  template <typename TX, typename TP, typename TT>
  TX BND
    ( const unsigned int ix, const TP*p, const TX*x, const TT&t,
      const unsigned int is )
    {
      return x[ix];
      //throw std::runtime_error("natural bounds undefined");
    }
  
private:
  //! @brief Default constructor (private methods to block default compiler methods)
  ODESTRUCT();
};

} // end namespace mc

#endif
