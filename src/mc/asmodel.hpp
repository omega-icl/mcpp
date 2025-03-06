// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_ASM Abstract Superposition Model Arithmetic for Factorable Functions
\author Yanlin Zha, Beno&icirc;t Chachuat

Superposition decomposition is a technique that underlies set-arithmetic for computing non-convex enclosures in the form of superposition of univariate estimators. 
Formally, to be modified 

The classes mc::ASModel and mc::ASVar provide an implementation of ASM arithmetic based on the operator/function overloading mechanism of C++. This makes ASM both simple and intuitive to compute, similar to computing function values in real arithmetics or function bounds in interval arithmetic (see \ref page_INTERVAL). mc::ASModel and mc::ASVar are templated in the type used to propagate the coefficients of various forms. 
In addition, mc::ASVar can be used as the template parameter of other available types in MC++; as well as types in <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A> for computing ASM of the partial derivatives or Taylor coefficients of a factorable function.


\section sec_ASM_use How do I compute an ASM of a factorable function? (to be modified)



Further exceptions may be thrown by the template class itself.


\section sec_ASM_refs References

- J. Su, Y. Zha, K. Wang, M.E. Villanueva, R. Paulen, B. Houska. <A HREF="https://doi.org/10.1016/j.ifacol.2019.06.124">Interval superposition arithmetic for guaranteed parameter estimation</A>, <I>IFAC-PapersOnLine</I> <B>52</B>(1):574-579, 2019.
- Y. Zha, M.E. Villanueva, B. Houska. <A HREF="https://arxiv.org/abs/1610.05862">Interval superposition arithmetic</A>, <I>ArXiv</I>, 1610.05862v2, 2018.
.
*/



#ifndef MC__ASMODEL_HPP
#define MC__ASMODEL_HPP

#define MC__ASMODEL_NOT_DEBUG_SHADOW
//#define TEST_MOVE

//#define FILIB__COMPUTATION_DEBUG
#include <iostream>
#include <iomanip> 
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <bitset>
#include <cassert>
#include "mcop.hpp"
#include "mcfunc.hpp"
#include "univarmodels.hpp" 

//#define MC__ASMODEL_TRACE
//#undef  MC__ASMODEL_DEBUG_PROD

namespace mc
{

const double MC__ASM_COMPUTATION_TOL(1e-15);
template <typename T> class ASVar;

//! @brief C++ class for interval superposition models of factorable functions: environment
////////////////////////////////////////////////////////////////////////
//! mc::ASModel is a C++ class for definition of interval superposition
//! model (ASM) environment. ASM propagation of factorable functions is
//! implemented via the C++ class mc::ASVar. The template parameter
//! corresponds to the type used to propagate the interval coefficients.
////////////////////////////////////////////////////////////////////////
template <typename T> 
class ASModel
////////////////////////////////////////////////////////////////////////
{
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const ASModel<U>& );

  template <typename U> friend class ASVar;

  template <typename U> friend ASVar<U> max
    ( ASVar<U> const&, ASVar<U> const& );
  template <typename U> friend ASVar<U> max
    ( ASVar<U> const&, double const& );
  template <typename U> friend ASVar<U> max
    ( ASVar<U> &&, double const& );
  // template <typename U> friend ASVar<U> min
  //   ( ASVar<U> const&, ASVar<U> const& );
  // template <typename U> friend ASVar<U> min
  //   ( ASVar<U> const&, double const& );
  // template <typename U> friend ASVar<U> min
  //   ( ASVar<U> &&, double const& );

  // template <typename U> friend ASVar<U> inv
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> inv
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> sqr
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> sqr
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> sqrt
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> sqrt
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> fabs
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> fabs
  //   ( ASVar<U> && );  
  template <typename U> friend ASVar<U> relu
    ( ASVar<U> const& );
  template <typename U> friend ASVar<U> relu
    ( ASVar<U> && );    
  // template <typename U> friend ASVar<U> exp
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> exp
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> log
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> log
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> xlog
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> xlog
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> cos
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> cos
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> sin
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> sin
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> pow
  //   ( ASVar<U> const&, int const& n );
  // template <typename U> friend ASVar<U> pow
  //   ( ASVar<U> &&, int const& n );
  // template <typename U> friend ASVar<U> pow
  //   ( ASVar<U> const&, double const& );
  // template <typename U> friend ASVar<U> pow
  //   ( ASVar<U> &&, double const&);
  // template <typename U> friend ASVar<U> tanh
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> tanh
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> affine_transform
  //   ( std::vector<ASVar<U> >const& , std::vector<double> const&, const double);
  // template <typename U> friend ASVar<U> affine_transform
  //   ( std::vector<ASVar<U> >&& , std::vector<double> const&, const double);    
  // template <typename U> friend ASVar<U> affine_transform
  //   ( std::vector<std::reference_wrapper<ASVar<U>> > , std::vector<double> const&, const double);    
  // template <typename U> friend ASVar<U> intersect
  //   ( ASVar<U> const& , U);
  // template <typename U> friend ASVar<U> intersect
  //   ( ASVar<U> && , U);

//  template <typename U> friend ASVar<U> mid 
//  (const ASVar<U>&, const ASVar<U>&, const double);
 private:

  //! @brief Number of partitions (to support up to 2^64 partitions)
  unsigned long long _ndiv;
  //! @brief Number of variables
  unsigned int _nvar;
  //! @brief Whether variables are defined or not
  std::vector<bool> _defvar;
  //! @brief Variable bounds
  std::vector<T> _bndvar;
  // Internal containing the parition sizes of each variable
  std::vector<double> _psize;  

  // Internal Intermediate containing lower bounds of coefficient matrix rows
  mutable std::vector<double> _L1;
  // Internal Intermediate containing lower bounds of coefficient matrix rows
  mutable std::vector<double> _L2;
  // Internal Intermediate containing upper bounds of coefficient matrix rows
  mutable std::vector<double> _U1;
  // Internal Intermediate containing upper bounds of coefficient matrix rows
  mutable std::vector<double> _U2;
  // Internal Intermediate containing reference of coefficient matrix rows
  mutable std::vector<double> _c1;
  // Internal Intermediate containing references of coefficient matrix rows
  mutable std::vector<double> _c2;
  // Internal Intermediate containing radius of coefficient matrix rows
  mutable std::vector<double> _r1;
  // Internal Intermediate containing radius of coefficient matrix rows
  mutable std::vector<double> _r2;
  // Internal Flag indicating whether the internal containers are updated
  mutable bool _IntmdtCntnrSeted;
  // temporary matrix 
  mutable std::vector<std::vector<double>> _COve;
  // temporary matrix 
  mutable std::vector<std::vector<double>> _CUnd;
  // temporary matrix 
  mutable std::vector<std::vector<T>> _COut;
  
 public:
#ifdef MC__ASM_DEBUG_EVAL_LOGGING 
  static std::size_t debugging_ctr;
#endif
  //! @brief Constructor of ASM with <a>nvar</a> variables and <a>ndiv</a> partitions
  ASModel
  ( const unsigned int& nvar, const unsigned int& ndiv )
  : _ndiv(ndiv), _nvar(nvar),_IntmdtCntnrSeted(false),_COve(0),_CUnd(0),_COut(0)
  {
    _bndvar.resize( _nvar );
    _defvar.resize( _nvar, false );
    
    _L1.resize( _nvar );
    _L2.resize( _nvar );
    _U1.resize( _nvar );
    _U2.resize( _nvar );
    _c1.resize( _nvar );
    _c2.resize( _nvar );
    _r1.resize( _nvar );
    _r2.resize( _nvar );
    _r2.resize( _nvar );
    _psize.resize( _nvar , 0. );
  }

  ~ASModel() 
  {
#ifdef FASM_LIFITIME_DEBUG     
    std::cout<< "ASM deleted, nvar = " <<_nvar <<std::endl;
#endif
//     _bndvar.resize( 0 );
//     _defvar.resize( 0 );
    
//     _L1.resize(0);
//     _L2.resize(0);
//     _U1.resize(0);
//     _U2.resize(0);
//     _c1.resize(0);
//     _c2.resize(0);
//     _r1.resize(0);
//     _r2.resize(0);
//     _r2.resize(0);
//     _psize.resize(0);
  }

  //! @brief Retrieve number of variables in ASM
  unsigned nvar
  ()
  const
  { return _nvar; };

  //! @brief Retrieve number of partitions in ASM
  unsigned ndiv
  ()
  const
  { return _ndiv; };

  //! @brief Retrieve partition sizes in ASM
  std::vector<double> const& psize
  ()
  const
  { return _psize; };

  std::vector<T> const& bndvar
  ()
  const
  { return _bndvar; };


 //! @brief Options of mc::ASModel
  static struct Options
  {
    //! @brief Constructor
    Options():
      ASYREM_USE(true), DCDEC_USE(true), SCALING_TYPE(FULL), INTERSECTION_USE(true), ENVEL_USE(true), ROOT_MAXIT(100), ROOT_TOL(1e-10), SLOPE_USE(false),SHADOW_USE(false),NSUB(16)
      {}

    //! @brief FASModel re-scaling option in binary product
    enum SCALING{
      NONE=0,	  //!< without using scaling
      PARTIAL,	//!< only re-scaling the radius of two mutiplicants(mutipliers) 
      FULL,	    //!< re-scaling the range of two mutiplicants(mutipliers) to [-1,1] 
      ADAPT     //!< adapted re-scaling
    };
    //! @brief Whether to use asymmetric inclusions for convex/concave terms as available
    bool ASYREM_USE;
    //! @brief Whether to use DC decomposition in product rule and composition rule with non-monotonic univariates
    bool DCDEC_USE;
    //! @brief Whether to use re-scaling in binary product, and which type is used 
    SCALING SCALING_TYPE; 
    //! @brief Whether to use intersection
    bool INTERSECTION_USE; 
    //! @brief Whether to use convex/concave envelopes of nonconvex terms as available  
    bool ENVEL_USE;
    //! @brief Maximal number of iterations in root search - Default: 100
    unsigned ROOT_MAXIT;
    //! @brief Termination tolerance in root search - Default: 1e-10
    double ROOT_TOL;
    //! @brief Whether to use slope-based enhancement - Default: false
    bool SLOPE_USE;    
    //! @brief Whether to use shadow enhancement - Default: false
    bool SHADOW_USE;   
    //! @brief Number of Subdivisions to Support PWC
    unsigned int NSUB;          
  } options;

  //! @brief Exceptions of mc::ASModel
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SCModel exception handling
    enum TYPE{
      DIV=1,	 //!< Division by zero scalar
      INV,	     //!< Inverse operation with zero in range
      LOG,	     //!< Log operation with non-positive numbers in range
      SQRT,	     //!< Square-root operation with negative numbers in range
      TAN,	     //!< Tangent operation with (k+1/2)·PI in range
      ROOT,      //!< Error during root search for obtaining the convex/concave envelope of a univariate term
      INTERN=-1, //!< Internal error
      INDEX=-2,  //!< Variable index out of range
      MODEL=-3,	 //!< Operation between variables belonging to different models
      UNDEF=-33  //!< Feature not yet implemented
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case DIV:
        return "mc::ASModel\t Division by zero scalar";
      case INV:
        return "mc::ASModel\t Inverse operation with zero in range";
      case LOG:
        return "mc::ASModel\t Log operation with non-positive numbers in range";
      case SQRT:
        return "mc::ASModel\t Square-root operation with negative numbers in range";
      case TAN:
        return "mc::ASModel\t Tangent operation with (k+1/2)·PI in range";
      case INDEX:
        return "mc::ASModel\t Variable index out of range";
      case MODEL:
        return "mc::ASModel\t Operation between variables belonging to different models not permitted";
      case UNDEF:
        return "mc::ASModel\t Feature not yet implemented";
      case INTERN:
      default:
        return "mc::ASModel\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

 private:

  //! @brief Return the bound of the ASM variable <a>var</a>
  T _B
  ( std::vector<UnivarPWL<T>> const& lst, const unsigned int rec=0 )
  const
  {
    assert( !lst.empty() );
    std::pair<double,double> B( 0., 0. );
    for( unsigned int i=0; i<_nvar; i++ ){
      if( lst[i].empty() ) continue;
      if( rec == 3 ){                     //  used for performaing intersection only;
        std::cout << "  Mode 3 in _B() has not yet been implemented "<< std::endl;
        throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); 
        //_L2[i] = _B_innerL( lst[i] );
        //_U2[i] = _B_innerU( lst[i] );
      }
      else{                               //  basic function
        //std::cout << "FOR ENDED " << i << std::endl;        
        std::pair<double,double> Brow = _B( lst[i] );
        if( rec == 1 ){
          _L1[i] = Brow.first; // Op<T>::l( Brow );
          _U1[i] = Brow.second; // Op<T>::u( Brow );
        }
        else if( rec == 2 ){
          _L2[i] = Brow.first; // Op<T>::l( Brow );
          _U2[i] = Brow.second; // Op<T>::u( Brow );
        }
        B.first += Brow.first;
        B.second += Brow.second;        
        // B += Brow;
      }
    } 
    //std::cout << "FOR ENDED" << std::endl;
    if( rec == 1 ) _IntmdtCntnrSeted = true;
    if(B.first > B.second){
      if(std::fabs(B.first - B.second) <= 1e2*MC__ASM_COMPUTATION_TOL){
        return T(B.second,B.first);        
      }
      else{
        std::cout << "  ERROR! THE lb of the ASM is greater than ub by " <<std::setprecision(18) << B.first - B.second << std::endl;
        throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); 
      }
    }
    return T(B.first,B.second);

    
  }

  std::pair<double,double> _B
  ( UnivarPWL<T> const& row )
  const
  {
    assert( !row.empty() );
    std::pair<double,double> iBnd( row.getBnd());
    return iBnd;
  }


  void _compute_COve(std::vector<UnivarPWL<T>> const& lst)   
  const
  { 
//    std::cout << "ove" << std::endl; 
    _COve.resize(_nvar); 
    //std::cout << "reized" << std::endl; 
    for (unsigned int i = 0; i < _nvar; i++){
      if (lst[i].empty()) continue;
        //std::cout << "to resize 0" << std::endl; 
        // _COve[i].resize(options.NSUB,0.);
        //std::cout << "resized" << std::endl; 
        //std::cout << "to set" << std::endl; 
        _COve[i] = lst[i].oveEst.get_PWC(options.NSUB);

    } 

  }  

  void _compute_CUnd(std::vector<UnivarPWL<T>> const& lst)
  const  
  { 
    //std::cout << "und" << std::endl; 
    _CUnd.resize(_nvar); 
    //std::cout << "reized" << std::endl; 
    for (unsigned int i = 0; i < _nvar; i++){
      if (lst[i].empty()) continue;
        //std::cout << "to resize 0" << std::endl; 
        // _CUnd[i].resize(options.NSUB,0.);
        //std::cout << "resized" << std::endl; 
        //std::cout << "to set" << std::endl; 
        _CUnd[i] = lst[i].undEst.get_PWC(options.NSUB);
        //std::cout << "set" << std::endl;  
    } 
  }
  
  void _compute_C(std::vector<UnivarPWL<T>> const& lst,unsigned int const ndep)
  const
  { 
    _compute_CUnd(lst);
    //std::cout << "und" << std::endl; 
    _compute_COve(lst);
    //std::cout << "ove" << std::endl; 

    double sum_r1 = 0.;
    for (unsigned int i=0; i<_nvar; i++){
      if (lst[i].empty()) continue;
      _r1[i] = _COve[i][0] - _CUnd[i][0];
      for (unsigned int j = 1; j < options.NSUB; j++){
        _r1[i] = std::min(_r1[i],_COve[i][j] - _CUnd[i][j]);
      }
      sum_r1 += _r1[i];
      //std::cout << i << std::endl;
    }

    if(sum_r1 < 0){
      std::cout << "Error in Display! sum_r1 = " << sum_r1 << std::endl;
            std::cout << "The domain is " << std::endl;
      for(unsigned int i=0; i<_nvar; i++ ){
        if (lst[i].empty()) continue;
        std::cout << _bndvar[i] <<std::endl;
        std::cout << std::defaultfloat << std::setprecision(2)  << "i: " << i << std::endl;
        std::cout << std::scientific << std::setprecision(12) << mc::Op<T>::l(_bndvar[i]) << " : " << mc::Op<T>::u(_bndvar[i]) << std::endl;
      }

      for( unsigned int i=0; i<_nvar; i++ ){
        std::cout << std::right << std::setw(5) << "row i = " << i << std::endl;
        for (unsigned int j = 0; j < options.NSUB; j++){
          std::cout << _COve[i][j] <<"    " << _CUnd[i][j] << "    ";
        }
        std::cout << std::endl;;   
      }
      for( unsigned int i=0; i<_nvar; i++ ){
        if (lst[i].empty()) continue;
        std::cout << std::right << std::setw(5) << "und " << i << std::endl;
        std::cout << lst[i].undEst <<std::endl;
        std::cout << std::right << std::setw(5) << "ove " << i << std::endl;
        std::cout << lst[i].oveEst <<std::endl; 
      }

      throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );  

	  
    } 
    sum_r1 = sum_r1/ndep;
    
    
    _COut.resize(_nvar);
    for (unsigned int i = 0; i < _nvar; i++){      
      if (lst[i].empty())  continue;
      _COut[i].resize(options.NSUB);
      // for (unsigned int j = 0; j < options.NSUB; j++)
      //   _COut[i][j] = T(0.);
      for (unsigned int j = 0; j < options.NSUB; j++){
        _COut[i][j] = T(_CUnd[i][j], _COve[i][j] - _r1[i] + sum_r1);
        
      }
    } 
  }

  
  void _compute_C(std::vector<double> const& lnr, double const & cst)
  const
  { 
    _COut.resize(_nvar);
    unsigned ndep = 0;

	// std::cout << std::right << std::setw(5) << "This ASVar is affine : " << std::endl;
    // std::cout << std::right << std::setw(5) << "cst" << std::setw(8);
    // for( unsigned int i=0; i<(lnr.size()); i++ ){
    //   if( lnr[i] == 0. ) continue;
    //   std::cout << "   " << _bndvar[i] << "  " << std::setw(8);   
    // }
    // std::cout << std::endl;
    // std::cout << std::right << std::setw(5) << cst << std::setw(8);
    // for( unsigned int i=0; i<(lnr.size()); i++ ){
    //   if( lnr[i] == 0. ) continue;
    //   std::cout << " + " << lnr[i] << " * x_" << i  << std::setw(8);   
    // }
    // std::cout << std::endl;

  
    for (unsigned int i = 0; i < _nvar; i++){
      if (lnr[i] == 0.){
        //_COut[i].resize(options.NSUB);
        //for (unsigned int j = 0; j < options.NSUB; j++)
        //  _COut[i][j] = T(0.);
        continue;
      } 
      ndep ++;
    }
    double tmp_cst = cst/ndep;
    for (unsigned int i = 0; i < _nvar; i++){
      if (lnr[i] == 0.) continue;
      else{
        _COut[i].resize(options.NSUB);
        double lb = mc::Op<T>::l(_bndvar[i]);
        double ub = mc::Op<T>::u(_bndvar[i]);
        double h = (ub - lb)/(double) options.NSUB;
        for (unsigned int j = 0; j < options.NSUB - 1; j++){
          ub = lb+h;
          _COut[i][j] = T(lb, ub)*lnr[i] + tmp_cst;
          lb = ub;
        }
        ub = mc::Op<T>::u(_bndvar[i]);
        _COut[i][options.NSUB-1] = T(lb, ub)*lnr[i] + tmp_cst;
      }
    } 

    // for( unsigned int i=0; i<(lnr.size()); i++ ){
    //   std::cout << std::right << std::setw(5) << "row i = " << i << std::setw(8);
    //   for (unsigned int j = 0; j < options.NSUB; j++){
    //     std::cout << _COut[i][j] << std::setw(8);
    //   }
    //   std::cout << std::endl;;   
    // }
    
    
  }

  void _compute_C(double const & cst)
  const
  { 
    _COut.resize(_nvar);
    _COut[0].resize(options.NSUB);
    for (unsigned int j = 0; j < options.NSUB; j++){
      _COut[0][j] = T(cst);
    }      
    // double tmp_cst = cst/(double) _nvar;
    // for (unsigned int i = 0; i < _nvar; i++){
    //   _COut[i].resize(options.NSUB);
    //   for (unsigned int j = 0; j < options.NSUB; j++)
    //     _COut[i][j] = T(tmp_cst);
    // } 
  }

  void _shadowInfo_init(std::vector<double> & shadow_info)
  const
  {
    shadow_info.resize(3);
    shadow_info[0]=0;  // number of shadow underestimator(s)
    shadow_info[1]=0;  // number of shadow overestimator(s) 
  }



  // template <typename PUNIV>
  // void _asym
  // ( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f,
  //   double const& zopt, bool const cvx )
  // const;
  // template <typename PUNIV>
  // void _asym
  // ( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f,
  //   double const& zopt, bool const cvx, T const& bnd )
  // const;

  // template <typename PUNIV, typename BNUIA>
  // void _asym_slope  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<double>>>& slope, const std::vector<double>& partitionSize, unsigned const& ndep, PUNIV const& f, BNUIA const& fDerv, double const& zopt, bool const cvx )
  // const;

  // template <typename PUNIV , typename BNUIA>
  // void _asym_slope  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<double>>>& slope, const std::vector<double>& partitionSize, unsigned const& ndep, PUNIV const& f, BNUIA const& fDerv, double const& zopt, bool const cvx, T const& bnd )
  // const;


  void _asym_relu (std::vector<UnivarPWL<T>>& lst, unsigned const& ndep )
  const;


  void _asym_relu_withShadow_old ( std::vector<UnivarPWL<T>>& lst,unsigned const& ndep,std::vector<UnivarPWL<T>>& shadow)
  const;

  void _asym_relu_withShadow ( std::vector<UnivarPWL<T>>& lst,unsigned const& ndep,std::vector<UnivarPWL<T>>& shadow, std::vector<double>& shadow_info)
  const;

  void _add_aggregate_shadow (std::vector<UnivarPWL<T>>& Alst, const std::vector<UnivarPWL<T>>& Blst,std::vector<UnivarPWL<T>>& Ashadow, const std::vector<UnivarPWL<T>>& Bshadow, 
  unsigned & Andep,std::vector<double>& Ashadow_info, const std::vector<double>& Bshadow_info)
  const;

  // void _asym_relu_withShadow_major ( std::vector<UnivarPWL<T>>& lst,unsigned const& ndep,std::vector<UnivarPWL<T>>& shadow, double lambda, double mu)
  // const;

  // void _asym_relu_withShadow_other ( std::vector<UnivarPWL<T>>& lst,unsigned const& ndep,std::vector<UnivarPWL<T>>& shadow, char mode)
  // const;


  // void _asym_slope_relu_shadow ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<double>>>& slope, const std::vector<double>& partitionSize, unsigned const& ndep, std::vector<std::vector<std::vector<double>>>& shadow, std::vector<std::vector<std::vector<double>>>& shadow_slope, std::vector<double> & shadow_info )
  // const;


  // void _add_aggregate_shadow ( std::vector<std::vector<T>>& Amat, const std::vector<std::vector<T>>& Bmat, std::vector<std::vector<std::vector<double>>>& Aslope, const std::vector<std::vector<std::vector<double>>>& Bslope,
  // std::vector<std::vector<std::vector<double>>>& Ashadow, const std::vector<std::vector<std::vector<double>>>& Bshadow, 
  // std::vector<std::vector<std::vector<double>>>& Ashadow_slope, const std::vector<std::vector<std::vector<double>>>& Bshadow_slope, 
  // std::vector<double>& Ashadow_info, const std::vector<double>& Bshadow_info,const std::vector<double>& partitionSize, unsigned const& ndep,T const& bndB)
  // const;



//   template <typename PUNIV>
//   void _asymDCdecNTC
//   ( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f,
//     double const& zopt, bool const cvx_ccv )
//   const;
//   template <typename PUNIV>
//   void _asymDCdecNTC
//   ( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f,
//     double const& zopt, bool const cvx_ccv, T const& bnd )
//   const;


//  template <typename PUNIV>
//   void _asymNTCsym
//   ( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f)
//   const;
//   template <typename PUNIV>
//   void _asymNTCsym
//   ( std::vector<std::vector<T>>& mat, unsigned const& ndep, PUNIV const& f, T const& bnd )
//   const;


//   template <typename PUNIV, typename BNUIA>
//   void _asymNTCsym_slope  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<double>>>& slope, const std::vector<double>& partitionSize, unsigned const& ndep, PUNIV const& f, BNUIA const& fDerv, double const& zopt, bool const cvx )
//   const;

//   template <typename PUNIV , typename BNUIA>
//   void _asymNTCsym_slope  ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<double>>>& slope, const std::vector<double>& partitionSize, unsigned const& ndep, PUNIV const& f, BNUIA const& fDerv, double const& zopt, bool const cvx, T const& bnd )
//   const;



  // void _inv
  // ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  // const;
  // void _sqr
  // ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  // const;
  // void _sqrt
  // ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  // const;
  // void _exp
  // ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  // const;
  // void _log
  // ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  // const;
  // void _xlog
  // ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  // const;
  // void _sin
  // ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  // const;
  // void _cos
  // ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
  // const;

//   void _cos_slope
//   ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<double>>>& slope, const std::vector<double>& partitionSize, unsigned const& ndep )
//   const;
//   void _exp
//   ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<double>>>& slope, const std::vector<double>& partitionSize, unsigned const& ndep )
//   const;

//   void _log
//   ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<double>>>& slope, const std::vector<double>& partitionSize, unsigned const& ndep )
//   const;
//   void _sqrt
//   ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<double>>>& slope, const std::vector<double>& partitionSize, unsigned const& ndep )
//   const;
//   void _inv
//   ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<double>>>& slope, const std::vector<double>& partitionSize, unsigned const& ndep )
//   const;

//   void _sqr
//   ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<double>>>& slope, const std::vector<double>& partitionSize, unsigned const& ndep )
//   const;

//   void _tanh
//   ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<double>>>& slope, const std::vector<double>& partitionSize, unsigned const& ndep )
//   const;


//  void _pow
//   ( std::vector<std::vector<T>>& mat,std::vector<std::vector<std::vector<double>>>& slope, const std::vector<double>& partitionSize, int const&iexp, unsigned const& ndep )
//   const;


//   void _prod
//   ( std::vector<std::vector<T>> const& mat1, std::vector<std::vector<T>> const& mat2,
//     std::vector<std::vector<T>>& mat3, unsigned& ndep3 )
//   const;
  
//   void _intersect
//   ( std::vector<std::vector<T>>& mat, unsigned const& ndep, T _rangebnd )
//   const;

//   void _intersect
//   ( std::vector<std::vector<T>>& mat, std::vector<std::vector<std::vector<double>>>& slope, const std::vector<double>& partitionSize, unsigned const& ndep, T _rangebnd )
//   const;

//   void _tanh
//   ( std::vector<std::vector<T>>& mat, unsigned const& ndep )
//   const;
//   void _pow
//   ( std::vector<std::vector<T>>& mat, int const&iexp, unsigned const& ndep )
//   const;

  // @unsure: the algorithm for computing asymmetric over-/under-estimators is applicable for sin and cos 
  // when they are concave or convex on the range of the input ASM  


  // @unsure: there is a similar algorithm which can refining the multiplicative result, 
  // by bounding x_iy_j+x_jy_i in 0.5*(x_i^2+y_i^2)*[-1,1] + 0.5*(x_j^2+y_j^2)*[-1,1]

  // static std::ostream& _dispmat
  // ( std::vector<std::vector<T>> const& mat, unsigned const len=3, std::ostream& out=std::cout );
  
  std::ostream& _dispvar
  ( std::vector<std::vector<T>> const& mat, unsigned const& ndep, const int& opt=0,
    std::ostream& out=std::cout )
  const;

  std::ostream& _dispvar
  ( std::vector<UnivarPWL<T>> const& lst, unsigned const& ndep, const int& opt=0,
    std::ostream& out=std::cout )
  const;

  // template <typename PUNIV>
  // double _goldsect
  // ( const double xL, const double xU, PUNIV const& f, const double TOL, const unsigned MAXIT )
  // const;
  
  // template <typename PUNIV>
  // double _goldsect_iter
  // ( unsigned& iter, const double a, const double fa, const double b,
  //   const double fb, const double c, const double fc, PUNIV const& f,
  //   const double TOL, const unsigned MAXIT )
  // const;


};

template <typename T> inline
typename ASModel<T>::Options ASModel<T>::options;


#ifdef MC__ASM_DEBUG_EVAL_LOGGING
template <typename T> inline
std::size_t ASModel<T>::debugging_ctr;
#endif

template <typename T>
std::ostream& operator<<
( std::ostream& out, const ASModel<T>& mod)
{
  out << std::endl;
  out << "ASM settings:" << std::endl;
  out << "   " << "no. variables:  " << mod._nvar << std::endl;
  out << "   " << "no. partitions: " << mod._ndiv << std::endl;
  out << "   " << "variable bounds: " << std::endl;
  for( unsigned int i=0; i<mod._nvar; i++ )
    out << "       " << i << ": " << (mod._defvar[i]? mod._bndvar[i]: "-") << std::endl;
  out << std::endl;
  return out;
}

// template <typename T>
// inline
// std::ostream& ASModel<T>::_dispmat
// ( std::vector<std::vector<T>> const& mat, unsigned const len, std::ostream& out )
// {
//   if( mat.empty() ) return out;
//   unsigned irow = 0;
//   for( auto&& row : mat ){
//     if( row.empty() ) continue;
//     out << std::right << std::setw(5) << irow++ <<": ";
//     unsigned icol = 0;
//     for( auto&&el : row ){
//       if( icol && !(icol%len) ) out << std::endl << "       ";
//       out << std::setw(0) << el;
//       icol++;
//     }
//     out << std::endl;
//   }
//   return out;
// }

template <typename T>
inline
std::ostream& ASModel<T>::_dispvar
( std::vector<std::vector<T>> const& mat, unsigned const& ndep, const int& opt,
  std::ostream& out )
const
{
  if( ndep > 2 ) return out;
  assert( !mat.empty() );
  out << std::scientific << std::setprecision(5) << std::right;
  
  if( ndep == 1 ){
    for( unsigned int i=0; i<_nvar; i++ ){
      if( mat[i].empty() ) continue;
      assert( _defvar[i] );
      double l = Op<T>::l(_bndvar[i]);
      double h = Op<T>::diam(_bndvar[i]) / (double)_ndiv;
      for( unsigned int j=0; j<_ndiv; j++, l+=h ){
        if( opt == 0 ){
          out << std::setw(14) << l   << std::setw(14) << Op<T>::l(mat[i][j]) << std::endl;
          out << std::setw(14) << l   << std::setw(14) << Op<T>::u(mat[i][j]) << std::endl;
          out << std::setw(14) << l+h << std::setw(14) << Op<T>::u(mat[i][j]) << std::endl;
          out << std::setw(14) << l+h << std::setw(14) << Op<T>::l(mat[i][j]) << std::endl;
          out << std::setw(14) << l   << std::setw(14) << Op<T>::l(mat[i][j]) << std::endl;
          out << std::endl;
        }
        else if( opt > 0 ){
          out << std::setw(14) << l   << std::setw(14) << Op<T>::u(mat[i][j]) << std::endl;
          out << std::setw(14) << l+h << std::setw(14) << Op<T>::u(mat[i][j]) << std::endl;
          out << std::endl;
        }
        else{
          out << std::setw(14) << l   << std::setw(14) << Op<T>::l(mat[i][j]) << std::endl;
          out << std::setw(14) << l+h << std::setw(14) << Op<T>::l(mat[i][j]) << std::endl;
          out << std::endl;
        }
      }
      break;
    }
  }
  
  else if( ndep == 2 ){
    for( unsigned int i1=0; i1<_nvar; i1++ ){
      if( mat[i1].empty() ) continue;
      assert( _defvar[i1] );
      //std::cout << "var1: " << _bndvar[i1] << std::endl;
      double l1 = Op<T>::l(_bndvar[i1]);
      double h1 = Op<T>::diam(_bndvar[i1]) / (double)_ndiv;
      for( unsigned int j1=0; j1<_ndiv; j1++, l1+=h1 ){
        for( unsigned int i2=i1+1; i2<_nvar; i2++ ){
          if( mat[i2].empty() ) continue;
          assert( _defvar[i2] );
          //std::cout << "var2: " << _bndvar[i2] << std::endl;
          double l2 = Op<T>::l(_bndvar[i2]);
          double h2 = Op<T>::diam(_bndvar[i2]) / (double)_ndiv;
          for( unsigned int j2=0; j2<_ndiv; j2++, l2+=h2 ){
            if( opt == 0 ){
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2    << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2    << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::endl << std::endl;
            }
            else if( opt > 0 ){
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2    << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::u(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::endl << std::endl;            
            }
            else{
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1+h1 << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2+h2 << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::setw(14) << l1    << std::setw(14) << l2    << std::setw(14) << Op<T>::l(mat[i1][j1]+mat[i2][j2]) << std::endl;
              out << std::endl << std::endl;            
            }
          }
          break;
        }
      }
      break;
    }
  }

  return out;
}


template <typename T>
inline
std::ostream& ASModel<T>::_dispvar
(std::vector<UnivarPWL<T>> const& lst, unsigned const& ndep, const int& opt,
  std::ostream& out )
const
{
  if( ndep > 2 ) return out;
  assert( !lst.empty() );
  out << std::scientific << std::setprecision(5) << std::right;


  if( ndep == 1 ){
    for( unsigned int i=0; i<_nvar; i++ ){
      if( lst[i].empty() ) continue;

      assert( _defvar[i] );

      const UnivarPWLE <double> & underEstmtr = lst[i].undEst;
      const UnivarPWLE <double> & overEstmtr  = lst[i].oveEst;   

      if( opt == 0 ){
        underEstmtr.display(out);
        overEstmtr.display(out);          
        // out << std::setw(14) << l   << std::setw(14) << Op<T>::l(mat[i][j]) << std::endl;
        // out << std::setw(14) << l   << std::setw(14) << Op<T>::u(mat[i][j]) << std::endl;
        // out << std::setw(14) << l+h << std::setw(14) << Op<T>::u(mat[i][j]) << std::endl;
        // out << std::setw(14) << l+h << std::setw(14) << Op<T>::l(mat[i][j]) << std::endl;
        // out << std::setw(14) << l   << std::setw(14) << Op<T>::l(mat[i][j]) << std::endl;
        // out << std::endl;
      }
      else if( opt > 0 ){
        overEstmtr.display(out);
        // if(overEstmtr.isCst()){
        //   if(overEstmtr.empty()){
        //     ;
        //   }
        //   else{
        //     out << std::setw(14) << overEstmtr.first[1] + overEstmtr.first[0]  << std::setw(14) << overEstmtr.second[0]  << std::endl;
        //     out << std::setw(14) << overEstmtr.first[2] + overEstmtr.first[0]  << std::setw(14) << overEstmtr.second[0]  << std::endl;
        //   }
        // }
        // else{
        //   for( unsigned int j=1; j < overEstmtr.second.size(); j++){
        //     out << std::setw(14) << overEstmtr.first[j] + overEstmtr.first[0]  << std::setw(14) << overEstmtr.second[j] + overEstmtr.second[0]  << std::endl;
        //   }          
        // }
      }
      else{
        underEstmtr.display(out);
        // if (underEstmtr.isCst()){
        //   if(underEstmtr.empty()){
        //     ;
        //   }
        //   else{
        //     out << std::setw(14) << underEstmtr.first[1] + underEstmtr.first[0]  << std::setw(14) << underEstmtr.second[0]  << std::endl;
        //     out << std::setw(14) << underEstmtr.first[2] + underEstmtr.first[0]  << std::setw(14) << underEstmtr.second[0]  << std::endl;
        //   }
        // }
        // else{
        //   for( unsigned int j=1; j< underEstmtr.second.size(); j++ ){
        //     out << std::setw(14) << underEstmtr.first[j] + underEstmtr.first[0] << std::setw(14) << underEstmtr.second[j] + underEstmtr.second[0] << std::endl;        
        //   }
        // }
      }
      break; // break the outer loop
    }
  }
  else if( ndep == 2 ){
    std::vector<unsigned int> indices;
    indices.reserve(2);     
    for( unsigned int i1=0; i1<_nvar; i1++ ){
      if( lst[i1].empty() ) continue;
      assert( _defvar[i1] );
      indices.push_back(i1);
      //std::cout << "var1: " << _bndvar[i1] << std::endl;
    }
    const UnivarPWLE <double> & underEstmtr = lst[indices[0]].undEst;
    const UnivarPWLE <double> & overEstmtr  = lst[indices[0]].oveEst;   

    if( opt == 0 ){
      underEstmtr.display(out,lst[indices[1]].undEst);
       overEstmtr.display(out,lst[indices[1]].oveEst);          
    }
    else if( opt > 0 ){
      overEstmtr.display(out,lst[indices[1]].oveEst);
    }
    else{
      underEstmtr.display(out,lst[indices[1]].undEst);
    }

  }
  return out;
}






template <typename T>
inline
void
ASModel<T>::_asym_relu
( std::vector<UnivarPWL<T>>& lst, unsigned const& ndep)
const
{

  assert( !lst.empty() );


  T bnd = _B( lst, 1 );

  // anchor points
  double lambda ( Op<T>::l(bnd) );
  double mu ( Op<T>::u(bnd) );
  double sum_r1( 0 );
  for( unsigned int i=0; i<_nvar; i++ ){
    if( lst[i].empty() ) continue;
    _U2[i] = lst[i].oveEst.get_lb();
    //_U1[i] = lst[i].oveEst.get_ub();  
    _r1[i] = ( _U1[i] - _U2[i] ); sum_r1+=_r1[i];
  }

      
  //UnivarPWL<T> _tmp(0.);
  double theta_i_times_mu = mu/ndep;
  for( unsigned int i=0; i<_nvar; i++ ){
    if( lst[i].empty() ) continue;   
    const double rowOffsetUnder = - _L1[i] + lambda;
    lst[i].undEst = relu( lst[i].undEst + rowOffsetUnder );

    if(sum_r1 != 0.){
      const double theta_i = _r1[i] / sum_r1;//(_c2[i] - _U2[i])/(C2 - sigma_u);
      theta_i_times_mu = theta_i*mu;
    }
    const double rowOffsetOver =  - _U1[i] + theta_i_times_mu;
    //std::cout << "i " << i << " rowOffsetOver " << rowOffsetOver <<std::endl;
    lst[i].oveEst = relu( lst[i].oveEst + rowOffsetOver );
  }
}



template <typename T>
inline
void
ASModel<T>::_asym_relu_withShadow_old
( std::vector<UnivarPWL<T>>& lst,unsigned const& ndep,std::vector<UnivarPWL<T>>& shadow)
const
{
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "relu_ws" << std::endl;
#endif
  assert( !lst.empty() );

  // Get the maxima of all components of the input underestimator, stored in _L2[i]
  double sigma_o = 0.;       // <- the maximum of the input underestimator
  for( unsigned int i=0; i<_nvar; i++ ){
    if( lst[i].empty() ) continue;
     _L2[i]= lst[i].undEst.get_ub();
     sigma_o += _L2[i];
  }

  T bnd = _B( lst, 1 );

  // anchor points
  double lambda ( Op<T>::l(bnd) );
  double mu ( Op<T>::u(bnd) );
  double sum_r1( mu - lambda );
  double C1( 0. ), C2(0.);
  for( unsigned int i=0; i<_nvar; i++ ){
    if( lst[i].empty() ) continue;
    _c1[i] = _L1[i]; C1 += _c1[i];
    _c2[i] = _U1[i]; C2 += _c2[i];     
    _r1[i] = ( _U1[i] - _L1[i] );
  }
  if(std::fabs(C1 - lambda) > 5e5*MC__ASM_COMPUTATION_TOL || std::fabs(C2 - mu) > 5e5*MC__ASM_COMPUTATION_TOL){
    std::cout << "numerical error in asym_slope_relu ws" << std::endl;
    std::cout << std::setprecision(18) << std::fabs(C1 - lambda) - MC__ASM_COMPUTATION_TOL << std::endl;
    std::cout << std::setprecision(18) << std::fabs(C2 - mu) - MC__ASM_COMPUTATION_TOL << std::endl;
    throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); 
  }
  else if(std::fabs(C1 - lambda) > 1e2*MC__ASM_COMPUTATION_TOL || std::fabs(C2 - mu) > 1e2*MC__ASM_COMPUTATION_TOL){
    C1 = std::min(C1,lambda);
    lambda = C1;
    C2 = std::max(C2,mu);
    mu = C1;
  }  
      
  const double shadow_global_offset = (1.0 - 1.0/((double) ndep))*sigma_o;
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "    shadow_global_offset = " << shadow_global_offset << std::endl;
#endif

  UnivarPWL<T> _tmp(0.);
  for( unsigned int i=0; i<_nvar; i++ ){
    if( lst[i].empty() ) continue;   
    if( isequal( _r1[i], 0. ) ){     
      lst[i] = _tmp;
      shadow[i] = _tmp;
      continue;
    }
    else{
      const double rowOffsetUnder = - _c1[i] + C1;    // - Li + lambda

      const double theta_i = _r1[i] / sum_r1;//(_c2[i] - _U2[i])/(C2 - sigma_u);
      const double theta_i_times_mu = theta_i*C2;
      const double rowOffsetOver =  - _c2[i] + theta_i_times_mu;  // -Ui + theta_i * mu
      lst[i].oveEst = relu( lst[i].oveEst + rowOffsetOver );     

      const double rowOffsetShadow = - _L2[i] + sigma_o;  // -sigma_i + sigma_o
      // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
      shadow[i].undEst = relu( lst[i].undEst + rowOffsetShadow ) - std::max(shadow_global_offset,0.);
         lst[i].undEst = relu( lst[i].undEst + rowOffsetUnder );
      shadow[i].undEst -= lst[i].undEst;


      // for( unsigned int j=0; j<_ndiv; j++ ){

      //     // compute shadow
      //     // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
      //     const double delta_l = std::fabs(slope[i][j][0]*partitionSize[i]); 
      //     const double loPtShadow = zL + rowOffsetShadow;//zL - _c2[i] + C2;
      //     const double upPtShadow = loPtShadow + delta_l;//zL - _c2[i] + C2;        
      //     const double shadowLoBnd = f( loPtShadow ) - El;// by monotonicity
      //     //const double shadowUpBnd = f( upPtShadow ) - f( loPtUnderEstmt + delta_l);
      //     shadow[0][i][j] = shadowLoBnd - f(shadow_global_offset);
      //     //shadow[2][i][j] = shadowUpBnd - f(shadow_global_offset);
            
      //     //if(shadowUpBnd > shadowLoBnd + 1e2*MC__ASM_COMPUTATION_TOL ){
      //       // The update of the derv needs special consideration
      //       const double upBndUnderEstmt = f(loPtUnderEstmt + delta_l);
      //       double slope0shadow_to_be_multiplied =  0.;//std::min(shadowFuncDervAtRightPt, shadowFuncDervAtLeftPt);
      //       if (shadowLoBnd == 0.)//f(loPtUnderEstmtShadow) < 1e2*MC__ASM_COMPUTATION_TOL)    // case 1 slope must be 0
      //         slope0shadow_to_be_multiplied = 0.;
      //       else if (f(loPtUnderEstmt+delta_l + 1e5*MC__ASM_COMPUTATION_TOL) == 0. ){         // case 2 slope shoud be 1
      //         slope0shadow_to_be_multiplied = 1.0;//std::min(fDerv(loPtShadow) - fDerv(loPtUnderEstmt),fDerv(upPtShadow) - fDerv(loPtUnderEstmt + delta_l));
      //         //std::cout << f(loPtUnderEstmt+delta_l + 1e2*MC__ASM_COMPUTATION_TOL) << "," << shadowFuncValueAtLeftPt << std::endl;
      //       }
      //       else if (f(loPtUnderEstmt-1e2*MC__ASM_COMPUTATION_TOL) == 0.){                    // case 3 concave part 
      //         if ( delta_l > 1e2 * MC__ASM_COMPUTATION_TOL)
      //           slope0shadow_to_be_multiplied = (f(upPtShadow) - upBndUnderEstmt - (f(loPtShadow) -f(loPtUnderEstmt)))/(delta_l);
      //         else 
      //           ;//slope0shadow_to_be_multiplied = fDerv(upPtShadow) - fDerv(loPtUnderEstmt + delta_l);
      //       }
      //       shadow_slope[0][i][j] = slope0shadow_to_be_multiplied * slope[i][j][0];
     

    }
  }

}


template <typename T>
inline
void
ASModel<T>::_asym_relu_withShadow
( std::vector<UnivarPWL<T>>& lst,unsigned const& ndep,std::vector<UnivarPWL<T>>& shadow, std::vector<double>& shadow_info)
const
{
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "    _asym_relu_withShadow" << std::endl;
  // for( unsigned int i=0; i<_nvar; i++ ){
  //   if( lst[i].empty() ) continue;
  //   if(shadow_info[1]>0) std::cout << "shadowO " << i << ": "<< shadow[i].oveEst << std::endl;
  //   std::cout << "lstO " << i << ": "<< lst[i].oveEst << std::endl;   
  //   if(shadow_info[0]>0) std::cout << "shadowU " << i << ": "<< shadow[i].undEst << std::endl;
  //   std::cout << "lstU " << i << ": "<< lst[i].undEst << std::endl;         
  // }  
#endif

  assert( !lst.empty() );

  
 
 // Main idea: shadow-enhanced abstract superposition models (SeASM) are of the form
 //   ove: = min\{z_1^o,z_2^o,\ldots,z_{n_o}^o\} where each z_k^o = \sum_{i=1}^{nvar} z_{k,i}^o(x_i) 
 //   und: = max\{z_1^u,z_2^u,\ldots,z_{n_u}^u\} where each z_k^u = \sum_{i=1}^{nvar} z_{k,i}^u(x_i) 
 //   We also let the z_1^o has the least maximum, and z_1^u has the greatest minimum 
 // Hence, the composition of SeASMs through relu is given by (the monotonicity)
 //   ove: =  min\{relu(z_1^o),relu(z_2^o),\ldots,relu(z_{n_o}^o)\}
 //   und: =  max\{relu(z_1^u),relu(z_2^u),\ldots,relu(z_{n_u}^o)\}

 // As the first draft, we only fix the maximum of shadow estimators to be 1


//   if(shadow_info.size()<2){
// #ifndef MC__ASMODEL_NOT_DEBUG_SHADOW 
//     std::cout << "    shadow_info.size()<2" << std::endl;
// #endif
//     _shadowInfo_init(shadow_info);
//   }

  // Step 1: make sure if the input var is greater than 0
  //         To this end, we need to get the lb of the var and the lb of the shadow underestimators (if any)
  double lambda(-DBL_MAX),mu(DBL_MAX);
  if(_IntmdtCntnrSeted){
    for( unsigned int i=0; i<_nvar; i++ ){
      if( lst[i].empty() ) continue;
      lambda += _L1[i];
      mu     += _U1[i];     
    }
    _IntmdtCntnrSeted = false;
  }
  else{
    T bnd  = _B( lst, 1 );
    lambda = Op<T>::l(bnd);
    mu     = Op<T>::u(bnd);
    _IntmdtCntnrSeted = false;
    // // To indicate numerical issues
    // double C1( 0. ), C2(0.);
    // for( unsigned int i=0; i<_nvar; i++ ){
    //   if( lst[i].empty() ) continue;
    //   _c1[i] = _L1[i]; C1 += _c1[i];
    //   _c2[i] = _U1[i]; C2 += _c2[i];     
    //   _r1[i] = ( _U1[i] - _L1[i] );
    // }
  
    // 
    // if(std::fabs(C1 - lambda) > 5e5*MC__ASM_COMPUTATION_TOL || std::fabs(C2 - mu) > 5e5*MC__ASM_COMPUTATION_TOL){
    //   std::cout << "numerical error in asym_slope_relu ws" << std::endl;
    //   std::cout << std::setprecision(18) << std::fabs(C1 - lambda) - MC__ASM_COMPUTATION_TOL << std::endl;
    //   std::cout << std::setprecision(18) << std::fabs(C2 - mu) - MC__ASM_COMPUTATION_TOL << std::endl;
    //   throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); 
    // }
    // else if(std::fabs(C1 - lambda) > 1e2*MC__ASM_COMPUTATION_TOL || std::fabs(C2 - mu) > 1e2*MC__ASM_COMPUTATION_TOL){
    //   C1 = std::min(C1,lambda);
    //   lambda = C1;
    //   C2 = std::max(C2,mu);
    //   mu = C1;
    // }  

  }
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "        Step 1: to see if relu doesn't truncate" << std::endl;
#endif

  if(shadow_info[0] == 0 && lambda > -MC__ASM_COMPUTATION_TOL){ // if there is no shadow underestimator, then we do nothing with the var
    // It would be benificial for debugging to check here if 
    // the lowerbounds of the overestimators are strictly greater than 0
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "            Step 1: relu doesn't truncate ACT" << std::endl;
#endif
    return ;
  }    
  else if(lambda > -MC__ASM_COMPUTATION_TOL){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "            Step 1: UND has SHADOW" << std::endl;
#endif    
    // At this time, we need to check if the lb of the shadow underestimator is greater than 0 as well, o.w.,
    // we need to get the shadow underestimator of the shadow underestimator

    // Get the lb of the shadow underestimator
    double lbShadowUnd = 0.;       
    for( unsigned int i=0; i<_nvar; i++ ){
      if( lst[i].empty() ) continue;
       lbShadowUnd += shadow[i].undEst.get_lb();
    }      
    if (lbShadowUnd > -MC__ASM_COMPUTATION_TOL){
      // It would be benificial for debugging to check here if 
      // the lowerbounds of the overestimators are strictly greater than 0     
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "            Step 1: relu doesn't truncate SHA" << std::endl;
#endif       
      return ; // if the lb of the shadow underestimator > 0, then we do nothing with the var
    }
    else{   // if the lb of the shadow underestimator < 0, then we take the shadow of it
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "            Step 1: relu truncate SHA" << std::endl;
#endif
      // Get the maximum of all components of shadow underestimator, stored in _L2[i]
      // This is to construct the shadow underestimator resulted from the shadow underestimator
      double sigma_o = 0.;       // <- the maximum of the shadow underestimator
      for( unsigned int i=0; i<_nvar; i++ ){
        if( lst[i].empty() ) continue;
         _L2[i]= shadow[i].undEst.get_ub();
         sigma_o += _L2[i];   
      }

      if (sigma_o < MC__ASM_COMPUTATION_TOL){  // if the shadow underestimator is useless
        for( unsigned int i=0; i<_nvar; i++ ){
          if( lst[i].empty() ) continue;
          shadow[i].undEst.clear();  // note that the label _under is flagged as true by default 
        }       
        shadow_info[0] = 0;
        return;
      }

      const double shadow_global_offset = (1.0 - 1.0/((double) ndep))*sigma_o;
      for( unsigned int i=0; i<_nvar; i++ ){
        if( lst[i].empty() ) continue;   
        else{
          const double rowOffsetShadow = - _L2[i] + sigma_o;  // -sigma_i + sigma_o
          // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
          shadow[i].undEst = relu( lst[i].undEst + rowOffsetShadow ) - shadow_global_offset;
          shadow[i].flag_nonEmpty();
        }
      }
      return ;
    }
  }
  


  // Step 2: since the lb of the var < 0, we process both the active and the shadow estimators
  //         To this end, we need to get the ub of the active UNDerestimator
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "        Step 2: to see how relu truncates" << std::endl;
#endif
  double sigma_o = 0.;       // <- the maximum of the active UNDerestimator
  for( unsigned int i=0; i<_nvar; i++ ){
    if( lst[i].empty() ) continue;
     _L2[i]= lst[i].undEst.get_ub();
     sigma_o += _L2[i];
  }      

  // Step 2.1: when there is no shadow UNDerestimator
  if(shadow_info[0] == 0){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "            Step 2: NO SHA " << std::endl;
#endif    
    if(sigma_o < 0.){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                Step 2: ubACT < 0" << std::endl;
#endif          
     for( unsigned int i=0; i<_nvar; i++ ){
      if( lst[i].empty() ) continue;
      lst[i].undEst.set_zero();
     }       
    }
    else{
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                Step 2: ubACT > 0" << std::endl;
#endif           
      const double shadow_global_offset = (1.0 - 1.0/((double) ndep))*sigma_o;
      
    
      //  UnivarPWL<T> _tmp(0.);
      for( unsigned int i=0; i<_nvar; i++ ){
        if( lst[i].empty() ) continue;   
        // if( isequal( _r1[i], 0. ) ){     
        //   lst[i] = _tmp;
        //   continue;
        // } // since the overestimator and the underestimator are not aligned, we cannot do this.

        const double rowOffsetUnder  = - _L1[i] + lambda;    // - Li + lambda
        const double rowOffsetShadow = - _L2[i] + sigma_o;  // -sigma_i + sigma_o
        // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
//        std::cout << lst[i].undEst + rowOffsetShadow  << std::endl;
        shadow[i].undEst = relu( lst[i].undEst + rowOffsetShadow ) - std::max(shadow_global_offset,0.);
           lst[i].undEst = relu( lst[i].undEst + rowOffsetUnder );
        shadow[i].flag_nonEmpty();
        //shadow[i].undEst -= lst[i].undEst;
//        std::cout << shadow[i].undEst  << std::endl;
      }
      shadow_info[0] = 1;
    }

    // update active OVErestimator
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                Step 2: update ACT OVE" << std::endl;
#endif     
    // Get the minima of all components of the input OVErestimator, stored in _U2[i]
    double sigma_u = 0.;       // <- the minimum of the active OVErestimator
    for( unsigned int i=0; i<_nvar; i++ ){
      if( lst[i].empty() ) continue;
       _U2[i]= lst[i].oveEst.get_lb();
       sigma_u += _U2[i];
       _r2[i] = ( _U1[i] - _U2[i] );
    }        
    double sum_r2( mu - sigma_u );       
    double theta_i_times_mu = mu/((double)ndep);
    for( unsigned int i=0; i<_nvar; i++ ){
      if( lst[i].empty() ) continue;   
      // if( isequal( _r1[i], 0. ) ){     
      //   lst[i] = _tmp;
      //   continue;
      // } // since the overestimator and the underestimator are not aligned, we cannot do this.
      if(sum_r2 != 0.){
        const double theta_i = _r2[i] / sum_r2;//(_c2[i] - _U2[i])/(C2 - sigma_u);
        theta_i_times_mu = theta_i*mu;
      }
      const double rowOffsetOver =  - _U1[i] + theta_i_times_mu;  // -Ui + theta_i * mu
      lst[i].oveEst = relu( lst[i].oveEst + rowOffsetOver );      // it should be note that we leave to relu(PWL) the determination of directly passing, 
                                                                  // instead of checking whether sigma_u > 0 here
    }

    // update shadow OVErestimator if any
    if (shadow_info[1] > 0){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                Step 2: update SHA OVE" << std::endl;
#endif           
      // shadow_info[1] = 0;
      double sigma_u = 0.;       // <- the minimum of the shadow overestimator
      mu = 0.;
      for( unsigned int i=0; i<_nvar; i++ ){
        if( lst[i].empty() ) continue;
        _U1[i]= shadow[i].oveEst.get_ub();
        _U2[i]= shadow[i].oveEst.get_lb();
        mu += _U1[i];
        sigma_u += _U2[i];
        _r2[i] = ( _U1[i] - _U2[i] );
      }        
      double sum_r2( mu - sigma_u ); 
      double theta_i_times_mu = mu/((double)ndep);
      for( unsigned int i=0; i<_nvar; i++ ){
        if( lst[i].empty() ) continue;   
        if(sum_r2 != 0.){
          const double theta_i = _r2[i] / sum_r2;//(_c2[i] - _U2[i])/(C2 - sigma_u);
          theta_i_times_mu = theta_i*mu;
        }
        const double rowOffsetOver =  - _U1[i] + theta_i_times_mu;  // -Ui + theta_i * mu
        //std::cout << "            I OVE \n" << shadow[i].oveEst << std::endl;
        shadow[i].oveEst = relu( shadow[i].oveEst + rowOffsetOver );      // it should be note that we leave to relu(PWL) the determination of directly passing, 
        //std::cout << "            I OVE \n" << shadow[i].oveEst << std::endl;
        shadow[i].flag_nonEmpty();
      }
    }   
    return ;    
  } // end if shadow_info[0] == 0 
  else{
  // Step 2.2: when there is a shadow UNDerestimator
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "            Step 2: with SHA " << std::endl;
#endif     
    // Step 2.2.1: 
    if(sigma_o <= MC__ASM_COMPUTATION_TOL){ // eliminate the active underestimator as it is truncated
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                Step 2: ACT ubUND <=0 " << std::endl;
#endif   

      for( unsigned int i=0; i<_nvar; i++ ){
        if( lst[i].empty() ) continue;
        lst[i].undEst.set_zero();
      }       
            
      // Recall the criteria of ACT and SHA: (1) lbACT >= lbSHA, (2) ubACT <= ubSHA, and (3) ubACT > ubSHA if lbACT = lbSHA,
      // so we know that both the ACT-SHA and the SHA-ACT must be zero, and then we only need to see whether
      // SHA-SHA is not zero 

      double sigma_o = 0.;       // <- the maximum of the SHA UNDerestimator
      for( unsigned int i=0; i<_nvar; i++ ){
        if( lst[i].empty() ) continue;
        _L2[i]= shadow[i].undEst.get_ub();
        sigma_o += _L2[i];
      }      
    
      if(sigma_o <= 0.){ // eliminate the shadow underestimator as it is useless
        for( unsigned int i=0; i<_nvar; i++ ){
          if( lst[i].empty() ) continue;
          shadow[i].undEst.clear();
        }
        shadow_info[0] = 0;          
      }
      else{  // update the shadow underestimator
        const double shadow_global_offset = std::max((1.0 - 1.0/((double) ndep))*sigma_o,0.);
        for( unsigned int i=0; i<_nvar; i++ ){
          if( lst[i].empty() ) continue;   
          const double rowOffsetShadow = - _L2[i] + sigma_o;  // -sigma_i + sigma_o
          // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
          shadow[i].undEst = relu( shadow[i].undEst + rowOffsetShadow ) - shadow_global_offset;
          shadow[i].flag_nonEmpty();
        }
      }
    }
    else{   
    // Step 2.2.2: 
    // in worst case we may need to find the best ACT and SHA from four new underestimators
    // we may also need to eliminate redundant SHA if that is strictly looser than the ACT
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                Step 2: ACT ubUND > 0 " << std::endl;
#endif 
      double sigma_oSHA = 0.;       // <- the maximum of the SHA UNDerestimator
      std::vector<double> L2(_nvar,0.);
      for( unsigned int i=0; i<_nvar; i++ ){
        if( lst[i].empty() ) continue;
        L2[i]= shadow[i].undEst.get_ub();
        sigma_oSHA += L2[i];
      }      

      // std::cout << "                Step 2: ACT ubUND > 0 sigma_oSHA = " << std::scientific << std::setprecision(14) << sigma_oSHA << std::endl;    
      // std::cout << "                Step 2: ACT ubUND > 0 sigma_o = " << std::scientific << std::setprecision(14) << sigma_o << std::endl;    

      if(sigma_oSHA > sigma_o + MC__ASM_COMPUTATION_TOL ){ // in this case, ACT-ACT with SHA-SHA, (at this time the lbSHA < lbACT) 
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                    Step 2: sigma_oSHA > sigma_o" << std::endl;
#endif 

        const double shadow_global_offset = (1.0 - 1.0/((double) ndep))*sigma_oSHA;

        for( unsigned int i=0; i<_nvar; i++ ){
          if( lst[i].empty() ) continue;
          const double rowOffsetUnder  = -_L1[i] + lambda;      // - Li + lambda
          const double rowOffsetShadow = - L2[i] + sigma_oSHA;  // -sigma_i + sigma_o
          // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
          shadow[i].undEst = relu( shadow[i].undEst + rowOffsetShadow ) - std::max(shadow_global_offset,0.);
             lst[i].undEst = relu(    lst[i].undEst + rowOffsetUnder );
        }        
      }
      else if(sigma_o > sigma_oSHA + MC__ASM_COMPUTATION_TOL){  // in this case ACT-ACT with ACT-SHA
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                    Step 2: sigma_o > sigma_oSHA" << std::endl;
#endif 
        const double shadow_global_offset = (1.0 - 1.0/((double) ndep))*sigma_o;

        for( unsigned int i=0; i<_nvar; i++ ){
          if( lst[i].empty() ) continue;
          const double rowOffsetUnder  = -_L1[i] + lambda;    // - Li + lambda
          const double rowOffsetShadow = -_L2[i] + sigma_o;   // -sigma_i + sigma_o
          // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
          shadow[i].undEst = relu( lst[i].undEst + rowOffsetShadow) - std::max(shadow_global_offset,0.);
             lst[i].undEst = relu( lst[i].undEst + rowOffsetUnder );
          shadow[i].flag_nonEmpty();   
        }  
      }
      else{ // in this case, we have either lbACT > lbSHA or lbACT = lbSHA. 
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                    Step 2: sigma_o approx sigma_oSHA" << std::endl;
#endif         
        // If lbACT > lbSHA, the new ACT is from ACT but the new SHA may from either SHA-SHA or ACT-SHA
        // We choose the NEW SHA with GREATER lb (by using the monotonicity of relu)

        double lambda_SHA = 0.;       // <- the minimum of the SHA UNDerestimator
        std::vector<double> L1(_nvar,0.);

        for( unsigned int i=0; i<_nvar; i++ ){
          if( lst[i].empty() ) continue;
          L1[i]= shadow[i].undEst.get_lb();
          lambda_SHA += L1[i];
        }     
        
        // std::cout << "                    Step 2: sigma_o approx sigma_oSHA lambda_SHA = " << std::scientific << std::setprecision(14) << lambda_SHA << std::endl;    
        // std::cout << "                    Step 2: sigma_o approx sigma_oSHA lambda = " << std::scientific << std::setprecision(14) << lambda << std::endl;    

        if(lambda > lambda_SHA - 5e-14){          
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                    Step 2: and lambda > lambda_SHA" << std::endl;
#endif                         
          double lbSHASHA = -((double ) ndep-1.)*sigma_oSHA;
          double lbACTSHA = -((double ) ndep-1.)*sigma_o;

          // std::cout << "                        lbSHASHA " << std::scientific << std::setprecision(14) << lbSHASHA << std::endl;  
          // std::cout << "                        lbACTSHA " << std::scientific << std::setprecision(14) << lbACTSHA << std::endl;  


          for( unsigned int i=0; i<_nvar; i++ ){
            if( lst[i].empty() ) continue;
            lbACTSHA += std::max(0., _L1[i] - _L2[i] + sigma_o);
            lbSHASHA += std::max(0.,  L1[i] -  L2[i] + sigma_oSHA);
          }     
          if(lbSHASHA >= lbACTSHA - 5e-14){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                    Step 2: and and lbSHASHA >= lbACTSHA " << std::scientific << std::setprecision(14) << lbSHASHA - lbACTSHA << std::endl;
  std::cout << "                        lbSHASHA " << std::scientific << std::setprecision(14) << lbSHASHA << std::endl;  
  std::cout << "                        lbACTSHA " << std::scientific << std::setprecision(14) << lbACTSHA << std::endl;  
#endif          
            const double shadow_global_offset = (1.0 - 1.0/((double) ndep))*sigma_oSHA;
    
            for( unsigned int i=0; i<_nvar; i++ ){
              if( lst[i].empty() ) continue;
              const double rowOffsetUnder  = -_L1[i] + lambda;    // - Li + lambda
              const double rowOffsetShadow = - L2[i] + sigma_oSHA;   // -sigma_i + sigma_o
              // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
              shadow[i].undEst = relu( shadow[i].undEst + rowOffsetShadow) - std::max(shadow_global_offset,0.);
                 lst[i].undEst = relu(    lst[i].undEst + rowOffsetUnder );
            }  
          }
          else{
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                    Step 2: and and lbSHASHA < lbACTSHA" << std::scientific << std::setprecision(14) << lbSHASHA - lbACTSHA << std::endl;
  std::cout << "                        lbSHASHA " << std::scientific << std::setprecision(14) << lbSHASHA << std::endl;  
  std::cout << "                        lbACTSHA " << std::scientific << std::setprecision(14) << lbACTSHA << std::endl;  
#endif                      
            const double shadow_global_offset = (1.0 - 1.0/((double) ndep))*sigma_o;
    
            for( unsigned int i=0; i<_nvar; i++ ){
              if( lst[i].empty() ) continue;
              const double rowOffsetUnder  = -_L1[i] + lambda;    // - Li + lambda
              const double rowOffsetShadow = -_L2[i] + sigma_o;   // -sigma_i + sigma_o
              // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
              shadow[i].undEst = relu( lst[i].undEst + rowOffsetShadow) - std::max(shadow_global_offset,0.);
                 lst[i].undEst = relu( lst[i].undEst + rowOffsetUnder );
              shadow[i].flag_nonEmpty();   
            }  
          }          
        }
        else{ // If lbACT = lbSHA, then ACT-ACT or SHA-ACT and ACT-SHA or SHA-SHA
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                    Step 2: and lambda < lambda_SHA" << std::endl;
#endif  
          double ubACTACT(0.);
          double ubSHAACT(0.);    
          double lbSHASHA = -((double ) ndep-1.)*sigma_oSHA;
          double lbACTSHA = -((double ) ndep-1.)*sigma_o;
          for( unsigned int i=0; i<_nvar; i++ ){
            if( lst[i].empty() ) continue;
            lbACTSHA = lbACTSHA + std::max(0., _L1[i] - _L2[i] + sigma_o);
            lbSHASHA = lbSHASHA + std::max(0.,  L1[i] -  L2[i] + sigma_oSHA);
            ubACTACT = ubACTACT + std::max(0., _L2[i] - _L1[i] + lambda);
            ubSHAACT = ubSHAACT + std::max(0.,  L2[i] -  L1[i] + lambda_SHA);
          }

          //const double shadow_global_offsetACT = (1.0 - 1.0/((double) ndep))*sigma_o;
          //const double shadow_global_offsetSHA = (1.0 - 1.0/((double) ndep))*sigma_oSHA;     

         

          if(lbSHASHA >= lbACTSHA - 5e-14 && ubACTACT >= ubSHAACT - 5e-14 ){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                    Step 2: and and lbSHASHA >= lbACTSHA && ubACTACT >= ubSHAACT" << std::endl;
#endif 
            const double shadow_global_offset = (1.0 - 1.0/((double) ndep))*sigma_oSHA;
    
            for( unsigned int i=0; i<_nvar; i++ ){
              if( lst[i].empty() ) continue;
              const double rowOffsetUnder  = -_L1[i] + lambda;    // - Li + lambda
              const double rowOffsetShadow = - L2[i] + sigma_oSHA;   // -sigma_i + sigma_o
              // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
              shadow[i].undEst = relu( shadow[i].undEst + rowOffsetShadow) - std::max(shadow_global_offset,0.);
                 lst[i].undEst = relu(    lst[i].undEst + rowOffsetUnder );
            }  
          }
          else if(lbSHASHA >= lbACTSHA - 5e-14){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                    Step 2: and and lbSHASHA >= lbACTSHA && ubACTACT < ubSHAACT" << std::endl;
#endif             
            const double shadow_global_offset = (1.0 - 1.0/((double) ndep))*sigma_oSHA;
    
            for( unsigned int i=0; i<_nvar; i++ ){
              if( lst[i].empty() ) continue;
              const double rowOffsetUnder  = - L1[i] + lambda_SHA;   // - Li + lambda
              const double rowOffsetShadow = - L2[i] + sigma_oSHA;   // -sigma_i + sigma_o
              // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
                 lst[i].undEst = relu( shadow[i].undEst + rowOffsetUnder );
              shadow[i].undEst = relu( shadow[i].undEst + rowOffsetShadow) - std::max(shadow_global_offset,0.);  
            }  
          }
          else if(ubACTACT >= ubSHAACT - 5e-14){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                    Step 2: and and lbSHASHA < lbACTSHA && ubACTACT >= ubSHAACT" << std::endl;
#endif             
            const double shadow_global_offset = (1.0 - 1.0/((double) ndep))*sigma_o;
    
            for( unsigned int i=0; i<_nvar; i++ ){
              if( lst[i].empty() ) continue;
              const double rowOffsetUnder  = -_L1[i] + lambda;    // - Li + lambda
              const double rowOffsetShadow = -_L2[i] + sigma_oSHA;   // -sigma_i + sigma_o
              // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing
              shadow[i].undEst = relu( lst[i].undEst + rowOffsetShadow) - std::max(shadow_global_offset,0.);
                 lst[i].undEst = relu( lst[i].undEst + rowOffsetUnder );
              shadow[i].flag_nonEmpty();   
            }  
          }
          else{
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                    Step 2: and and lbSHASHA < lbACTSHA && ubACTACT < ubSHAACT" << std::endl;
#endif            
            const double shadow_global_offset = (1.0 - 1.0/((double) ndep))*sigma_o;
    
            for( unsigned int i=0; i<_nvar; i++ ){
              if( lst[i].empty() ) continue;
              const double rowOffsetUnder  = - L1[i] + lambda_SHA;    // - Li + lambda
              const double rowOffsetShadow = -_L2[i] + sigma_o;       // -sigma_i + sigma_o
              // Note that the function max(zi - zi^U + mu) - max(zi - zi^L + lambda) is no longer convex but still nondecreasing              
                 lst[i].undEst = relu(    lst[i].undEst + rowOffsetShadow) - std::max(shadow_global_offset,0.);
              shadow[i].undEst = relu( shadow[i].undEst + rowOffsetUnder );
              std::swap(shadow[i].undEst,lst[i].undEst);
            }  
          }          

        }
            
      } // end the case of either lbACT > lbSHA or lbACT = lbSHA. 
    
    }

    // update both the active and shadow overestimator (if any)
    // update the actOve 
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                Step 2: update ACT OVE " << std::endl;
#endif         
    double sigma_u = 0.;       // <- the minimum of the active overestimator
    for( unsigned int i=0; i<_nvar; i++ ){
      if( lst[i].empty() ) continue;
      // _U1[i]= lst[i].oveEst.get_ub(); as it has been set. 
      _U2[i]= lst[i].oveEst.get_lb();
      sigma_u += _U2[i];
      _r2[i] = ( _U1[i] - _U2[i] );
    }        
    double sum_r2( mu - sigma_u ); 
    double theta_i_times_mu = mu/((double) ndep);
    for( unsigned int i=0; i<_nvar; i++ ){
      if( lst[i].empty() ) continue;   
      if(sum_r2 != 0.){
        const double theta_i = _r2[i] / sum_r2;//(_c2[i] - _U2[i])/(C2 - sigma_u);
        theta_i_times_mu = theta_i*mu;
      }
      const double rowOffsetOver =  - _U1[i] + theta_i_times_mu;  // -Ui + theta_i * mu
      lst[i].oveEst = relu( lst[i].oveEst + rowOffsetOver );      // it should be note that we leave to relu(PWL) the determination of directly passing, 
    }

      
    // update the shaOve 
    if (shadow_info[1] > 0){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW
  std::cout << "                Step 2: update SHA OVE " << std::endl;
#endif               
      double sigma_u = 0.;       // <- the minimum of the shadow overestimator
      mu = 0.;
      for( unsigned int i=0; i<_nvar; i++ ){
        if( lst[i].empty() ) continue;
        _U1[i]= shadow[i].oveEst.get_ub();
        _U2[i]= shadow[i].oveEst.get_lb();
        mu += _U1[i];
        sigma_u += _U2[i];
        _r2[i] = ( _U1[i] - _U2[i] );
      }        
      double sum_r2( mu - sigma_u ); 
      double theta_i_times_mu = mu/((double) ndep);
      for( unsigned int i=0; i<_nvar; i++ ){
        if( lst[i].empty() ) continue;   
        if(sum_r2 != 0.){
          const double theta_i = _r2[i] / sum_r2;//(_c2[i] - _U2[i])/(C2 - sigma_u);
          theta_i_times_mu = theta_i*mu;
        }
        const double rowOffsetOver =  - _U1[i] + theta_i_times_mu;  // -Ui + theta_i * mu
        shadow[i].oveEst = relu( shadow[i].oveEst + rowOffsetOver );      // it should be note that we leave to relu(PWL) the determination of directly passing, 
        shadow[i].flag_nonEmpty();
      }
    }       
    return ;
  }

}




template <typename T>
inline
void
ASModel<T>::_add_aggregate_shadow
( std::vector<UnivarPWL<T>>& Alst, const std::vector<UnivarPWL<T>>& Blst,
  std::vector<UnivarPWL<T>>& Ashadow, const std::vector<UnivarPWL<T>>& Bshadow, 
  unsigned & Andep,std::vector<double>& Ashadow_info, const std::vector<double>& Bshadow_info)
const
{
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "    start agrt" << std::endl;
  
  // for( unsigned int i=0; i<_nvar; i++ ){
  //   if( Alst[i].empty() ) continue;
  //   if(Ashadow_info[1]>0) std::cout << "AshadowO " << i << ": "<< Ashadow[i].oveEst << std::endl;
  //   std::cout << "AlstO " << i << ": "<< Alst[i].oveEst << std::endl;   
  //   if(Ashadow_info[0]>0) std::cout << "AshadowU " << i << ": "<< Ashadow[i].undEst << std::endl;
  //   std::cout << "AlstU " << i << ": "<< Alst[i].undEst << std::endl;         
  // }

  // for( unsigned int i=0; i<_nvar; i++ ){
  //   if( Blst[i].empty() ) continue;
  //   if(Bshadow_info[1]>0) std::cout << "BshadowO " << i << ": "<< Bshadow[i].oveEst << std::endl;
  //   std::cout << "BlstO " << i << ": "<< Blst[i].oveEst << std::endl;   
  //   if(Bshadow_info[0]>0) std::cout << "BshadowU " << i << ": "<< Bshadow[i].undEst << std::endl;
  //   std::cout << "BlstU " << i << ": "<< Blst[i].undEst << std::endl;         
  // }
#endif
  // Step 1: addition of active estimators
  unsigned ndep = Andep;
  std::vector<UnivarPWL<T>> ABlst(_nvar);
  bool cmplmtyFlag = true;
  for( unsigned int i=0; i<_nvar; i++ ){
    if( !Alst[i].empty() && !Blst[i].empty() ){
      // std::cout << "    SHADOW A+B " << i << std::endl;
      // std::cout << "      SHADOW A ove " << Alst[i].oveEst << std::endl; 
      // std::cout << "      SHADOW B ove " << Blst[i].oveEst << std::endl; 
      // std::cout << "      SHADOW A und " << Alst[i].undEst << std::endl; 
      // std::cout << "      SHADOW b und " << Blst[i].undEst << std::endl; 
      ABlst[i] = Alst[i] + Blst[i];
      cmplmtyFlag = false;
    }

    else if( !Alst[i].empty() ){
      ABlst[i] = Alst[i];
    }
    else if( !Blst[i].empty() ){
      ndep++;                // add dependency
      ABlst[i] = Blst[i];    // copy entire row
    }
    
  }
#ifdef MC__ASM_DEBUG  
  for( unsigned int i=0; i<_nvar; i++ ){
    if( ABlst[i].empty()) continue;
    if(ABlst[i].oveEst.get_lb()> 1e7){
      std::cout << ABlst[i].oveEst.get_lb() << std::endl;
      std::cout << ABlst[i].oveEst << std::endl;
    }  
    if(ABlst[i].oveEst.get_ub()> 1e7){
      std::cout << ABlst[i].oveEst.get_ub() << std::endl;
      std::cout << ABlst[i].oveEst << std::endl;
    } 
    if(ABlst[i].undEst.get_lb()> 1e7){
      std::cout << ABlst[i].undEst.get_lb() << std::endl;
      std::cout << ABlst[i].undEst << std::endl;    
    }  
    if(ABlst[i].undEst.get_ub()> 1e7){
      std::cout << ABlst[i].undEst.get_ub() << std::endl;
      std::cout << ABlst[i].undEst << std::endl;    
    }       
  }
    for( unsigned int i=0; i<_nvar; i++ ){
    if( ABlst[i].empty()) continue;
    assert(ABlst[i].oveEst.get_lb()< 1e7);
    assert(ABlst[i].oveEst.get_ub()< 1e7);    
    assert(ABlst[i].undEst.get_lb()< 1e7);    
    assert(ABlst[i].undEst.get_ub()< 1e7);  
  }
#endif

#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "    ACT added" << std::endl;
  // std::cout << "        AundSHA: " << Ashadow_info[0] << std::endl;
  // std::cout << "        BundSHA: " << Ashadow_info[1] << std::endl;
  // std::cout << "        AoveSHA: " << Bshadow_info[0] << std::endl;
  // std::cout << "        BoveSHA: " << Bshadow_info[1] << std::endl;
#endif

  // Step 2: prepare for cross addition
  const unsigned int AundSHA = static_cast<unsigned int>(Ashadow_info[0]);
  const unsigned int BundSHA = static_cast<unsigned int>(Bshadow_info[0]);
  const unsigned int AoveSHA = static_cast<unsigned int>(Ashadow_info[1]);
  const unsigned int BoveSHA = static_cast<unsigned int>(Bshadow_info[1]);
  //  // process the aggreagation mode
  // if (){     
  //   std::cout << "    error in aggregating in preprocessing addition" << std::endl; 
  //   throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );   
  // }
  bool undEst2BUpdated(AundSHA + BundSHA); 
  bool oveEst2BUpdated(AoveSHA + BoveSHA);

#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "        AundSHA: " << AundSHA << std::endl;
  std::cout << "        BundSHA: " << BundSHA << std::endl;
  std::cout << "        AoveSHA: " << AoveSHA << std::endl;
  std::cout << "        BoveSHA: " << BoveSHA << std::endl;      
#endif

#ifdef MC__ASM_DEBUG_LOGGING
  MC__ASM_DEBUG_LOGGER.open("MC__ASM_DEBUG_LOGGER.txt", std::ios::app); // append
  if (!MC__ASM_DEBUG_LOGGER.is_open()) throw std::runtime_error("Could not open file: 'MC__ASM_DEBUG_LOGGER.txt'" );
//  if (!MC__ASM_DEBUG_LOGGER.is_open()) std::cerr << "Open log file failed! in shadow aggr addition " << std::endl;
  else{
    MC__ASM_DEBUG_LOGGER << "    SHADOW AGGR" << std::endl;
    MC__ASM_DEBUG_LOGGER << "        AundSHA: " << AundSHA << std::endl;
    MC__ASM_DEBUG_LOGGER << "        BundSHA: " << BundSHA << std::endl;
    MC__ASM_DEBUG_LOGGER << "        AoveSHA: " << AoveSHA << std::endl;
    MC__ASM_DEBUG_LOGGER << "        BoveSHA: " << BoveSHA << std::endl;   
    MC__ASM_DEBUG_LOGGER << "        A:" << std::endl;
    for( unsigned int i=0; i<_nvar; i++ ){
      if( Alst[i].empty()) continue;
      MC__ASM_DEBUG_LOGGER << "            " << Alst[i].oveEst;
      MC__ASM_DEBUG_LOGGER << "            " << Alst[i].undEst; 
      if(AoveSHA)  MC__ASM_DEBUG_LOGGER << "            shadow ove " << Ashadow[i].oveEst;  
      if(AundSHA)  MC__ASM_DEBUG_LOGGER << "            shadow und " << Ashadow[i].undEst;
    }
    MC__ASM_DEBUG_LOGGER << "        B:" << std::endl;
    for( unsigned int i=0; i<_nvar; i++ ){
      if( Blst[i].empty()) continue;
      MC__ASM_DEBUG_LOGGER << "            " << Blst[i].oveEst;
      MC__ASM_DEBUG_LOGGER << "            " << Blst[i].undEst;    
      if(BoveSHA)  MC__ASM_DEBUG_LOGGER << "            shadow ove " << Bshadow[i].oveEst;  
      if(BundSHA)  MC__ASM_DEBUG_LOGGER << "            shadow und " << Bshadow[i].undEst;          
    }
  }  
  MC__ASM_DEBUG_LOGGER.close();

#ifdef MC__ASM_DEBUG_EVAL_LOGGING  
  MC__ASM_EVAL_LOGGER_N.open("MC__ASM_EVAL_LOGGER_N.txt", std::ios::app); // append
  if (!MC__ASM_EVAL_LOGGER_N.is_open()) throw std::runtime_error("Could not open file: 'MC__ASM_EVAL_LOGGER_N.txt'" );
//  if (!MC__ASM_EVAL_LOGGER_N.is_open()) std::cerr << "Open log file failed! in shadow aggr addition " << std::endl;
  else{
    // MC__ASM_EVAL_LOGGER_N << std::setw(24) << debugging_ctr
    //                       << std::setw(24) << AoveSHA 
    //                       << std::setw(24) << AundSHA
    //                       << std::setw(24) << BoveSHA 
    //                       << std::setw(24) << BundSHA << std::endl;      
    MC__ASM_EVAL_LOGGER_N << debugging_ctr << ","
                          << AoveSHA << ","
                          << AundSHA << ","
                          << BoveSHA << ","
                          << BundSHA << "," << std::endl;                                                                           
    for( unsigned int i=0; i<_nvar; i++ ){
      if( Alst[i].empty()) continue;
      Alst[i].oveEst.write2log(MC__ASM_EVAL_LOGGER_N);
      Alst[i].undEst.write2log(MC__ASM_EVAL_LOGGER_N);
      if(AoveSHA)  Ashadow[i].oveEst.write2log(MC__ASM_EVAL_LOGGER_N);
      if(AundSHA)  Ashadow[i].undEst.write2log(MC__ASM_EVAL_LOGGER_N);
    }
    for( unsigned int i=0; i<_nvar; i++ ){
      if( Blst[i].empty()) continue;       
      Blst[i].oveEst.write2log(MC__ASM_EVAL_LOGGER_N);
      Blst[i].undEst.write2log(MC__ASM_EVAL_LOGGER_N);
      if(BoveSHA)  Bshadow[i].oveEst.write2log(MC__ASM_EVAL_LOGGER_N);
      if(BundSHA)  Bshadow[i].undEst.write2log(MC__ASM_EVAL_LOGGER_N);
    }
  }  
  MC__ASM_EVAL_LOGGER_N.close();
  debugging_ctr ++;
#endif
#endif

  // Shortcut for early return if cmplmtyFlag == true
  if(cmplmtyFlag){
    ;
  }


  double ubMaxAllUnd = -DBL_MAX; 
  unsigned int indMaxUbAllUnd = -1;  
  double lbMaxAllUnd = -DBL_MAX; 
  unsigned int indMaxLbAllUnd = 0;    
  std::vector<std::vector<UnivarPWLE<double>>> AactBshaUnd(BundSHA);
  std::vector<std::vector<UnivarPWLE<double>>> AshaBactUnd(AundSHA); 
  std::vector<std::vector<std::vector<UnivarPWLE<double>>>> AshaBshaUnd(AundSHA);

#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "    preparison has been done for cross addition" << std::endl;
#endif


  // Step 3: process the underestimators
  if (undEst2BUpdated){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
    std::cout << "        undEst2BUpdated" << std::endl;
#endif

    // Step 3.1: cross addition for undEst
    // Add Alst.undEst to all BShadow.undEst
    //std::vector<std::vector<UnivarPWL<T>>> AactBshaUnd(BundSHA);
    for(unsigned int iB = 0; iB < BundSHA; iB++){
      //std::cout << "            AactBshaUnd" << std::endl;      
      AactBshaUnd[iB].resize(_nvar);      
      for( unsigned int i=0; i<_nvar; i++ ){
        if( !Alst[i].empty() && !Blst[i].empty() )
          AactBshaUnd[iB][i] = Alst[i].undEst + Bshadow[i].undEst;
        else if( !Alst[i].empty() ){
          AactBshaUnd[iB][i] = Alst[i].undEst;
        }          
        else if( !Blst[i].empty() ){
          AactBshaUnd[iB][i] = Bshadow[i].undEst;    // copy entire row
        }
      }
    }
    // Add Blst.undEst to all AShadow.undEst
    //std::vector<std::vector<UnivarPWL<T>>> AshaBactUnd(AundSHA);    
    for(unsigned int iA = 0; iA < AundSHA; iA++){
      //std::cout << "            AshaBactUnd" << std::endl;
      AshaBactUnd[iA].resize(_nvar);      
      for( unsigned int i=0; i<_nvar; i++ ){
        if( !Alst[i].empty() && !Blst[i].empty() )
          AshaBactUnd[iA][i] = Blst[i].undEst + Ashadow[i].undEst;
        else if( !Alst[i].empty() ){
          AshaBactUnd[iA][i] = Ashadow[i].undEst;    // copy entire row
        }
        else if( !Blst[i].empty() ){
          AshaBactUnd[iA][i] = Blst[i].undEst;
        }        
      }      
    }
 

    // Add AShadow.undEst to all BShadow.undEst
    //std::vector<std::vector<std::vector<UnivarPWL<T>>>> AshaBshaUnd(AundSHA);
    for(unsigned int iA = 0; iA < AundSHA; iA++){                
      AshaBshaUnd[iA].resize(BundSHA);  
      //std::cout << "            AshaBshaUnd: " << BundSHA << std::endl;
      for(unsigned int iB = 0; iB < BundSHA; iB++){
        AshaBshaUnd[iA][iB].resize(_nvar);  
        for(unsigned int i = 0; i < _nvar; i++){
          if( !Alst[i].empty() && !Blst[i].empty() )
            AshaBshaUnd[iA][iB][i] = Ashadow[i].undEst + Bshadow[i].undEst;
          else if( !Alst[i].empty() ){
            AshaBshaUnd[iA][iB][i] = Ashadow[i].undEst;    // copy entire row
          }                  
          else if( !Blst[i].empty() ){
            AshaBshaUnd[iA][iB][i] = Bshadow[i].undEst;    // copy entire row
          }
        }
      }        
    }

    // Step 3.2: select the best ones: 
    // new act: the unique greatest minimum or the greatest minimum with the greatest maximum
    // new sha: the greatest maximum, while in fact the best shadow is the one provides the most info supplementing the act 
    std::vector<double> minimaUnd((AundSHA+1)*(BundSHA+1),0.);
    std::vector<double> maximaUnd((AundSHA+1)*(BundSHA+1),0.);        
    unsigned int cnt = 0;
    for(unsigned int i = 0; i < _nvar; i++){
      if(ABlst[i].empty()) continue;
      minimaUnd[0] += ABlst[i].undEst.get_lb();
      maximaUnd[0] += ABlst[i].undEst.get_ub();
      for(unsigned int iA = 0; iA < AundSHA; iA++){
        minimaUnd[iA+1] += AshaBactUnd[iA][i].get_lb();
        maximaUnd[iA+1] += AshaBactUnd[iA][i].get_ub();
      }
      for(unsigned int iB = 0; iB < BundSHA; iB++){
        minimaUnd[iB+AundSHA+1] += AactBshaUnd[iB][i].get_lb();
        maximaUnd[iB+AundSHA+1] += AactBshaUnd[iB][i].get_ub();
      }
      cnt = 0;
      for(unsigned int iA = 0; iA < AundSHA; iA++){
        for(unsigned int iB = 0; iB < BundSHA; iB++){
          minimaUnd[cnt+BundSHA+AundSHA+1] += AshaBshaUnd[iA][iB][i].get_lb();
          maximaUnd[cnt+BundSHA+AundSHA+1] += AshaBshaUnd[iA][iB][i].get_ub();
          cnt ++;
        }
      }          
    }

    lbMaxAllUnd  = minimaUnd[0]; 
    //indMaxLbAllUnd = 0; 
    //std::cout << "         minimaUnd[0] : maximaUnd[0] " << std::scientific << std::setprecision(12) << minimaUnd[0] << " : " << maximaUnd[0] << std::endl;     
    for(unsigned int i = 1; i < minimaUnd.size(); i++){
    // Note that this is a heuristic for choosing the resultant ACT and SHA
      //std::cout << "         minimaUnd[i] : maximaUnd[i] " << std::scientific << std::setprecision(12) << minimaUnd[i] << " : " << maximaUnd[i] << std::endl;
      if (lbMaxAllUnd < minimaUnd[i] - 5e-14){
        lbMaxAllUnd  = minimaUnd[i];
        indMaxLbAllUnd = i;
      }
      //else if(lbMaxAllUnd == minimaUnd[i] && maximaUnd[indMaxLbAllUnd] < maximaUnd[i]){
      else if(std::fabs(lbMaxAllUnd-minimaUnd[i]) < 5e-14 && maximaUnd[indMaxLbAllUnd] + 5e-14 < maximaUnd[i]){ //else if(isequal(lbMaxAllUnd,minimaUnd[i]) && maximaUnd[indMaxLbAllUnd] + 5e-14 < maximaUnd[i]){          
        indMaxLbAllUnd = i;
      }           
    }

    //ubMaxAllUnd  = -DBL_MAX ;  
    for(unsigned int i = 0; i < maximaUnd.size(); i++){
    // Note that this is a heuristic for choosing the resultant ACT and SHA
      //std::cout << "maximaUnd "  << i << "  " << std::scientific << std::setprecision(18) << maximaUnd[i] << std::endl;
      if(i == indMaxLbAllUnd) continue;
      if (ubMaxAllUnd <= maximaUnd[i] - 5e-14){
        ubMaxAllUnd  = maximaUnd[i];
        indMaxUbAllUnd = i;
      }   
    }    
#ifdef MC__ASM_DEBUG  
    assert(indMaxUbAllUnd < maximaUnd.size());    
#endif
  }

#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "    UND has been processed" << std::endl;
#endif


  // Step 4: process the overestimators

  double ubMinAllOve = DBL_MAX; 
  unsigned int indMinUbAllOve = 0;  
  double lbMinAllOve = DBL_MAX; 
  unsigned int indMinLbAllOve = -1;  
  std::vector<std::vector<UnivarPWLE<double>>> AactBshaOve(BoveSHA);
  std::vector<std::vector<UnivarPWLE<double>>> AshaBactOve(AoveSHA);
  std::vector<std::vector<std::vector<UnivarPWLE<double>>>> AshaBshaOve(AoveSHA);
  if (oveEst2BUpdated){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
    std::cout << "        oveEst2BUpdated" << std::endl;
#endif
    
    // Step 4.1: cross addition for oveEst
    // Add Alst.oveEst to all BShadow.oveEst
    //std::vector<std::vector<UnivarPWL<T>>> AactBshaOve(BoveSHA);
    for(unsigned int iB = 0; iB < BoveSHA; iB++){
      AactBshaOve[iB].resize(_nvar);      
      for( unsigned int i=0; i<_nvar; i++ ){
        if( !Alst[i].empty() && !Blst[i].empty() )
          AactBshaOve[iB][i] = Alst[i].oveEst + Bshadow[i].oveEst;
        else if( !Alst[i].empty() ){
          AactBshaOve[iB][i] = Alst[i].oveEst;
        }          
        else if( !Blst[i].empty() ){
          AactBshaOve[iB][i] = Bshadow[i].oveEst;    // copy entire row
        }
      }
    }
    // Add Blst.oveEst to all AShadow.oveEst
    //std::vector<std::vector<UnivarPWL<T>>> AshaBactOve(AoveSHA);    
    for(unsigned int iA = 0; iA < AoveSHA; iA++){
      AshaBactOve[iA].resize(_nvar);      
      for( unsigned int i=0; i<_nvar; i++ ){
        if( !Alst[i].empty() && !Blst[i].empty() )
          AshaBactOve[iA][i] = Blst[i].oveEst + Ashadow[i].oveEst;
        else if( !Alst[i].empty() ){
          AshaBactOve[iA][i] = Ashadow[i].oveEst;    // copy entire row
        }
        else if( !Blst[i].empty() ){
          AshaBactOve[iA][i] = Blst[i].oveEst;
        }        
      }      
    }


    // Add AShadow.oveEst to all BShadow.oveEst
    //std::vector<std::vector<std::vector<UnivarPWL<T>>>> AshaBshaOve(AoveSHA);
    for(unsigned int iA = 0; iA < AoveSHA; iA++){
      AshaBshaOve[iA].resize(BoveSHA);  
      for(unsigned int iB = 0; iB < BoveSHA; iB++){
        AshaBshaOve[iA][iB].resize(_nvar);  
        for(unsigned int i = 0; i < _nvar; i++){
          if( !Alst[i].empty() && !Blst[i].empty() )
            AshaBshaOve[iA][iB][i] = Ashadow[i].oveEst + Bshadow[i].oveEst;
          else if( !Alst[i].empty() ){
            AshaBshaOve[iA][iB][i] = Ashadow[i].oveEst;    // copy entire row
          }                  
          else if( !Blst[i].empty() ){
            AshaBshaOve[iA][iB][i] = Bshadow[i].oveEst;    // copy entire row
          }
        }
      }        
    }

    // Step 4.2: select the best ones 
    // new act: the unique greatest minimum or the greatest minimum with the greatest maximum
    // new sha: the greatest maximum, while in fact the best shadow is the one provides the most info supplementing the act 
    std::vector<double> minimaOve((AoveSHA+1)*(BoveSHA+1),0.);
    std::vector<double> maximaOve((AoveSHA+1)*(BoveSHA+1),0.);        
    unsigned int cnt = 0;
    for(unsigned int i = 0; i < _nvar; i++){
      if(ABlst[i].empty()) continue;
      minimaOve[0] += ABlst[i].oveEst.get_lb();
      maximaOve[0] += ABlst[i].oveEst.get_ub();
      for(unsigned int iA = 0; iA < AoveSHA; iA++){
        minimaOve[iA+1] += AshaBactOve[iA][i].get_lb();
        maximaOve[iA+1] += AshaBactOve[iA][i].get_ub();
      }
      for(unsigned int iB = 0; iB < BoveSHA; iB++){
        minimaOve[iB+AoveSHA+1] += AactBshaOve[iB][i].get_lb();
        maximaOve[iB+AoveSHA+1] += AactBshaOve[iB][i].get_ub();
      }
      cnt = 0;
      for(unsigned int iA = 0; iA < AoveSHA; iA++){
        for(unsigned int iB = 0; iB < BoveSHA; iB++){
          minimaOve[cnt+BoveSHA+AoveSHA+1] += AshaBshaOve[iA][iB][i].get_lb();
          maximaOve[cnt+BoveSHA+AoveSHA+1] += AshaBshaOve[iA][iB][i].get_ub();
          cnt ++;
        }
      }          
    }
    
    ubMinAllOve  = maximaOve[0]; 
    //indMinUbAllOve = 0;       
    //std::cout << "maximaOve "  << 0 << "  " << std::scientific << std::setprecision(12) << maximaOve[0] << std::endl;        
    for(unsigned int i = 1; i < maximaOve.size(); i++){
    // Note that this is a heuristic for choosing the resultant ACT and SHA
 //     std::cout << "maximaOve "  << i << "  " << std::scientific << std::setprecision(12) << maximaOve[i] << " ubMinAllOve - maximaOve[i] " << ubMinAllOve - maximaOve[i] << std::endl;      
 //     std::cout << isequal(ubMinAllOve,maximaOve[i]) << std::endl;  
      if (ubMinAllOve > maximaOve[i] + 5e-14){
        ubMinAllOve = maximaOve[i];
        indMinUbAllOve = i;
//        std::cout << "indMinUbAllOve set at "  << i << "  " << std::scientific << std::setprecision(12) << maximaOve[i] << " ubMinAllOve - maximaOve[i] " << ubMinAllOve - maximaOve[i] << std::endl;      
      }   
      else if(std::fabs(ubMinAllOve-maximaOve[i]) < 5e-14 && minimaOve[indMinUbAllOve] - 5e-14 > minimaOve[i]){ // else if(isequal(ubMinAllOve,maximaOve[i]) && minimaOve[indMinUbAllOve] - 5e-14 > minimaOve[i]){    //      //else if(ubMinAllOve == maximaOve[i] && minimaOve[indMinUbAllOve] > minimaOve[i]){
        indMinUbAllOve = i;
//        std::cout << "minimaOve "  << i << "  " << std::scientific << std::setprecision(12) << minimaOve[i] << " minimaOve[indMinUbAllOve]- " << minimaOve[indMinUbAllOve] - minimaOve[i] << std::endl;    
      }     
    }
  
    //lbMinAllOve  = DBL_MAX; 
    for(unsigned int i = 0; i < minimaOve.size(); i++){
    // Note that this is a heuristic for choosing the resultant ACT and SHA
      //std::cout << "minimaOve "  << i << "  " << std::scientific << std::setprecision(12) << minimaOve[i] << std::endl;
      if(i == indMinUbAllOve) continue;
      //std::cout << "        lbMinAllOve - minimaOve[i] "   << std::scientific << std::setprecision(12) << lbMinAllOve - minimaOve[i] << std::endl; 
      if (lbMinAllOve >= minimaOve[i] + 5e-14){
        lbMinAllOve = minimaOve[i];
        indMinLbAllOve = i;
        //std::cout << "        set indMinLbAllOve "  << i << std::endl;
      }
    }
    //std::cout << "       minimaOve 1 "  << std::scientific << std::setprecision(12) << minimaOve[1] << std::endl;
    //std::cout << "       minimaOve 3 "  << std::scientific << std::setprecision(12) << minimaOve[1] << std::endl;  
    //std::cout << "    indMinLbAllOve " << indMinLbAllOve << std::endl;
#ifdef MC__ASM_DEBUG      
    assert(indMinLbAllOve < minimaOve.size());   
#endif
  }

#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "    OVE has been processed" << std::endl;
  std::cout << "    indMinLbAllOve : indMinUbAllOve " << indMinLbAllOve << " : " << indMinUbAllOve << std::endl;
  std::cout << "       lbMinAllOve : ubMinAllOve "    << lbMinAllOve    << " : " << ubMinAllOve << std::endl;
  std::cout << "    indMaxLbAllUnd : indMaxUbAllUnd " << indMaxLbAllUnd << " : " << indMaxUbAllUnd << std::endl;  
  std::cout << "       lbMaxAllUnd : ubMaxAllUnd "    << lbMaxAllUnd    << " : " << ubMaxAllUnd << std::endl;
#endif

#ifdef MC__ASM_DEBUG_LOGGING
  MC__ASM_DEBUG_LOGGER.open("MC__ASM_DEBUG_LOGGER.txt", std::ios::app); // append 模式
  if (!MC__ASM_DEBUG_LOGGER.is_open()){
    std::cerr << "Open log file failed! in shadow aggr addition " << std::endl;
  }
  MC__ASM_DEBUG_LOGGER << "    SHADOW AGGR" << std::endl;
  MC__ASM_DEBUG_LOGGER << "    indMinLbAllOve : indMinUbAllOve " << indMinLbAllOve << " : " << indMinUbAllOve << std::endl;
  MC__ASM_DEBUG_LOGGER << "       lbMinAllOve : ubMinAllOve "    << lbMinAllOve    << " : " << ubMinAllOve << std::endl;
  MC__ASM_DEBUG_LOGGER << "    indMaxLbAllUnd : indMaxUbAllUnd " << indMaxLbAllUnd << " : " << indMaxUbAllUnd << std::endl;  
  MC__ASM_DEBUG_LOGGER << "       lbMaxAllUnd : ubMaxAllUnd "    << lbMaxAllUnd    << " : " << ubMaxAllUnd << std::endl;
  MC__ASM_DEBUG_LOGGER.close();
#endif
  // Step 5: assembly the output
  Andep = ndep;
  if(!undEst2BUpdated && !oveEst2BUpdated){

#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "      assembly the output: no update" << std::endl;
#endif

    Alst.swap(ABlst);
    // in this case, both the ove and und do not have shadow
    // hence we do not need to update shadow_info
    return;
  }
  else if(!undEst2BUpdated){    
    Alst.swap(ABlst);
    //Ashadow_info[1] = 0;
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "      assembly the output: Update OVE" << std::endl;
#endif
    Ashadow_info[1] = 1;
    // first: rewrite the swap, second, as we steal the estmators
    // for( unsigned int i=0; i<_nvar; i++ ){
    //   if(ABlst[i].empty()) continue;
    //   std::swap(Alst[i].undEst,ABlst[i].undEst);
    //   Alst[i].flag_nonEmpty();
    // }

    // process the indices indMinUbAllOve and indMinLbAllOve to get to pointer
    // Note that we have known that AoveSHA == 0 and BoveSHA == 0 do not hold together
    // Also note that ABlst has been swapped to Alst
    std::vector<UnivarPWLE<double>> sha2Bset(0);
    std::vector<UnivarPWLE<double>> act2Bset(0);  
    bool setAct = false;  
    bool shouldSwapActSha = false;     
    if(AoveSHA > 0 && BoveSHA == 0){
      // in this case AshaBactOve[0] is chosen when indMinUbAllOve or indMinLbAllOve == 1
      std::swap(sha2Bset,AshaBactOve[0]);
      if (indMinLbAllOve == 0){
        shouldSwapActSha = true;
      }
    }
    else if(AoveSHA == 0){
      // in this case AactBshaOve[0] is chosen when indMinUbAllOve or indMinLbAllOve == 1
      std::swap(sha2Bset,AactBshaOve[0]);
      if (indMinLbAllOve == 0){
        shouldSwapActSha = true;
      }
    }
    else{
      if (indMinLbAllOve == 0){
        switch( indMinUbAllOve ){
          case 1: std::swap(sha2Bset,AshaBactOve[0]); break;
          case 2: std::swap(sha2Bset,AactBshaOve[0]); break;
          case 3: std::swap(sha2Bset,AshaBshaOve[0][0]); break;
        }                
        shouldSwapActSha = true;
      }
      else{
        switch( indMinLbAllOve ){
          case 1: std::swap(sha2Bset,AshaBactOve[0]); break;
          case 2: std::swap(sha2Bset,AactBshaOve[0]); break;
          case 3: std::swap(sha2Bset,AshaBshaOve[0][0]); break;
        }           
        
        switch( indMinUbAllOve ){
          case 0: break;
          case 1: std::swap(act2Bset,AshaBactOve[0]); setAct = true; break;
          case 2: std::swap(act2Bset,AactBshaOve[0]); setAct = true; break;
          case 3: std::swap(act2Bset,AshaBshaOve[0][0]); setAct = true; break;
        }   
      }
    } 

    for( unsigned int i=0; i<_nvar; i++ ){
      if(Alst[i].empty()) continue;
      std::swap(Ashadow[i].oveEst,sha2Bset[i]);
      if(shouldSwapActSha) std::swap(Alst[i].oveEst,Ashadow[i].oveEst);
      else if(setAct) std::swap(Alst[i].oveEst,act2Bset[i]);
      Ashadow[i].flag_nonEmpty();
    }

    // for( unsigned int i=0; i<_nvar; i++ ){
    //   if(Alst[i].empty()) continue;
    //   // if act + act then sha+sha, or act + sha with sha + act

    //   // if(indMinUbAllOve == 0){         
    //   //   std::swap(Ashadow[i].oveEst,AshaBshaOve[0][0][i].oveEst);
    //   // } 
    //   // if(indMinUbAllOve == 0){         
    //   //   std::swap(Alst[i].oveEst,ABlst[i].oveEst);
    //   // }      

    //   if(AoveSHA == 1 && indMinUbAllOve == 1){
    //     std::swap(Alst[i].oveEst,AshaBactOve[0][i].oveEst);
    //   }
    //   else if(indMinUbAllOve == 1 || indMinUbAllOve == 2){
    //     std::swap(Alst[i].oveEst,AactBshaOve[0][i].oveEst);
    //   }
    //   else if(indMinUbAllOve == 3){ // AoveSHA ==1 && BoveSHA == 1 && 
    //     std::swap(Alst[i].oveEst,AshaBshaOve[0][0][i].oveEst);
    //   }            

    //   if(indMinLbAllOve == 0){         
    //     std::swap(Ashadow[i].oveEst,ABlst[i].oveEst);
    //   } 
    //   else if(AoveSHA == 1 && indMinLbAllOve == 1){
    //     std::swap(Ashadow[i].oveEst,AshaBactOve[0][i].oveEst);
    //   }
    //   else if(indMinLbAllOve == 1 || indMinLbAllOve == 2){
    //     std::swap(Ashadow[i].oveEst,AactBshaOve[0][i].oveEst);
    //   }
    //   else if(indMinLbAllOve == 3){
    //     std::swap(Ashadow[i].oveEst,AshaBshaOve[0][0][i].oveEst);
    //   }                  


    // }  

  }
  else if(!oveEst2BUpdated){ 
    Alst.swap(ABlst);
    //Ashadow_info[0] = 0;
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "      assembly the output: step5 Update UND" << std::endl;
#endif
    Ashadow_info[0] = 1;

    // process the indices indMinUbAllUnd and indMinLbAllUnd to get to pointer
    // Note that we have known that AundSHA == 0 and BundSHA == 0 do not hold together
    // Also note that ABlst has been swapped to Alst
    std::vector<UnivarPWLE<double>> sha2Bset(0);
    std::vector<UnivarPWLE<double>> act2Bset(0);      
    bool shouldSwapActSha = false;
    bool setAct = false;     
    if(AundSHA > 0 && BundSHA == 0){
      // in this case AshaBactUnd[0] is chosen when indMaxLbAllUnd or indMaxUbAllUnd == 1
      std::swap(sha2Bset,AshaBactUnd[0]);
      if (indMaxUbAllUnd == 0){
        shouldSwapActSha = true;
      }
    }
    else if(AundSHA == 0){
      // in this case AactBshaUnd[0] is chosen when indMaxLbAllUnd or indMaxUbAllUnd == 1
      std::swap(sha2Bset,AactBshaUnd[0]);
      if (indMaxUbAllUnd == 0){
        shouldSwapActSha = true;
      }
    }
    else{     
      if (indMaxUbAllUnd == 0){
        switch( indMaxLbAllUnd ){
          case 1: std::swap(sha2Bset,AshaBactUnd[0]); break;
          case 2: std::swap(sha2Bset,AactBshaUnd[0]); break;
          case 3: std::swap(sha2Bset,AshaBshaUnd[0][0]); break;
        }                
        shouldSwapActSha = true;
      }
      else{
        switch( indMaxUbAllUnd ){
          case 1: std::swap(sha2Bset,AshaBactUnd[0]); break;
          case 2: std::swap(sha2Bset,AactBshaUnd[0]); break;
          case 3: std::swap(sha2Bset,AshaBshaUnd[0][0]); break;
        }          
        switch( indMaxLbAllUnd ){
          case 0: break;
          case 1: std::swap(act2Bset,AshaBactUnd[0]); setAct = true; break;
          case 2: std::swap(act2Bset,AactBshaUnd[0]); setAct = true; break;
          case 3: std::swap(act2Bset,AshaBshaUnd[0][0]); setAct = true; break;
        }            
      }
    } 

    for( unsigned int i=0; i<_nvar; i++ ){
      if(Alst[i].empty()) continue;
      std::swap(Ashadow[i].undEst,sha2Bset[i]);
      if(shouldSwapActSha) std::swap(Alst[i].undEst,Ashadow[i].undEst);
      else if(setAct) std::swap(Alst[i].undEst,act2Bset[i]);
      Ashadow[i].flag_nonEmpty();     
    }


    // // first: rewrite the swap, second, as we steal the estmators
    // // for( unsigned int i=0; i<_nvar; i++ ){
    // //   if(ABlst[i].empty()) continue;
    // //   std::swap(Alst[i].oveEst,ABlst[i].oveEst);
    // //   Alst[i].flag_nonEmpty();      
    // // }

    // for( unsigned int i=0; i<_nvar; i++ ){
    //   if(Alst[i].empty()) continue;
    //   // if act + act then sha+sha, or act + sha with sha + act

    //   // if(indMaxLbAllUnd == 0){         
    //   //   std::swap(Ashadow[i].undEst,AshaBshaUnd[0][0][i].undEst);
    //   // } 
    //   // if(indMaxLbAllUnd == 0){         
    //   //   std::swap(Alst[i].undEst,ABlst[i].undEst);
    //   // } 

             
    //   if(AundSHA == 1 && indMaxLbAllUnd == 1){
    //     std::swap(Alst[i].undEst,AshaBactUnd[0][i].undEst);
    //   }
    //   else if(indMaxLbAllUnd == 1 || indMaxLbAllUnd == 2){
    //     std::swap(Alst[i].undEst,AactBshaUnd[0][i].undEst);
    //   }
    //   else if(indMaxLbAllUnd == 3){ // AundSHA ==1 && BundSHA == 1 && 
    //     std::swap(Alst[i].undEst,AshaBshaUnd[0][0][i].undEst);
    //   }            
    //   std::cout << "indMaxLbAllUnd = " << indMaxLbAllUnd << std::endl;

    // }  

    //   if(indMaxUbAllUnd == 0){         
    //     std::swap(Ashadow[i].undEst,ABlst[i].undEst);
    //   } 
    //   else if(AundSHA == 1 && indMaxUbAllUnd == 1){
    //     std::swap(Ashadow[i].undEst,AshaBactUnd[0][i].undEst);
    //   }
    //   else if(indMaxUbAllUnd == 1 || indMaxUbAllUnd == 2){
    //     std::swap(Ashadow[i].undEst,AactBshaUnd[0][i].undEst);
    //   }
    //   else if(indMaxUbAllUnd == 3){
    //     std::swap(Ashadow[i].undEst,AshaBshaUnd[0][0][i].undEst);
    //   }
    //   std::cout << "indMaxUbAllUnd = " << indMaxUbAllUnd << std::endl;     


  }
  else{
    Alst.swap(ABlst);
    //Ashadow_info[0] = 0;
    //Ashadow_info[1] = 0;
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "      assembly the output: step5 UPDATE ALL" << std::endl;
#endif
    Ashadow_info[0] = 1;
    Ashadow_info[1] = 1;



  //   for( unsigned int i=0; i<_nvar; i++ ){
  //     if(ABlst[i].empty()) continue;
  //     Alst[i].flag_nonEmpty();
  //     // if act + act then sha+sha, or act + sha with sha + act

  //     if(indMinUbAllOve == 0){         
  //       std::swap(Alst[i].oveEst,ABlst[i].oveEst);
  //     } 
  //     if(AoveSHA == 1 && indMinUbAllOve == 1){
  //       std::swap(Alst[i].oveEst,AshaBactOve[0][i].oveEst);
  //     }
  //     else if(indMinUbAllOve == 1 || indMinUbAllOve == 2){
  //       std::swap(Alst[i].oveEst,AactBshaOve[0][i].oveEst);
  //     }
  //     else if(indMinUbAllOve == 3){ // AoveSHA ==1 && BoveSHA == 1 && 
  //       std::swap(Alst[i].oveEst,AshaBshaOve[0][0][i].oveEst);
  //     }            
  // //std::cout << "          assembly the output: over ACT finished" << std::endl;
  //     if(indMinLbAllOve == 0){         
  //       std::swap(Ashadow[i].oveEst,ABlst[i].oveEst);
  //     } 
  //     else if(AoveSHA == 1 && indMinLbAllOve == 1){
  //       std::swap(Ashadow[i].oveEst,AshaBactOve[0][i].oveEst);
  //     }
  //     else if(indMinLbAllOve == 1 || indMinLbAllOve == 2){
  //       std::swap(Ashadow[i].oveEst,AactBshaOve[0][i].oveEst);
  //     }
  //     else if(indMinLbAllOve == 3){
  //       std::swap(Ashadow[i].oveEst,AshaBshaOve[0][0][i].oveEst);
  //     }                  
  //   }  

    // process the indices indMinUbAllOve and indMinLbAllOve to get to pointer
    // Note that we have known that AoveSHA == 0 and BoveSHA == 0 do not hold together
    // Also note that ABlst has been swapped to Alst
    std::vector<UnivarPWLE<double>> sha2Bset(0);
    std::vector<UnivarPWLE<double>> act2Bset(0);   
    bool setAct = false;  
    bool shouldSwapActSha = false;     
    if(AoveSHA > 0 && BoveSHA == 0){
      // in this case AshaBactOve[0] is chosen when indMinUbAllOve or indMinLbAllOve == 1
      std::swap(sha2Bset,AshaBactOve[0]);
      if (indMinLbAllOve == 0){
        shouldSwapActSha = true;
      }
    }
    else if(AoveSHA == 0){
      // i n this case AactBshaOve[0] is chosen when indMinUbAllOve or indMinLbAllOve == 1
      std::swap(sha2Bset,AactBshaOve[0]);
      if (indMinLbAllOve == 0){
        shouldSwapActSha = true;
      }
    }
    else{
      if (indMinLbAllOve == 0){
//        std::cout << "          indMinUbAllOve " << indMinUbAllOve << std::endl;      
        switch( indMinUbAllOve ){
          case 1: std::swap(sha2Bset,AshaBactOve[0]); break;
          case 2: std::swap(sha2Bset,AactBshaOve[0]); break;
          case 3: std::swap(sha2Bset,AshaBshaOve[0][0]); break;
        }                
        shouldSwapActSha = true;
      }
      else{
//        std::cout << "          indMinLbAllOve " << indMinLbAllOve << std::endl;
        switch( indMinLbAllOve ){
          case 1: std::swap(sha2Bset,AshaBactOve[0]); break;
          case 2: std::swap(sha2Bset,AactBshaOve[0]); break;
          case 3: std::swap(sha2Bset,AshaBshaOve[0][0]); break;
        }           
        
        switch( indMinUbAllOve ){
          case 0: break;
          case 1: std::swap(act2Bset,AshaBactOve[0]); setAct = true; break;
          case 2: std::swap(act2Bset,AactBshaOve[0]); setAct = true; break;
          case 3: std::swap(act2Bset,AshaBshaOve[0][0]); setAct = true; break;
        }   
      }
    } 

#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "          before swap" << std::endl;
#endif
//  std::cout << "          sha2Bset.size() " << sha2Bset.size() << std::endl;
    for( unsigned int i=0; i<_nvar; i++ ){
      if(Alst[i].empty()) continue;
//      std::cout << "          before swap 1" << std::endl;
      std::swap(Ashadow[i].oveEst,sha2Bset[i]);
//      std::cout << "          before swap 2" << std::endl;
      if(shouldSwapActSha){std::swap(Alst[i].oveEst,Ashadow[i].oveEst);}
      else if(setAct){ std::swap(Alst[i].oveEst,act2Bset[i]);}
//      std::cout << "          before swap 7" << std::endl;
      Ashadow[i].flag_nonEmpty();
//      std::cout << "          before swap 8" << std::endl;
    }


#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "          assembly the output: over finished" << std::endl;
#endif

    // process the indices indMinUbAllUnd and indMinLbAllUnd to get to pointer
    // Note that we have known that AundSHA == 0 and BundSHA == 0 do not hold together
    // Also note that ABlst has been swapped to Alst

    std::vector<UnivarPWLE<double>> sha2BsetUnd(0);
    std::vector<UnivarPWLE<double>> act2BsetUnd(0);  
    bool shouldSwapActShaUnd = false;
    bool setActUnd = false;     
    if(AundSHA > 0 && BundSHA == 0){
      // in this case AshaBactUnd[0] is chosen when indMaxLbAllUnd or indMaxUbAllUnd == 1
      std::swap(sha2BsetUnd,AshaBactUnd[0]);
      if (indMaxUbAllUnd == 0){
        shouldSwapActShaUnd = true;
      }
    }
    else if(AundSHA == 0){
      // in this case AactBshaUnd[0] is chosen when indMaxLbAllUnd or indMaxUbAllUnd == 1
      std::swap(sha2BsetUnd,AactBshaUnd[0]);
      if (indMaxUbAllUnd == 0){
        shouldSwapActShaUnd = true;
      }
    }
    else{     
      if (indMaxUbAllUnd == 0){
  // std::cout << "          indMaxLbAllUnd " << indMaxLbAllUnd << std::endl;      
        switch( indMaxLbAllUnd ){
          case 1: std::swap(sha2BsetUnd,AshaBactUnd[0]); break;
          case 2: std::swap(sha2BsetUnd,AactBshaUnd[0]); break;
          case 3: std::swap(sha2BsetUnd,AshaBshaUnd[0][0]); break;
        }                
        shouldSwapActShaUnd = true;
      }
      else{
  // std::cout << "          indMaxUbAllUnd " << indMaxUbAllUnd << std::endl;
        switch( indMaxUbAllUnd ){
          case 1: std::swap(sha2BsetUnd,AshaBactUnd[0]); break;
          case 2: std::swap(sha2BsetUnd,AactBshaUnd[0]); break;
          case 3: std::swap(sha2BsetUnd,AshaBshaUnd[0][0]); break;
        }          
        switch( indMaxLbAllUnd ){
          case 0: break;
          case 1: std::swap(act2BsetUnd,AshaBactUnd[0]); setActUnd = true; break;
          case 2: std::swap(act2BsetUnd,AactBshaUnd[0]); setActUnd = true; break;
          case 3: std::swap(act2BsetUnd,AshaBshaUnd[0][0]); setActUnd = true; break;
        }            
      }
    } 

#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "          before swap" << std::endl;
#endif
  //std::cout << "          sha2BsetUnd.size() " << sha2BsetUnd.size() << std::endl;

    for( unsigned int i=0; i<_nvar; i++ ){
      if(Alst[i].empty()) continue;
      std::swap(Ashadow[i].undEst,sha2BsetUnd[i]);
      if(shouldSwapActShaUnd) std::swap(Alst[i].undEst,Ashadow[i].undEst);
      else if(setActUnd) std::swap(Alst[i].undEst,act2BsetUnd[i]);      
    }


  //   for( unsigned int i=0; i<_nvar; i++ ){
  //     if(ABlst[i].empty()) continue;
  //     // if act + act then sha+sha, or act + sha with sha + act

  //     if(indMaxLbAllUnd == 0){         
  //       std::swap(Alst[i].undEst,ABlst[i].undEst);
  //     } 
  //     else if(AundSHA == 1 && indMaxLbAllUnd == 1){
  //       std::swap(Alst[i].undEst,AshaBactUnd[0][i].undEst);
  //     }
  //     else if(indMaxLbAllUnd == 1 || indMaxLbAllUnd == 2){
  //       std::swap(Alst[i].undEst,AactBshaUnd[0][i].undEst);
  //     }
  //     else if(indMaxLbAllUnd == 3){ // AundSHA ==1 && BundSHA == 1 && 
  //       std::swap(Alst[i].undEst,AshaBshaUnd[0][0][i].undEst);
  //     }            
  // //std::cout << "          assembly the output: under ACT finished" << std::endl;
  //     if(indMaxUbAllUnd == 0){         
  //       std::swap(Ashadow[i].undEst,ABlst[i].undEst);
  //     } 
  //     else if(AundSHA == 1 && indMaxUbAllUnd == 1){
  //       std::swap(Ashadow[i].undEst,AshaBactUnd[0][i].undEst);
  //     }
  //     else if(indMaxUbAllUnd == 1 || indMaxUbAllUnd == 2){
  //       std::swap(Ashadow[i].undEst,AactBshaUnd[0][i].undEst);
  //     }
  //     else if(indMaxUbAllUnd == 3){
  //       std::swap(Ashadow[i].undEst,AshaBshaUnd[0][0][i].undEst);
  //     }                  
  //   }  

#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "          assembly the output: under finished" << std::endl;
#endif

  }
  // std::cout << "i " << 0 <<" " << Ashadow[0].empty() << std::endl; 
  // std::cout << "i " << 1 <<" " << Ashadow[1].empty() << std::endl;
  
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
  std::cout << "    ADD AGGR FINISHED" << std::endl;
#endif

}


//! @brief C++ class for interval superposition models of factorable functions: variable
////////////////////////////////////////////////////////////////////////
//! mc::ASVar is a C++ class for propagation of interval superposition
//! models (ASM) through factorable functions. The template parameter
//! corresponds to the type used to propagate the interval coefficients.
////////////////////////////////////////////////////////////////////////
template <typename T> 
class ASVar
{
  template <typename U> friend class ASModel;

  template <typename U> friend std::ostream& operator<<
    ( std::ostream &, ASVar<U> const& );

  template <typename U> friend ASVar<U> operator+
    ( ASVar<U> const& );
  template <typename U> friend ASVar<U> operator+
    ( ASVar<U> const&, ASVar<U> const& );
  template <typename U> friend ASVar<U> operator+
    ( double const&, ASVar<U> const& );
  template <typename U> friend ASVar<U> operator+
    ( ASVar<U> const&, double const& );


  template <typename U> friend ASVar<U> operator+
    ( ASVar<U> && );
  template <typename U> friend ASVar<U> operator+
    ( ASVar<U> const&, ASVar<U> && );
  template <typename U> friend ASVar<U> operator+
    ( ASVar<U> &&, ASVar<U> const& );    
  template <typename U> friend ASVar<U> operator+
    ( ASVar<U> &&, ASVar<U> && );    

  template <typename U> friend ASVar<U> operator+
    ( double const&, ASVar<U> && );
  template <typename U> friend ASVar<U> operator+
    ( ASVar<U> &&, double const& );


  template <typename U> friend ASVar<U> operator-
    ( ASVar<U> const& );
  template <typename U> friend ASVar<U> operator-
    ( ASVar<U> const&, ASVar<U> const& );
  template <typename U> friend ASVar<U> operator-
    ( double const&, ASVar<U> const& );
  template <typename U> friend ASVar<U> operator-
    ( ASVar<U> const&, double const& );


  template <typename U> friend ASVar<U> operator-
    ( ASVar<U> && );

  template <typename U> friend ASVar<U> operator-
    ( ASVar<U> const&, ASVar<U> && );
  template <typename U> friend ASVar<U> operator-
    ( ASVar<U> &&, ASVar<U> const& );
  template <typename U> friend ASVar<U> operator-
    ( ASVar<U> &&, ASVar<U> && );

  template <typename U> friend ASVar<U> operator-
    ( double const&, ASVar<U> && );
  template <typename U> friend ASVar<U> operator-
    ( ASVar<U> &&, double const& );



  template <typename U> friend ASVar<U> operator*
    ( ASVar<U> const&, ASVar<U> const& );
  template <typename U> friend ASVar<U> operator*
    ( double const&, ASVar<U> const& );
  template <typename U> friend ASVar<U> operator*
    ( ASVar<U> const&, double const& );
  template <typename U> friend ASVar<U> operator*
    ( double const&, ASVar<U> && );
  template <typename U> friend ASVar<U> operator*
    ( ASVar<U> &&, double const& );


  template <typename U> friend ASVar<U> operator/
    ( ASVar<U> const&, ASVar<U> const& );
  template <typename U> friend ASVar<U> operator/
    ( double const&, ASVar<U> const& );
  template <typename U> friend ASVar<U> operator/
    ( ASVar<U> const&, double const& );

  template <typename U> friend ASVar<U> max
    ( ASVar<U> const&, ASVar<U> const& );
  template <typename U> friend ASVar<U> max
    ( ASVar<U> const&, double const& );
  template <typename U> friend ASVar<U> max
    ( ASVar<U> &&, double const& );
  // template <typename U> friend ASVar<U> min
  //   ( ASVar<U> const&, ASVar<U> const& );
  // template <typename U> friend ASVar<U> min
  //   ( ASVar<U> const&, double const& );
  // template <typename U> friend ASVar<U> min
  //   ( ASVar<U> &&, double const& );

  // template <typename U> friend ASVar<U> inv
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> inv
  //   ( ASVar<U>&& );
  // template <typename U> friend ASVar<U> sqr
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> sqr
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> sqrt
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> sqrt
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> fabs
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> fabs
  //   ( ASVar<U> && );  
  template <typename U> friend ASVar<U> relu
    ( ASVar<U> const& );
  template <typename U> friend ASVar<U> relu
    ( ASVar<U> && );    
  // template <typename U> friend ASVar<U> exp
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> exp
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> log
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> log
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> xlog
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> xlog
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> sin
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> sin
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> cos
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> cos
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> tan
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> tan
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> tanh
  //   ( ASVar<U> const& );
  // template <typename U> friend ASVar<U> tanh
  //   ( ASVar<U> && );
  // template <typename U> friend ASVar<U> pow
  //   ( ASVar<U> const&, int const& n );
  // template <typename U> friend ASVar<U> pow
  //   ( ASVar<U> &&, int const& n );
  // template <typename U> friend ASVar<U> pow
  //   ( ASVar<U> const&, double const& );
  // template <typename U> friend ASVar<U> pow
  //   ( ASVar<U> &&, double const&);    
  // template <typename U> friend ASVar<U> intersect
  //   ( ASVar<U> const& , U);
  // template <typename U> friend ASVar<U> intersect
  //   ( ASVar<U> && , U);  
  // template <typename U> friend ASVar<U> affine_transform
  //   ( std::vector<ASVar<U> >const& , std::vector<double> const&, const double);
  // template <typename U> friend ASVar<U> affine_transform
  //   ( std::vector<ASVar<U> >&& , std::vector<double> const&, const double);
  // template <typename U> friend ASVar<U> affine_transform
  //   ( std::vector<std::reference_wrapper<ASVar<U>> > , std::vector<double> const&, const double);        
  // template <typename U> friend ASVar<U> cheb
  //   ( ASVar<U> const&, unsigned int const& n );
  // template <typename U> friend ASVar<U> cheb
  //   ( ASVar<U> &&, unsigned int const& n );

  // template <typename U> friend ASVar<U> mid 
  // (const ASVar<U>&, const ASVar<U>&, const double);

 private:
 
  //! @brief Pointer to interval superposition model
  mutable ASModel<T>* _mod;
  //! @brief Number of partitions
  unsigned long long _ndiv;                   // To support up to 2^64 partitions
  //! @brief Number of variables
  unsigned int _nvar;
  //! @brief Number of dependencies
  unsigned int _ndep;
  //! @brief estimator list
  std::vector<UnivarPWL<T>> _lst;

  //! @brief Constant value (in case _mod is nullptr)
  double _cst;
  //! @brief Constant value (in case _lst is empty)
  std::pair<std::vector<double>,bool> _lnr;
 
  //! @brief Variable bound
  mutable std::pair<T,bool> _bnd;

  //! @brief shadow estimator list
  std::vector<UnivarPWL<T>> _shadow;  
  //! @brief shadow estimator information
  mutable std::vector<double> _shadow_info;
  //! @brief upper bound/constant overestimator
  double _oveCut;
  //! @brief lower bound/constant underestimator
  double _undCut;

  // The hierarchical priority is as follows:
  // if(!mod), then only _cst is valid; 
  // else if(_lnr.second), then only _cst and _lnr are valid; 
  // else if(!_mod->options.SHADOW_USE), then only _lst and _bnd are valid
  

 public:
 
  ASVar<T>& operator+=
    ( ASVar<T> const& );
  ASVar<T>& operator+=
    ( double const& );
 
  ASVar<T>& operator-=
    ( ASVar<T> const& );
  ASVar<T>& operator-=
    ( double const& );
 
  ASVar<T>& operator*=
    ( ASVar<T> const& );
  ASVar<T>& operator*=
    ( double const& );
 
  ASVar<T>& operator/=
    ( ASVar<T> const& );
  ASVar<T>& operator/=
    ( double const& );

  ASVar<T>& operator=
  ( double const& cst )
  {

  // The hierarchical priority:
  // if(!mod), then only _cst is valid; 
  // else if(_lnr.second), then only _cst and _lnr are valid; 
  // else if(!_mod->options.SHADOW_USE), then only _lst and _bnd are valid

    _mod = nullptr;
    _cst = cst;
 
    // Privious implementation
    /*******************************/    
    // _mod = nullptr;
    // _nvar = 0;
    // _ndiv = 0;
    // _ndep = 0;
    // _lst.clear(); // may not resize _lst
    // _cst = cst;
    // _bnd = std::make_pair( 0., false );
    // if (_mod->options.SHADOW_USE){
    //     _shadow.clear(); // added for allowing shadow enhancement
    //     _shadow.resize(0);    
    //     _shadow_info.clear();    
    //     _shadow_info.resize(3,0);
    //     _oveCut = cst;
    //     _undCut = cst;
    // }        

    return *this;
  }

  ASVar<T>& operator=
  ( ASVar<T> const& var )
  {
#ifdef MC__ASMODEL_TRACE
    std::cerr << "-- ASVar<T>& operator= ( ASVar<T> const& )\n";
#endif

  // The hierarchical priority:
  // if(!mod), then only _cst is valid; 
  // else if(_lnr.second), then only _cst and _lnr are valid; 
  // else if(!_mod->options.SHADOW_USE), then only _lst and _bnd are valid

    if( this == &var )
      return *this;

    _mod = var._mod;
    if( !_mod ){
      _cst = var._cst;
    }
    else if(var._lnr.second){
      _cst = var._cst;
      _lnr = var._lnr;
    }
    else{
      _lnr = std::make_pair(std::vector<double>(0),false);
      _nvar = var._nvar;
      _ndiv = var._ndiv;
      _ndep = var._ndep;
      _lst = var._lst;
      _bnd = var._bnd;     

      if(_mod->options.SHADOW_USE){
        _shadow = var._shadow;
        _shadow_info = var._shadow_info;
        _oveCut = var._oveCut;
        _undCut = var._undCut;  
      } 

    }
    


    // Privious implementation
    /*******************************/    
    // if( this == &var )
    //   return *this;
    // _mod = var._mod;
    // if( !_mod ) _cst = var._cst;
    // _nvar = var._nvar;
    // _ndiv = var._ndiv;
    // _ndep = var._ndep;
#ifdef TEST_MOVE    
    std::cout << "  Copy Assignment" << std::endl;
#endif
    // _lst = var._lst;
    // _bnd = var._bnd;
    // _shadow = var._shadow;
    // _shadow_info = var._shadow_info;
    // _oveCut = var._oveCut;
    // _undCut = var._undCut;    
    return *this;
  }

  ASVar<T>& operator=
  ( ASVar<T> && var )
  {
#ifdef MC__ASMODEL_TRACE
    std::cerr << "-- ASVar<T>& operator= ( ASVar<T> && )\n";
#endif

  // The hierarchical priority:
  // if(!mod), then only _cst is valid; 
  // else if(_lnr.second), then only _cst and _lnr are valid; 
  // else if(!_mod->options.SHADOW_USE), then only _lst and _bnd are valid

    if( this == &var )
      return *this;

    _mod = var._mod;
    if( !_mod ){
      _cst = var._cst;
    }
    else if(var._lnr.second){
      _cst = var._cst;
      _lnr = std::move(var._lnr);
    }
    else{
      _lnr = std::make_pair(std::vector<double>(0),false);
      _nvar = var._nvar;
      _ndiv = var._ndiv;
      _ndep = var._ndep;
      _lst = std::move(var._lst);
      _bnd = std::move(var._bnd);     

      if(_mod->options.SHADOW_USE){
        _shadow_info = std::move(var._shadow_info);        
        _shadow = std::move(var._shadow);
        _oveCut = var._oveCut;
        _undCut = var._undCut;  
      } 
      
    }
    


    // Privious implementation
    /*******************************/  
    // if( this == &var )
    //   return *this;
    // _mod = var._mod;
    // if( !_mod ) _cst = var._cst;

    // _nvar = var._nvar;
    // _ndiv = var._ndiv;
    // _ndep = var._ndep;
#ifdef TEST_MOVE
    std::cout << "  Move Assignment" << std::endl;
#endif
    // _lst = std::move(var._lst);
    // _bnd = std::move(var._bnd);
    // _shadow = std::move(var._shadow);  
    // _shadow_info = std::move(var._shadow_info);
    // _oveCut = var._oveCut;
    // _undCut = var._undCut;      
    return *this;
  }

  ASVar
  ( ASModel<T>* const mod )
  : _mod(mod), _ndiv(mod->_ndiv), _nvar(mod->_nvar), _ndep(0), _lst(_nvar),_cst(0.),_lnr(0,false),_bnd(0.,false),
    _shadow(0),_shadow_info(3,0),_oveCut(DBL_MAX),_undCut(-DBL_MAX)
  {
 
    if (_mod->options.SHADOW_USE){
        _shadow.resize(_nvar);    
    }    
  }

  ASVar
  ( ASModel<T>* const mod, unsigned int ndx, T const& bnd )
  : _mod(mod), _cst(0),_lnr(std::vector<double>(mod->_nvar,0.),true)
  {
    if( ndx >= mod->_nvar )
      throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::INDEX );
   
  //_cst(0.5*(Op<T>::u(bnd) + Op<T>::l(bnd)))

    _lnr.first[ndx] = 1.0;
    //_lst[ndx].debug_check_overNunder_flags();
    _mod->_defvar[ndx] = true;
    _mod->_bndvar[ndx] = bnd;
  
  } 
  


  // ASVar
  // ( ASModel<T>* const mod, unsigned int ndx, T const& bnd )
  // : _mod(mod), _ndiv(mod->_ndiv), _nvar(mod->_nvar), _ndep(1), _lst(_nvar), _bnd(bnd,true),
  //   _shadow(0),_shadow_info(3,0),_oveCut(Op<T>::u(bnd)),_undCut(Op<T>::l(bnd))
  // {
  //   if( ndx >= _nvar )
  //     throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::INDEX );
   
  //   _lst[ndx] = UnivarPWL(bnd);
  //   //_lst[ndx].debug_check_overNunder_flags();
  //   _mod->_defvar[ndx] = true;
  //   _mod->_bndvar[ndx] = bnd;


  //   if (_mod->options.SHADOW_USE){
  //       _shadow.resize(_nvar);    
  //   }    
  
  // } 


  ASVar
  ( double const& cst=0. )
  : _mod(nullptr), _cst(cst)
  {}
 
  // ASVar
  // ( double const& cst=0. )
  // : _mod(nullptr), _ndiv(0), _nvar(0), _ndep(0), _cst(cst), _lnr(0,false), _bnd(0.,false),
  //   _shadow(_nvar),_shadow_info(3,0),_oveCut(cst),_undCut(cst)
  // {}

  ASVar
  ( ASVar<T> const& var )
  : _mod(var._mod)
  {
#ifdef MC__ASMODEL_TRACE
    std::cerr << "-- ASVar( ASVar<T> const& )\n";
#endif

  // The hierarchical priority:
  // if(!mod), then only _cst is valid; 
  // else if(_lnr.second), then only _cst and _lnr are valid; 
  // else if(!_mod->options.SHADOW_USE), then only _lst and _bnd are valid

    if( this == &var )
      return;
    else if( !_mod ){
      _cst = var._cst;
    }
    else if(var._lnr.second){
      _cst = var._cst;
      _lnr = var._lnr;
    }
    else{
      _lnr = std::make_pair(std::vector<double>(0),false);
      _nvar = var._nvar;
      _ndiv = var._ndiv;
      _ndep = var._ndep;
      _lst = var._lst;
      _bnd = var._bnd;     

      if(_mod->options.SHADOW_USE){
        _shadow_info = var._shadow_info;        
        _shadow = var._shadow;
        _oveCut = var._oveCut;
        _undCut = var._undCut;  
      } 
      
    }
    
  }


//   ASVar
//   ( ASVar<T> const& var )
//   : _mod(var._mod), _ndiv(var._ndiv), _nvar(var._nvar), _ndep(var._ndep),
//     _oveCut(var._oveCut),_undCut(var._undCut)
//   {
// #ifdef MC__ASMODEL_TRACE
//     std::cerr << "-- ASVar( ASVar<T> const& )\n";
// #endif
//     if( this == &var ) return;
//     if( !_mod ) _cst = var._cst;
// #ifdef TEST_MOVE
//     std::cout << "  Copy Constructor" << std::endl;
// #endif
//     _lst = var._lst;
//     _bnd = var._bnd;
//     _shadow = var._shadow;
//     _shadow_info = var._shadow_info;
//   }


  ASVar
  ( ASVar<T> && var )
  : _mod(var._mod)
  {
  // The hierarchical priority:
  // if(!mod), then only _cst is valid; 
  // else if(_lnr.second), then only _cst and _lnr are valid; 
  // else if(!_mod->options.SHADOW_USE), then only _lst and _bnd are valid    
#ifdef MC__ASMODEL_TRACE
    std::cerr << "-- ASVar( ASVar<T> && var )\n";
#endif
    if( this == &var ) return;
#ifdef TEST_MOVE    
    std::cout << "  Move Constructor" << std::endl;
#endif
    else if( !_mod ){
      _cst = var._cst;
    }
    else if(var._lnr.second){
      _cst = var._cst;
      _lnr = std::move(var._lnr);
    }
    else{
      _lnr = std::make_pair(std::vector<double>(0),false);
      _nvar = var._nvar;
      _ndiv = var._ndiv;
      _ndep = var._ndep;
      _lst = std::move(var._lst);
      _bnd = std::move(var._bnd);     

      if(_mod->options.SHADOW_USE){
        _shadow_info = std::move(var._shadow_info);        
        _shadow = std::move(var._shadow);
        _oveCut = var._oveCut;
        _undCut = var._undCut;  
      } 
     
    }
     
  }





//   ASVar
//   ( ASVar<T> && var )
//   : _mod(var._mod), _ndiv(var._ndiv), _nvar(var._nvar), _ndep(var._ndep),
//     _oveCut(var._oveCut),_undCut(var._undCut)
//   {
// #ifdef MC__ASMODEL_TRACE
//     std::cerr << "-- ASVar( ASVar<T> && var )\n";
// #endif
//     if( this == &var ) return;
//     if( !_mod ) _cst = var._cst;
// #ifdef TEST_MOVE    
//     std::cout << "  Move Constructor" << std::endl;
// #endif

//     _lst = std::move(var._lst);
//     _bnd = std::move(var._bnd);
//     _shadow = std::move(var._shadow); 
//     _shadow_info = std::move(var._shadow_info);       
//   }



  ASVar
  ( ASVar<T> const& var, const double mtpr ) // multiplier
  : _mod(var._mod), _ndiv(var._ndiv), _nvar(var._nvar), _ndep(var._ndep), _lst(var._nvar),_lnr(0,false), _bnd(0.,true),
    _shadow(_nvar),_shadow_info(3,0),_oveCut(DBL_MAX),_undCut(-DBL_MAX)
  {
#ifdef MC__ASMODEL_TRACE
    std::cerr << "-- ASVar( ASVar<T> const&, const double )\n";
#endif
#ifdef TEST_MOVE
    std::cout << "Copy and Scale Constructor" << std::endl;
#endif

    if( this == &var ) return;
    if( !_mod ) {_cst = var._cst*mtpr; return;}
    else if(var._lnr.second){
      std::cout << "ASVar constructor for mtpr*var has not yet been defined " << std::endl;
      throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::INTERN );
    }

    // if( !_ndep )
    //   throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::INTERN );


    for( unsigned int i=0; i<_nvar; i++ ){
      if( var._lst[i].empty() ) continue;
      //_lst[i].reserve(_ndiv);
      _lst[i] = UnivarPWL(var._lst[i],mtpr);
    }
    
    if( var._bnd.second ){
      _bnd.first = var._bnd.first * mtpr;
    } 
    else{
      _bnd.second = false;
    }

    //_shadow = var._shadow;
    //_shadow_slope = var._shadow_slope;
    //_shadow_info = var._shadow_info;


    if (_mod->options.SHADOW_USE){
      std::cout << "shadow has not yet been implemented" << std::endl;
      throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::INTERN );
      //_shadow_info[0] = 1.;
    }    


  if (_mod->options.SHADOW_USE && _ndep > 1){
    std::cout << "shadow has not yet been implemented" << std::endl;
    throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::INTERN );
    }

  }



  ~ASVar
  () 
  {
#ifdef FASM_LIFITIME_DEBUG 
    std::cout<< "ASV delated, nvar = " <<_nvar <<std::endl;
#endif
  }


  ASModel<T>* mod
  ()
  const
  {
    return _mod;
  }    

  ASVar<T>& set
  ( ASModel<T>* const mod )
  {
    _mod = mod;
    _nvar = mod->_nvar;
    _ndiv = mod->_ndiv;
    _ndep = 0;
    _lst.clear();
    _lst.resize( _nvar );
    _lnr = std::make_pair(std::vector<double>(0), false );
    _bnd = std::make_pair( 0., false );
    if (_mod->options.SHADOW_USE){
      _shadow.clear();
      _shadow.resize(_nvar);
      _shadow_info.resize(3,0);
      _oveCut =  DBL_MAX;
      _undCut = -DBL_MAX;
    }    

    return *this;
  }
  

  // ASVar<T>& set
  // ( ASModel<T>* const mod )
  // {
  //   _mod = mod;
  //   _nvar = mod->_nvar;
  //   _ndiv = mod->_ndiv;
  //   _ndep = 0;
  //   _lst.clear();
  //   _lst.resize( _nvar );
  //   _bnd = std::make_pair( 0., false );
  //   if (_mod->options.SHADOW_USE){
  //     _shadow.clear();
  //     _shadow.resize(_nvar);
  //     _oveCut =  DBL_MAX;
  //     _undCut = -DBL_MAX;
  //   }    

  //   return *this;
  // }




  ASVar<T>& set
  ( ASModel<T>* const mod, unsigned int ndx, T const& bnd )
  {
    _mod = mod;
    _nvar = mod->_nvar;
    if( ndx >= _nvar )
      throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::INDEX );

    _cst = 0.; // 0.5*(Op<T>::u(bnd) + Op<T>::l(bnd))
    _lnr = std::make_pair(std::vector<double>(_nvar,0.), true );
    _lnr.first[ndx] = 1.0;
    _mod->_defvar[ndx] = true;
    _mod->_bndvar[ndx] = bnd;

    // _ndiv = mod->_ndiv;
    // _ndep = 1;
    // _lst.clear();
    // _lst.resize( _nvar );

    // _bnd = std::make_pair( bnd, true );

    // if (_mod->options.SHADOW_USE){
    //   _shadow_info.resize(3,0);
    //   _shadow.clear();
    //   _shadow.resize(_nvar);
    //   _oveCut = Op<T>::u(bnd);
    //   _undCut = Op<T>::l(bnd);     
    // }
    return *this;
  }


  // ASVar<T>& set
  // ( ASModel<T>* const mod, unsigned int ndx, T const& bnd )
  // {
  //   _mod = mod;
  //   _nvar = mod->_nvar;
  //   if( ndx >= _nvar )
  //     throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::INDEX );

  //   _ndiv = mod->_ndiv;
  //   _ndep = 1;
  //   _lst.clear();
  //   _lst.resize( _nvar );
  //   _lst[ndx] = UnivarPWL(bnd);
  //   _mod->_defvar[ndx] = true;
  //   _mod->_bndvar[ndx] = bnd;
  //   _bnd = std::make_pair( bnd, true );

  //   if (_mod->options.SHADOW_USE){
  //     _shadow.clear();
  //     _shadow.resize(_nvar);
  //     _oveCut = Op<T>::u(bnd);
  //     _undCut = Op<T>::l(bnd);     
  //   }
  //   return *this;
  // }

  unsigned int get_ASVar
  ()
  const
  { 
  // The hierarchical priority:
  // if(!mod), then only _cst is valid; 
  // else if(_lnr.second), then only _cst and _lnr are valid; 
  // else if(!_mod->options.SHADOW_USE), then only _lst and _bnd are valid    
    if(!_mod) return 1;
    else if(_lnr.second) return 2;
    else{
      if(!_mod->options.SHADOW_USE)
        return 3;
      else if(_shadow_info[0] + _shadow_info[1] == 0 )
        return 3;
      else
        return 4;
    }
  }


  double const& get_cst
  ()
  const
  { return _cst; }

  std::pair<std::vector<double>,bool> const & get_lnr_debug
  ()
  const  
  { return _lnr; }


  std::vector<double> const & get_lnr
  ()
  const  
  { return _lnr.first; }


  std::vector<UnivarPWL<T>> const& get_lst
  ()
  const
  { return _lst; }


  std::vector<UnivarPWL<T>> const& get_shadow
  ()
  const  
  { return _shadow; }
  

  std::vector<double> const& get_shadow_info
  ()
  const
  { return _shadow_info; }




  std::vector<std::vector<double>> const& COve()   
  const
  { 
    return _mod->_COve; 
  }

  std::vector<std::vector<double>> const& CUnd()
  const
  { 

    return _mod->_CUnd; 
  }

  void get_C() 
  const
  { 
    if(!_mod) 
      _mod->_compute_C(_cst);
    else if(_lnr.second){
      _mod->_compute_C(_lnr.first,_cst);
    }
    else
      _mod->_compute_C(_lst,_ndep);
    return ;
  }



  std::vector<std::vector<T>> const& C() 
  const
  { 
//    _mod->_compute_C(_lst,_ndep);
    return _mod->_COut; 
  }


  double const& cst
  ()
  const
  { return _cst; }

  unsigned int const& ndep
  ()
  const
  { return _ndep; }

  unsigned int& ndep
  ()
  { return _ndep; }

  double ub
  ()
  const
  { return Op<T>::u(bound()); }

  double lb
  ()
  const
  { return Op<T>::l(bound()); }

  T B
  ()
  const
  { return bound(); }

  T bound
  ()
  const
  {
  // The hierarchical priority:
  // if(!mod), then only _cst is valid; 
  // else if(_lnr.second), then only _cst and _lnr are valid; 
  // else if(!_mod->options.SHADOW_USE), then only _lst and _bnd are valid    
    //std::cout << "in bound " << _cst << std::endl;
    if(!_mod) return T(_cst);
    else if(_lnr.second){
      T bnd(_cst);
      //std::cout << "cst " << _cst << std::endl;
      for(unsigned int i = 0; i <(_mod->_nvar); i++){
        if( _lnr.first[i] == 0. ) continue;
        bnd += (_mod->_bndvar[i]) * _lnr.first[i];
        //std::cout << "coef "<< i << " " <<_lnr.first[i] << std::endl;
      }
      _bnd = std::make_pair(bnd,true);
      return bnd;
    } 
    else{
      if( _bnd.second ) return _bnd.first;
      _bnd.first = ( _mod? _mod->_B( _lst ): _cst );
      _bnd.second = true;
      return _bnd.first;
    }

  }

   T eval
  ( double const* const point )
  const  
  {
    assert( point );
    if( !_mod ) return _cst;
    else if(_lnr.second){
      //std::cout << "This ASVar is linear on its domain" << std::endl;
      double val = _cst; 
      for( unsigned int i=0; i<(_mod->_nvar); i++ ){
        if( _lnr.first[i] == 0. ) continue;
        val += point[i]*_lnr.first[i]; 
      }
      return T(val);
    } 
    else{    
      std::pair<double,double> val( 0.,0. );
      std::pair<double,double> valSha( 0.,0. );    
      for( unsigned int i=0; i<_nvar; i++ ){
        if( _lst[i].empty() ) continue;
  
        double lb = _lst[i].undEst.eval(point[i]); 
        double ub = _lst[i].oveEst.eval(point[i]);
        if(_mod->options.SHADOW_USE){
          if(_shadow_info[0] > 0)
            valSha.first  += _shadow[i].undEst.eval(point[i]); 
          if(_shadow_info[1] > 0)
            valSha.second += _shadow[i].oveEst.eval(point[i]); 
        }
      // double l = _lst[i].undEst.lbVar();
      // double u = _lst[i].oveEst.ubVar();
      // if(point[i] < l || point[i] > u){
      //   std::cout<<"l: "<<l<<" p: "<<point[i]<<" u: "<<u<<std::endl;
      //   const double _eps (1e-14);
      //   if (point[i] - u >= _eps || point[i] - l <= -_eps)
      //     assert( point[i] >= l && point[i] <= u );
      // }
        val.first += lb;
        val.second += ub;      
      }  
      if(!(_mod->options.SHADOW_USE)){
          _shadow_info.resize(3,0);
      }
      double lb = _shadow_info[0] > 0? std::max(val.first, valSha.first): val.first;
      double ub = _shadow_info[1] > 0? std::min(val.second,valSha.second):val.second;
      if (ub - lb > 0.)      
        return T(lb,ub);    
      else if (lb - ub < 1e2*MC__ASM_COMPUTATION_TOL)
        return T(ub,lb);
      else{ 
        std::cout <<std::setprecision(18)<< "Eval ERROR! ub - lb = " << ub - lb <<  " at ";
        for( unsigned int i=0; i<_nvar; i++ ){
          if( _lst[i].empty() ) continue;
          std::cout << " i: " << point[i];
        }
        std::cout << std::endl;  
        throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );
      }
    }
  }

  std::ostream& display
  ( bool toDisplayBoxes = false, const int& opt=0, std::ostream& out=std::cout )
  const
  {
    if( _mod && (!_lnr.second) ){
      if (toDisplayBoxes){
        //std::vector<std::vector<T>> mat;
        //if(_mod->options.SHADOW_USE)
        //  _mod->_compute_C(_shadow,_ndep);
        //else 
          _mod->_compute_C(_lst,_ndep);
        _mod->_dispvar( _mod->_COut, _ndep, opt, out );
      }
      else{
        _mod->_dispvar( _lst, _ndep, opt, out );
      }
    } 

    return out;
  }
  

  void debug_check_overNunder_flags()
  const
  {
    if(!_mod){ 
      std::cout << "this ASVar is constant" << std::endl;
      return ;
    }
    else if(_lnr.second){
      std::cout << "this ASVar is affine" << std::endl;      
      return ;
    }
    for( unsigned int i=0; i<_nvar; i++ ){
      if( _lst[i].empty() ) continue;
      std::cout << "active " << i <<std::endl;
      _lst[i].debug_check_overNunder_flags();
      if(_mod->options.SHADOW_USE){
        if(_shadow_info[0] > 0 || _shadow_info[1] > 0 )
        std::cout << "shadow " << i <<std::endl;
        _shadow[i].debug_check_overNunder_flags(_shadow_info);
      }    
    }
  }


 private:

  ASVar
  ( ASModel<T> const* const mod )
  : _mod(mod), _ndiv(mod->_ndiv), _nvar(mod->_nvar ), _lst(_nvar), _bnd(0.,false),_shadow(_nvar),_shadow_info(3,0)
  { }

};



template <typename T>  
std::ostream& operator<<
( std::ostream& out, ASVar<T> const& var )
{
  
  // The hierarchical priority:
  // if(!mod), then only _cst is valid; 
  // else if(_lnr.second), then only _cst and _lnr are valid; 
  // else if(!_mod->options.SHADOW_USE), then only _lst and _bnd are valid    
  if(!var._mod){
    out << std::right << std::setw(5) << "This ASVar is constant " << var._cst << std::endl;
  }
  else if(var._lnr.second){
    out << std::right << std::setw(5) << "This ASVar is affine : " << std::endl;
    out << std::right << std::setw(5) << "cst" << std::setw(8);
    for( unsigned int i=0; i<(var._lnr.first.size()); i++ ){
      if( var._lnr.first[i] == 0. ) continue;
      out << "   " << var._mod->bndvar()[i] << "  " << std::setw(8);   
    }
    out << std::endl;
    out << std::right << std::setw(5) << var._cst << std::setw(8);
    for( unsigned int i=0; i<(var._lnr.first.size()); i++ ){
      if( var._lnr.first[i] == 0. ) continue;
      out << " + " << var._lnr.first[i] << " * x_" << i  << std::setw(8);   
    }
    out << std::endl;
  }
  else{
    auto const& lst = var._lst;
    for( unsigned int i=0; !lst.empty() && i<var._nvar; i++ ){
      if( lst[i].empty() ) continue;
      out << std::right << std::setw(5) << "Var No." << i <<": " << std::endl;
      out << lst[i].undEst << std::endl;
      out << lst[i].oveEst;
      // for( unsigned int j=0; j<var._ndiv; j++ ){
      //   if( j && !(j%3) ) out << std::endl << "       ";
      //   out << std::setw(0) << mat[i][j];
      // }
      out << std::endl;
    }

    if(var._mod->options.SHADOW_USE){
      auto const& sha = var._shadow;
      out << std::right << std::setw(5) << "The shadow estimators: (" << var._shadow_info[0] << ",    " << var._shadow_info[1] << ")" << std::endl;
      for( unsigned int i=0; !sha.empty() && i<var._nvar; i++ ){
        if( sha[i].empty() ) continue;
        out << std::right << std::setw(5) << "Var No." << i <<": " << std::endl;
        if( var._shadow_info[0] > 0)
          out << sha[i].undEst << std::endl;
        if( var._shadow_info[1] > 0)
          out << sha[i].oveEst;
        // for( unsigned int j=0; j<var._ndiv; j++ ){
        //   if( j && !(j%3) ) out << std::endl << "       ";
        //   out << std::setw(0) << mat[i][j];
        // }
        out << std::endl;
      }
    }
    //std::cout << "to compute the bound" << std::endl;
    out << std::right << std::setw(7) << "B: " << var.B() << std::endl;    

  }
  

  return out;
}

template <typename T>
inline
ASVar<T> operator+
( ASVar<T> const& var )
{
  return var;
}

template <typename T>
inline
ASVar<T> operator+
( ASVar<T> && varIn )
{
  ASVar<T> var( std::move(varIn) );  
  return var;
}


template <typename T>
inline
ASVar<T>& ASVar<T>::operator+=
( double const& cst )
{
  if( cst == 0. ){
    return *this;
  }
  if( !_mod ){
    _cst += cst;
    return *this;
  }
  else if (_lnr.second){
    _cst += cst;
    return *this;
  }
  else if( !_ndep )
    throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::INTERN );
  else{
    for( unsigned int i=0; i<_nvar; i++ ){
      if( _lst[i].empty() ) continue;
      
      _lst[i] += cst / (double)_ndep;
      if (_mod->options.SHADOW_USE){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW    
        std::cout << "SHADOW ADDITION cst" << std::endl;
#endif
        if (_shadow_info[0] > 0 ){
          _shadow[i].undEst += cst / (double)_ndep;
  
        }
        if (_shadow_info[1] > 0){
          _shadow[i].oveEst += cst / (double)_ndep;
        }
      }      
    }
    if( _bnd.second ) _bnd.first += cst;
  }

  return *this;
}

template <typename T>
inline
ASVar<T>& ASVar<T>::operator+=
( ASVar<T> const& var )
{
  if( !_mod && !var._mod ){
    _cst += var._cst;
    return *this;
  }
  else if( !_mod ){
    const double copy_cst = _cst;
    *this = var;
    *this += copy_cst;
    return *this;
  }
  else if( !var._mod ){
    *this += var._cst;
    return *this;
  }
  else if( _mod != var._mod )
    throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::MODEL );

  else if (_lnr.second && var._lnr.second){
    _cst += var._cst;
    for(unsigned int i = 0; i < (_mod->_nvar); i++)  // _lnr.first.size()  == _mod->_nvar
      _lnr.first[i] += var._lnr.first[i];
    return *this;
  }
  else if(var._lnr.second){
    for(unsigned int i = 0; i < _nvar; i++){
      if(var._lnr.first[i] ==0.) continue; 
      if (_lst[i].empty()){
        _lst[i] = UnivarPWL((_mod->_bndvar[i]),var._lnr.first[i]);
        _ndep ++;
        if(_mod->options.SHADOW_USE){
          if (_shadow_info[0] > 0)
            _shadow[i].undEst = _lst[i].undEst;
          if (_shadow_info[1] > 0)
            _shadow[i].oveEst = _lst[i].oveEst;
        }
      }
      else{
        if(_mod->options.SHADOW_USE){
          UnivarPWL _temp((_mod->_bndvar[i]),var._lnr.first[i]);
          _lst[i] += _temp; 
          if (_shadow_info[0] > 0)
            _shadow[i].undEst = _temp.undEst;
          if (_shadow_info[1] > 0)
            _shadow[i].oveEst = _temp.oveEst;
        }
        else{
          _lst[i] += UnivarPWL((_mod->_bndvar[i]),var._lnr.first[i]);           
        }
      }
    }    
    if(var._cst != 0){
      _cst = var._cst/ (double) _ndep;    
      for(unsigned int i = 0; i < _nvar; i++){
        if(_lst[i].empty()) continue;
        _lst[i] += _cst; 
        if(_mod->options.SHADOW_USE){
           if (_shadow_info[0] > 0)
            _shadow[i].undEst += _cst;
          if (_shadow_info[1] > 0)
            _shadow[i].oveEst += _cst;         
        }
      }    
    }
    _bnd.second = false;
    return *this;
  }
  else if(_lnr.second){
    _ndiv = _mod->_ndiv;
    _nvar = _mod->_nvar;
    _lst.resize(_nvar);
    _ndep = 0;
    _bnd = std::make_pair(T(_cst),true);
    for(unsigned int i = 0; i < _nvar; i++){
      if(_lnr.first[i] ==0.) continue; 
      _bnd.first += (_mod->_bndvar[i])*_lnr.first[i];   
      _lst[i] = UnivarPWL((_mod->_bndvar[i]),_lnr.first[i]);  
      _ndep ++; 
    }
    if(_cst != 0){
      _cst = _cst/ (double) _ndep;    
      for(unsigned int i = 0; i < _nvar; i++){
        if(_lnr.first[i] == 0.) continue; 
        _lst[i] += _cst; 
      }    
    }
    _lnr.second = false;
    _lnr.first.resize(0);
    
    if(_mod->options.SHADOW_USE){
      _shadow_info.resize(3,0.);
      _shadow.resize(_nvar);
      _oveCut = 0;
      _undCut = 0;
    }
  }


  // ASVar
  // ( ASModel<T>* const mod, unsigned int ndx, T const& bnd )
  // : _mod(mod), _ndiv(mod->_ndiv), _nvar(mod->_nvar), _ndep(1), _lst(_nvar), _bnd(bnd,true),
  //   _shadow(0),_shadow_info(3,0),_oveCut(Op<T>::u(bnd)),_undCut(Op<T>::l(bnd))
  // {
  //   if( ndx >= _nvar )
  //     throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::INDEX );
   
  //   _lst[ndx] = UnivarPWL(bnd);
  //   //_lst[ndx].debug_check_overNunder_flags();
  //   _mod->_defvar[ndx] = true;
  //   _mod->_bndvar[ndx] = bnd;


  //   if (_mod->options.SHADOW_USE){
  //       _shadow.resize(_nvar);    
  //   }    
  
  // } 


  if (_mod->options.SHADOW_USE){

#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW    
    std::cout << "SHADOW ADDITION" << std::endl;
    std::cout << "    shadowinfo A UND vs OVE: " << _shadow_info[0] << " " << _shadow_info[1] << std::endl;
    std::cout << "    shadowinfo B UND vs OVE: " << var._shadow_info[0] << " " << var._shadow_info[1] << std::endl;
#endif
    if (_ndep > 1 || var._ndep > 1){
      //std::cout << "SHADOW FLAG CHECK before" << std::endl;
      //debug_check_overNunder_flags(); 
      _mod->_add_aggregate_shadow( _lst, var._lst, _shadow, var._shadow,_ndep, _shadow_info, var._shadow_info);
      //std::cout << "SHADOW FLAG CHECK later" << std::endl;
      //debug_check_overNunder_flags(); 
      _bnd.second = false;
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW    
    std::cout << "SHADOW ADDITION RETURN" << std::endl;
#endif
      return *this; 
    }
  }

  for( unsigned int i=0; i<_nvar; i++ ){
    if( !_lst[i].empty() && !var._lst[i].empty() )
      _lst[i] += var._lst[i];
    else if( !var._lst[i].empty() ){
      _ndep++;                // add dependency
      _lst[i] = var._lst[i];  // copy entire row
    }
  }
  _bnd.second = false;

  return *this;
}

template <typename T>
inline
ASVar<T> operator+
( ASVar<T> const& var1, ASVar<T> const& var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst + var2._cst );
  if( var1._mod ){
    ASVar<T> var3( var1 );
    var3 += var2;
    return var3;
  }
  ASVar<T> var3( var2 );
    var3 += var1;
    return var3;
}


template <typename T>
inline
ASVar<T> operator+
( ASVar<T> const& var1, ASVar<T> && var2 )
{
  //std::cout << "+ 1move" << std::endl;
  if( !var1._mod && !var2._mod ) 
    return( var1._cst + var2._cst );
   
  if( var2._mod ){
    ASVar<T> var3( std::move(var2) );
    var3 += var1;
    return var3;
  }
  ASVar<T> var3( var1 );
  var3 += var2;
  return var3;
}


template <typename T>
inline
ASVar<T> operator+
( ASVar<T> && var1, ASVar<T> const& var2 )
{
  //std::cout << "+ 2move" << std::endl;
  if( !var1._mod && !var2._mod ) 
    return( var1._cst + var2._cst );
   
  if( var1._mod ){
    ASVar<T> var3( std::move(var1) );
    var3 += var2;
    return var3;
  }
  ASVar<T> var3( var2 );
  var3 += var1;
  return var3;
}


template <typename T>
inline
ASVar<T> operator+
( ASVar<T> && var1, ASVar<T> && var2 )
{
  //std::cout << "+ 3move" << std::endl;
  if( !var1._mod && !var2._mod ) 
    return( var1._cst + var2._cst );
  if( var1._mod ){
    ASVar<T> var3( std::move(var1) );
    var3 += var2;
    return var3;
  }
  ASVar<T> var3( std::move(var2) );
    var3 += var1;
    return var3;
}



template <typename T>
inline
ASVar<T> operator+
( ASVar<T> const& var1, double const& cst2 )
{
  if( !var1._mod ) 
    return( var1._cst + cst2 );
  ASVar<T> var3( var1 );
  var3 += cst2;
  return var3;
}



template <typename T>
inline
ASVar<T> operator+
( double const& cst1, ASVar<T> const& var2 )
{
  if( !var2._mod ) 
    return( cst1 + var2._cst );
  ASVar<T> var3( var2 );
  var3 += cst1;
  return var3;
}


template <typename T>
inline
ASVar<T> operator+
( ASVar<T> && var1, double const& cst2 )
{
  if( !var1._mod ) 
    return( var1._cst + cst2 );
  ASVar<T> var3( std::move(var1) );
  var3 += cst2;
  return var3;
}


template <typename T>
inline
ASVar<T> operator+
( double const& cst1, ASVar<T> && var2 )
{
  if( !var2._mod ) 
    return( cst1 + var2._cst );
  ASVar<T> var3( std::move(var2) );
  var3 += cst1;
  return var3;
}



template <typename T> 
inline
ASVar<T> operator-
( ASVar<T> const& var )
{
  ASVar<T> var2( var );
  if( !var2._mod ){
    var2._cst *= -1;
  }
  else if (var2._lnr.second){
    var2._cst = -var2._cst;
    for(unsigned int i = 0; i < (var2._lnr.first.size()); i++)  // _lnr.first.size()  == _mod->_nvar
      var2._lnr.first[i] = -var2._lnr.first[i];
  }
  else if( !var2._ndep )
    throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::INTERN );
  else{ 
    auto&& lst2 = var2._lst;
    for( unsigned int i=0; i<var2._nvar; i++ ){
      if( lst2[i].empty() ) continue;
      lst2[i] *= -1;
    
    
      if (var2._mod->options.SHADOW_USE){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW    
        std::cout << "SHADOW NEG" << std::endl;
#endif
        if (var2._shadow_info[0] > 0 ){
          var2._shadow[i].undEst *= -1;
        }
        if (var2._shadow_info[1] > 0){
          var2._shadow[i].oveEst *= -1;
        }
        std::swap(var2._shadow[i].undEst,var2._shadow[i].oveEst);
      }  
    }
  
    if(var2._mod->options.SHADOW_USE){
      std::swap(var2._shadow_info[0],var2._shadow_info[1]);  
    }
    auto&& bnd = var2._bnd;
    if( bnd.second ) bnd.first *= -1;
  }
  
  return var2;
}





template <typename T> 
inline
ASVar<T> operator-
( ASVar<T> && var )
{
  ASVar<T> var2( std::move(var) );
  if( !var2._mod ){
    var2._cst *= -1;
  }
  else if (var2._lnr.second){
    var2._cst = -var._cst;
    for(unsigned int i = 0; i < (var2._lnr.first.size()); i++)  // _lnr.first.size()  == _mod->_nvar
      var2._lnr.first[i] = -var2._lnr.first[i];
  }
  else if( !var2._ndep )
    throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::INTERN );
  else{ 
    auto&& lst2 = var2._lst;
    for( unsigned int i=0; i<var2._nvar; i++ ){
      if( lst2[i].empty() ) continue;
      lst2[i] *= -1;
    
      if (var2._mod->options.SHADOW_USE){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW    
      std::cout << "SHADOW NEG" << std::endl;
#endif
        if (var2._shadow_info[0] > 0 ){
          var2._shadow[i].undEst *= -1;
        }
        if (var2._shadow_info[1] > 0){
          var2._shadow[i].oveEst *= -1;
        }
        std::swap(var2._shadow[i].undEst,var2._shadow[i].oveEst);
      }  
    }
  
    if(var2._mod->options.SHADOW_USE){
      std::swap(var2._shadow_info[0],var2._shadow_info[1]);   
    }  

  auto&& bnd = var2._bnd;
  if( bnd.second ) bnd.first *= -1;
  }
  return var2;
}






template <typename T>
inline
ASVar<T>& ASVar<T>::operator-=
( double const& cst )
{
  *this += -cst;
  return *this;
}

template <typename T>
inline
ASVar<T>& ASVar<T>::operator-=
( ASVar<T> const& var )
{
  if( !_mod && !var._mod ){
    _cst -= var._cst;
    return *this;
  }
  else if( !_mod ){
    const double copy_cst = _cst;
    *this = -var;
    *this += copy_cst;
    return *this;
  }
  else if( !var._mod ){
    *this -= var._cst;
    return *this;
  }
  else if( _mod != var._mod )
    throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::MODEL );
  else if (_lnr.second && var._lnr.second){
    _cst -= var._cst;
    for(unsigned int i = 0; i < (_mod->_nvar); i++)  // _lnr.first.size()  == _mod->_nvar
      _lnr.first[i] -= var._lnr.first[i];
    return *this;
  }
  else if(var._lnr.second){
    for(unsigned int i = 0; i < _nvar; i++){
      if(var._lnr.first[i] ==0. ) continue;
      if (_lst[i].empty()){
        _lst[i] = UnivarPWL((_mod->_bndvar[i]),(-var._lnr.first[i]));
        _ndep ++;
        if(_mod->options.SHADOW_USE){
          if (_shadow_info[0] > 0)
            _shadow[i].undEst = _lst[i].undEst;
          if (_shadow_info[1] > 0)
            _shadow[i].oveEst = _lst[i].oveEst;
        }
      }
      else{
        if(_mod->options.SHADOW_USE){
          UnivarPWL _temp((_mod->_bndvar[i]),(-var._lnr.first[i]));
          _lst[i] += _temp; 
          if (_shadow_info[0] > 0)
            _shadow[i].undEst = _temp.undEst;
          if (_shadow_info[1] > 0)
            _shadow[i].oveEst = _temp.oveEst;
        }
        else{
          _lst[i] += UnivarPWL((_mod->_bndvar[i]),(-var._lnr.first[i]));           
        }
      }
    }    
    if(var._cst != 0){
      _cst = -var._cst/ (double) _ndep;    
      for(unsigned int i = 0; i < _nvar; i++){
        if(_lst[i].empty()) continue;
        _lst[i] += _cst; 
        if(_mod->options.SHADOW_USE){
           if (_shadow_info[0] > 0)
            _shadow[i].undEst += _cst;
          if (_shadow_info[1] > 0)
            _shadow[i].oveEst += _cst;         
        }
      }    
    }
    _bnd.second = false;
    return *this;
  }
  else if(_lnr.second){
    _ndiv = _mod->_ndiv;
    _nvar = _mod->_nvar;
    _lst.resize(_nvar);
    _ndep = 0;
    _bnd = std::make_pair(T(_cst),true);
    for(unsigned int i = 0; i < _nvar; i++){
       if(_lnr.first[i] ==0.) continue;
      _bnd.first += (_mod->_bndvar[i])*_lnr.first[i];   
      _lst[i] = UnivarPWL((_mod->_bndvar[i]),_lnr.first[i]);
      _ndep ++; 
    }
    if(_cst != 0){
      _cst = _cst/ (double) _ndep;    
      for(unsigned int i = 0; i < _nvar ; i++){
        if(_lnr.first[i] ==0.) continue;
        _lst[i] += _cst; 
      }    
    }
    _lnr.second = false;
    _lnr.first.resize(0);
    
    if(_mod->options.SHADOW_USE){
      _shadow_info.resize(3,0.);
      _shadow.resize(_nvar);
      _oveCut = 0;
      _undCut = 0;
    }
  }
  else if( !_ndep )
    throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::INTERN );



  if(_mod->options.SHADOW_USE){
    std::cout << "operator -= error" << std::endl;
    throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );
  }

  for( unsigned int i=0; i<_nvar; i++ ){
    if( !_lst[i].empty() && !var._lst[i].empty() )
      _lst[i] -= var._lst[i];
    else if( !var._lst[i].empty() ){
      _ndep++;            // add dependency
      _lst[i] = -var._lst[i];
    }
  }
  _bnd.second = false;
  return *this;

}

template <typename T>
inline
ASVar<T> operator-
( ASVar<T> const& var1, ASVar<T> const& var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst - var2._cst );
  if( var1._mod ){
    ASVar<T> var3( var1 );
    var3 += (-var2);
    return var3;
  }
  ASVar<T> var3( -var2 );
    var3 += var1;
    return var3;
}


template <typename T>
inline
ASVar<T> operator-
( ASVar<T> const& var1, ASVar<T> && var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst - var2._cst );
  if( var2._mod ){
    ASVar<T> var3( std::move(-var2) );
    var3 += var1;  
    return var3;
  }
  ASVar<T> var3( var1 );
  var3 += (-var2); 
  return var3;
}


template <typename T>
inline
ASVar<T> operator-
( ASVar<T> && var1, ASVar<T> const& var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst - var2._cst );
  if( var1._mod ){
    ASVar<T> var3( std::move(var1) );
    var3 += (-var2); 
    return var3;
  }
  ASVar<T> var3( -var2 );
  var3 += var1;
  return var3;
}


template <typename T>
inline
ASVar<T> operator-
( ASVar<T> && var1, ASVar<T> && var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst - var2._cst );
  if( var1._mod ){
    ASVar<T> var3( std::move(var1) );
    var3 += (-var2);  
    return var3;
  }
  ASVar<T> var3( std::move(-var2) );
  var3 += var1;
  return var3;
}


template <typename T>
inline
ASVar<T> operator-
( ASVar<T> && var1, double const& cst2 )
{
  if( !var1._mod ) 
    return( var1._cst - cst2 );
  ASVar<T> var3( std::move(var1) );
  var3 += -cst2;
  return var3;
}


template <typename T>
inline
ASVar<T> operator-
( double const& cst1, ASVar<T> && var2 )
{
  if( !var2._mod ) 
    return( cst1 - var2._cst );
  ASVar<T> var3( std::move(-var2) );
  var3 += cst1;
  return var3;
}



template <typename T>
inline
ASVar<T> operator-
( ASVar<T> const& var1, double const& cst2 )
{
  if( !var1._mod ) 
    return( var1._cst - cst2 );
  ASVar<T> var3( var1 );
  var3 += -cst2;
  return var3;
}

template <typename T>
inline
ASVar<T> operator-
( double const& cst1, ASVar<T> const& var2 )
{
  if( !var2._mod ) 
    return( cst1 - var2._cst );
  ASVar<T> var3( -var2 );
  var3 += cst1;
  return var3;
}

template <typename T>
inline
ASVar<T>& ASVar<T>::operator*=
( double const& cst )
{
  // if( cst == 0. ){
  //   *this = 0.;
  //   return *this;
  // }
  if( isequal(cst,0.) ){
    *this = 0.;
    return *this;
  }  
  if( cst == 1. ){
    return *this;
  }
  if( !_mod ){
    _cst *= cst;
    return *this;
  }
  else if (_lnr.second){
    _cst *= cst;
    for(unsigned int i = 0; i < (_mod->_nvar); i++)  // _lnr.first.size()  == _mod->_nvar
      _lnr.first[i] *= cst;
    return *this;
  }
  else if( !_ndep )
    throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::INTERN );
  else{
    for( unsigned int i=0; i<_nvar; i++ ){
      if( _lst[i].empty() ) continue;
//      std::cout << "  before ACT und " << i << _lst[i].undEst << std::endl;
//      std::cout << "  before ACT ove " << i << _lst[i].oveEst << std::endl;        
      _lst[i] *= cst;
      //std::cout << "  after ACT und " << i << _lst[i].undEst << std::endl;
      //std::cout << "  after ACT ove " << i << _lst[i].oveEst << std::endl;      
      if (_mod->options.SHADOW_USE){
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW    
        std::cout << "SHADOW MULTIPLY" << std::endl;
#endif
        if (_shadow_info[0] > 0 ){
          _shadow[i].undEst *= cst;
            //std::cout << "  SHA und " << i << _shadow[i].undEst << std::endl;
        }
        if (_shadow_info[1] > 0){
          _shadow[i].oveEst *= cst;
//std::cout << "  SHA ove " << i << _shadow[i].oveEst << std::endl;
        }
        if(cst < 0.){
          std::swap(_shadow[i].undEst,_shadow[i].oveEst);
        }
      }  
    }
  
    if(_mod->options.SHADOW_USE){
      if(cst < 0.) std::swap(_shadow_info[0],_shadow_info[1]);
    }  
  
    if( _bnd.second ) _bnd.first *= cst;
    return *this;
  }
}

template <typename T>
inline
ASVar<T>& ASVar<T>::operator*=
( ASVar<T> const& var )
{
  if( !_mod && !var._mod ){
    _cst *= var._cst;
    return *this;
  }
  if( !_mod ){
    const double copy_cst = _cst;
    *this = var;
    *this *= copy_cst;
    return *this;
  }
  if( !var._mod ){
    *this *= var._cst;
    return *this;
  }
  if( _mod != var._mod )
    throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::MODEL );

  throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );

  return *this;
}

template <typename T>
inline
ASVar<T> operator*
( ASVar<T> const& var1, ASVar<T> const& var2 )
{
  if( !var1._mod && !var2._mod ) 
    return( var1._cst * var2._cst );
  if( !var2._mod ){
    ASVar<T> var3( var1 );
    var3 *= var2._cst;
    return var3;
  }
  if( !var1._mod ){
    ASVar<T> var3( var2 );
    var3 *= var1._cst;
    return var3;
  }
  throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );
  ASVar<T> var3( var2 );
  var3 *= var1;
  return var3;
  
}

template <typename T>
inline
ASVar<T> operator*
( ASVar<T> const& var1, double const& cst2 )
{
  if( !var1._mod ) 
    return( var1._cst * cst2 );
    ASVar<T> var3( var1 );
    var3 *= cst2;
    return var3;
}

template <typename T>
inline
ASVar<T> operator*
( double const& cst1, ASVar<T> const& var2 )
{
  if( !var2._mod ) 
    return( cst1 * var2._cst );
  ASVar<T> var3( var2 );
  var3 *= cst1;
  return var3;

}


template <typename T>
inline
ASVar<T> operator*
( ASVar<T> && var1, double const& cst2 )
{
  if( !var1._mod ) 
    return( var1._cst * cst2 );
  ASVar<T> var3( std::move(var1) );
  var3 *= cst2;
  return var3;
}


template <typename T>
inline
ASVar<T> operator*
( double const& cst1, ASVar<T> && var2 )
{
  if( !var2._mod ) 
    return( cst1 * var2._cst );
  ASVar<T> var3( std::move(var2) );
  var3 *= cst1;
  return var3;
}

template <typename T>
inline
ASVar<T>& ASVar<T>::operator/=
( double const& cst )
{
  if( cst == 0. )
    throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::DIV );
  *this *= inv( cst );
  return *this;
}

template <typename T>
inline
ASVar<T>& ASVar<T>::operator/=
( ASVar<T> const& var )
{
   throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );
  return *this;
}


template <typename T>
inline
ASVar<T> operator/
( ASVar<T> const& var1, double const& cst2 )
{
  ASVar<T> var3( var1 );
  var3 /= cst2;
  return var3;
}


template <typename T>
inline
ASVar<T> operator/
( ASVar<T> const& var1, ASVar<T> const& var2 )
{
  throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );
  return var1;
}


template <typename T>
inline
ASVar<T> operator/
( double const& cst1, ASVar<T> const& var2 )
{
  throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );
  return var2;
}



template <typename T>
inline
ASVar<T> relu
( ASVar<T> const& var )
{
#ifdef MC__ASM_DEBUG_TRACE  
  std::cout << " SD relu &" << std::endl;	
#endif
  if( !var._mod )
    return std::max(var._cst,0.);
  //std::cout << "  relu" << std::endl;
  if(var._mod->options.SHADOW_USE)
  {
    //std::cout << "in relu" << std::endl;     
    var._mod->_IntmdtCntnrSeted = false;
    T lazyBnd = var.bound();
    if(Op<T>::u(lazyBnd) < MC__ASM_COMPUTATION_TOL){
      ASVar<T> var2( 0. );  // note that the _bnd will be set to (0,false) in the constructor
      //std::cout << "end relu" << std::endl; 
      return var2;
    }  

  // if(Op<T>::l(lazyBnd) > -MC__ASM_COMPUTATION_TOL){
  //   ASVar<T> var2( var );
  //   var2._bnd.second = false;
  //   var2._mod->_IntmdtCntnrSeted = false;
  //   return var2;
  // } 

    ASVar<T> var2( var ); 
    if(var2._lnr.second){
      var2._ndiv = var2._mod->_ndiv;
      var2._nvar = var2._mod->_nvar;
      var2._lst.resize(var2._nvar);
      var2._ndep = 0;
      for(unsigned int i = 0; i < var2._nvar; i++){
        if(var2._lnr.first[i] == 0.) continue;
        var2._lst[i] = UnivarPWL((var2._mod->_bndvar[i]),var2._lnr.first[i]); 
        var2._ndep ++; 
      }
      if(var2._cst != 0){
        double __cst = var2._cst/ (double) var2._ndep;    
        for(unsigned int i = 0; i < var2._nvar; i++){
          if(var2._lnr.first[i] == 0.) continue;
          var2._lst[i] += __cst; 
        }    
      }
      var2._lnr.second = false;
      var2._lnr.first.resize(0);
    
      var2._shadow_info.resize(3,0.);
      var2._shadow.resize(var2._nvar);
      var2._oveCut = DBL_MAX;
      var2._undCut = -DBL_MAX;
    }



#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
    std::cout << " in values " << std::endl;
    // for( unsigned int i=0; i<var._nvar; i++ ){
    //   if( var2._lst[i].empty() ) continue;
    //   if(var2._shadow_info[1]>0) std::cout << "shadowO " << i << ": "<< var2._shadow[i].oveEst << std::endl;
    //   std::cout << "lstO " << i << ": "<< var2._lst[i].oveEst << std::endl;   
    //   if(var2._shadow_info[0]>0) std::cout << "shadowU " << i << ": "<< var2._shadow[i].undEst << std::endl;
    //   std::cout << "lstU " << i << ": "<< var2._lst[i].undEst << std::endl;         
    // }
    std::cout << "_asym_relu_withShadow: UND vs OVE " << var2._shadow_info[0] << " " << var2._shadow_info[1] << std::endl;
#endif    
    var2._mod->_asym_relu_withShadow( var2._lst, var2._ndep, var2._shadow, var2._shadow_info);
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
    std::cout << "_asym_relu_withShadow: UND vs OVE " << var2._shadow_info[0] << " " << var2._shadow_info[1] << std::endl;
    std::cout << " out values " << std::endl;
    // for( unsigned int i=0; i<var._nvar; i++ ){  
    //   if( var2._lst[i].empty() ) continue;
    //   if(var2._shadow_info[1]>0) std::cout << "shadowO " << i << ": "<< var2._shadow[i].oveEst << std::endl;
    //   std::cout << "lstO " << i << ": "<< var2._lst[i].oveEst << std::endl;   
    //   if(var2._shadow_info[0]>0) std::cout << "shadowU " << i << ": "<< var2._shadow[i].undEst << std::endl;
    //   std::cout << "lstU " << i << ": "<< var2._lst[i].undEst << std::endl;         
    // }  
#endif
    var2._bnd.second = false;
//    var2.get_C();
#ifdef MC__ASM_DEBUG_TRACE 
   //std::cout << "end relu" << std::endl;
#endif    
    return var2;
  } 

  //var._mod->_IntmdtCntnrSeted = false;
  T lazyBnd = var.bound();
  //var._mod->_IntmdtCntnrSeted = false;
  if(Op<T>::l(lazyBnd) > -MC__ASM_COMPUTATION_TOL){
    ASVar<T> var2( var );
    var2._bnd.second = false;
    return var2;
  } 
  
  if(Op<T>::u(lazyBnd) < MC__ASM_COMPUTATION_TOL){
    ASVar<T> var2( 0. );  // note that the _bnd will be set to (0,false) in the constructor
    var2._bnd.second = false;
    return var2;
  }  

  //std::cout << "  _relu input" << std::endl;
  ASVar<T> var2( var ); 
    if(var2._lnr.second){
    var2._ndiv = var2._mod->_ndiv;
    var2._nvar = var2._mod->_nvar;
    var2._lst.resize(var2._nvar);
    var2._ndep = 0;
    for(unsigned int i = 0; i < var2._nvar; i++){
      if(var2._lnr.first[i] == 0. ) continue;
      var2._lst[i] = UnivarPWL((var2._mod->_bndvar[i]),var2._lnr.first[i]); 
      var2._ndep ++; 
    }
    if(var2._cst != 0){
      double __cst = var2._cst/ (double) var2._ndep;    
      for(unsigned int i = 0; i < var2._nvar; i++){
        if(var2._lnr.first[i] == 0. ) continue;
        var2._lst[i] += __cst; 
      }    
    }
    var2._lnr.second = false;
    var2._lnr.first.resize(0);
  }

  var2._mod->_asym_relu( var2._lst,var2._ndep);
  var2._bnd.second = false;
  return var2;
  

}

template <typename T>
inline
ASVar<T> relu
( ASVar<T> && var )
{
#ifdef MC__ASM_DEBUG_TRACE  
  std::cout << "  SD relu &&" << std::endl;	
#endif
  if( !var._mod )
    return std::max(var._cst,0.);

  if (var._mod->options.SHADOW_USE)
  {
    //std::cout << "in relu" << std::endl;     
    var._mod->_IntmdtCntnrSeted = false;
    T lazyBnd = var.bound();
    if(Op<T>::u(lazyBnd) < MC__ASM_COMPUTATION_TOL){
      ASVar<T> var2( 0. );  // note that the _bnd will be set to (0,false) in the constructor
      //std::cout << "end relu" << std::endl; 
      return var2;
    }  

    // if(Op<T>::l(lazyBnd) > -MC__ASM_COMPUTATION_TOL){
    //   var._bnd.second = false;
    //   var._mod->_IntmdtCntnrSeted = false;
    //   return var;
    // } 
    
    if(var._lnr.second){
      var._ndiv = var._mod->_ndiv;
      var._nvar = var._mod->_nvar;
      var._lst.resize(var._nvar);
      var._ndep = 0;
      for(unsigned int i = 0; i < var._nvar; i++){
        if(var._lnr.first[i] == 0.) continue;
        var._lst[i] = UnivarPWL((var._mod->_bndvar[i]),var._lnr.first[i]); 
        var._ndep ++; 
      }
      if(var._cst != 0){
        double __cst = var._cst/ (double) var._ndep;    
        for(unsigned int i = 0; i < var._nvar; i++){
          if(var._lnr.first[i] == 0.) continue;
          var._lst[i] += __cst; 
        }    
      }
      var._lnr.second = false;
      var._lnr.first.resize(0);
    
      var._shadow_info.resize(3,0.);
      var._shadow.resize(var._nvar);
      var._oveCut = DBL_MAX;
      var._undCut = -DBL_MAX;
    }


#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
    std::cout << " in values " << std::endl;
    // for( unsigned int i=0; i<var._nvar; i++ ){
    //   if( var._lst[i].empty() ) continue;
    //   if(var._shadow_info[1]>0) std::cout << "shadowO " << i << ": "<< var._shadow[i].oveEst << std::endl;
    //   std::cout << "lstO " << i << ": "<< var._lst[i].oveEst << std::endl;   
    //   if(var._shadow_info[0]>0) std::cout << "shadowU " << i << ": "<< var._shadow[i].undEst << std::endl;
    //   std::cout << "lstU " << i << ": "<< var._lst[i].undEst << std::endl;         
    // }
    std::cout << "_asym_relu_withShadow: UND vs OVE " << var._shadow_info[0] << " " << var._shadow_info[1] << std::endl;
#endif
    var._mod->_asym_relu_withShadow( var._lst, var._ndep, var._shadow, var._shadow_info);
#ifndef MC__ASMODEL_NOT_DEBUG_SHADOW  
    std::cout << "_asym_relu_withShadow: UND vs OVE " << var._shadow_info[0] << " " << var._shadow_info[1] << std::endl;
    std::cout << " out values " << std::endl;
    //for( unsigned int i=0; i<var._nvar; i++ ){
    //   if( var._lst[i].empty() ) continue;
    //   if(var._shadow_info[1]>0) std::cout << "shadowO " << i << ": "<< var._shadow[i].oveEst << std::endl;
    //   std::cout << "lstO " << i << ": "<< var._lst[i].oveEst << std::endl;   
    //   if(var._shadow_info[0]>0) std::cout << "shadowU " << i << ": "<< var._shadow[i].undEst << std::endl;
    //   std::cout << "lstU " << i << ": "<< var._lst[i].undEst << std::endl;         
    //}  
#endif
    var._bnd.second = false;
#ifdef MC__ASM_DEBUG_TRACE 
   //std::cout << "end relu" << std::endl;
#endif
   //var.get_C();    
    return var;
  }
  //var._mod->_IntmdtCntnrSeted = false;
  T lazyBnd = var.bound();
  //var._mod->_IntmdtCntnrSeted = false;
  if(Op<T>::l(lazyBnd)> -MC__ASM_COMPUTATION_TOL){
    var._bnd.second = false;
    return var;
  } 
  
  if(Op<T>::u(lazyBnd) < MC__ASM_COMPUTATION_TOL){
    var = 0.; 
    //ASVar<T> var2( 0. );
    var._bnd.second = false;
    return var;
  }  

    if(var._lnr.second){
    var._ndiv = var._mod->_ndiv;
    var._nvar = var._mod->_nvar;
    var._lst.resize(var._nvar);
    var._ndep = 0;
    for(unsigned int i = 0; i < var._nvar; i++){
      if(var._lnr.first[i] == 0.) continue;
      var._lst[i] = UnivarPWL((var._mod->_bndvar[i]),var._lnr.first[i]); 
      var._ndep ++; 
    }
    if(var._cst != 0){
      double __cst = var._cst/ (double) var._ndep;    
      for(unsigned int i = 0; i < var._nvar ; i++){
        if(var._lnr.first[i] == 0.) continue;
        var._lst[i] += __cst; 
      }    
    }
    var._lnr.second = false;
    var._lnr.first.resize(0);
  }


  //std::cout << "  _relu input" << std::endl;
  var._mod->_asym_relu( var._lst,var._ndep);
  var._bnd.second = false;
  return var;

}

template <typename T>
inline
ASVar<T>
max
( ASVar<T> const& var1, ASVar<T> const& var2 )
{
  if( !var1._mod && !var2._mod ) 
    return std::max( var1._cst, var2._cst );
  if( !var2._mod ) 
    return max( var1, var2._cst );
  if( !var1._mod ) 
    return max( var2, var1._cst );
  return max( var1 - var2, 0. ) + var2;
}

template <typename T>
inline
ASVar<T>
max
( ASVar<T> const& var1, double const& cst2 )
{
  //std::cout<< "max" << std::endl;
  if(isequal(cst2,0.)){
    return relu(var1);
  }  

  if( !var1._mod ) 
    return std::max( var1._cst, cst2 );
  //std::cout << "exec max:" << cst2 <<std::endl;


/*  T lazyBnd = var1.bound();
  
  if(Op<T>::l(lazyBnd)> cst2){
    ASVar<T> var2( var1 );
    var2._bnd.second = false;
    return var2;
  } 
  
  if(Op<T>::u(lazyBnd) < cst2){
    ASVar<T> var2( cst2 );
    return var2;
  }  


  ASVar<T> var2( var1 );
  auto const& f = [=]( const double& x ){return std::max( x, cst2 ); };//Op<T>::max( x, cst2 ); }; // why not  return std::max( x, cst2 ); };
  if (var2._mod->options.SLOPE_USE ){ 
    auto const& fDerv = [=]( const double& x ){ return x > cst2? 1.: 0.; };  
    var2._mod->_asym_slope( var2._mat, var2._slope, var2._mod->_psize, var2._ndep, f,fDerv, -DBL_MAX, true );
    var2._bnd.second = false;
    return var2;
  }           
  var2._mod->_asym( var2._mat, var2._ndep, f, cst2, true ); // convex term, min cst2
  var2._bnd.second = false;
  return var2;
*/
  throw typename mc::ASModel<T>::Exceptions( mc::ASModel<T>::Exceptions::UNDEF );
  return var1;
}

template <typename T>
inline
ASVar<T>
max
( ASVar<T> && var1, double const& cst2 )
{
  //std::cout<< "max" << std::endl;

  if(isequal(cst2,0.)){
    return relu(std::move(var1));
  }

  if( !var1._mod ) 
    return std::max( var1._cst, cst2 );

/*  T lazyBnd = var1.bound();
  
  if(Op<T>::l(lazyBnd)> cst2){
    var1._bnd.second = false;
    return var1;
  } 
  
  if(Op<T>::u(lazyBnd) < cst2){
    var1 = cst2;
    return var1;
  }  

  auto const& f = [=]( const double& x ){ return std::max( x, cst2 ); };
  if (var1._mod->options.SLOPE_USE){
    auto const& fDerv = [=]( const double& x ){ return x > cst2? 1.: 0.; }; 
    var1._mod->_asym_slope( var1._mat, var1._slope, var1._mod->_psize, var1._ndep, f,fDerv, -DBL_MAX, true );
    var1._bnd.second = false;
    return var1;
  }  
  var1._mod->_asym( var1._mat, var1._ndep, f, cst2, true ); // convex term, min cst2
  var1._bnd.second = false;
  return var1;
*/
  throw typename mc::ASModel<T>::Exceptions( mc::ASModel<T>::Exceptions::UNDEF );
  return var1;
}

template <typename T>
inline
ASVar<T>
max
( double const& cst1, ASVar<T> const& var2 )
{
  return max( var2, cst1 );
}

template <typename T>
inline
ASVar<T>
max
( double const& cst1, ASVar<T> && var2 )
{
  return max( var2, cst1 );
}



}//namespace mc
#include "mcfadbad.hpp"

namespace fadbad
{

//! @brief Specialization of the structure fadbad::Op for use of the type mc::ASVar of MC++ as a template parameter of the classes fadbad::F, fadbad::B and fadbad::T of FADBAD++
template< typename T > struct Op< mc::ASVar<T> >
{ 
  typedef mc::ASVar<T> ASV;
  typedef double Base;
  static Base myInteger( const int i ) { return Base(i); }
  static Base myZero() { return myInteger(0); }
  static Base myOne() { return myInteger(1);}
  static Base myTwo() { return myInteger(2); }
  static double myPI() { return mc::PI; }
  static ASV myPos( const ASV& x ) { return  x; }
  static ASV myNeg( const ASV& x ) { return -x; }
  template <typename U> static ASV& myCadd( ASV& x, const U& y ) { return x+=y; }
  template <typename U> static ASV& myCsub( ASV& x, const U& y ) { return x-=y; }
  template <typename U> static ASV& myCmul( ASV& x, const U& y ) { return x*=y; }
  template <typename U> static ASV& myCdiv( ASV& x, const U& y ) { return x/=y; }
  static ASV myInv( const ASV& x ) { return mc::inv( x ); }
  static ASV mySqr( const ASV& x ) { return mc::sqr( x ); }
  template <typename X, typename Y> static ASV myPow( const X& x, const Y& y ) { return mc::pow( x, y ); }
  //static ASV myCheb( const ASV& x, const unsigned n ) { return mc::cheb( x, n ); }
  static ASV mySqrt( const ASV& x ) { throw typename mc::ASModel<T>::Exceptions( mc::ASModel<T>::Exceptions::UNDEF );  }//{ return mc::sqrt(x); }
  static ASV myLog( const ASV& x ) { throw typename mc::ASModel<T>::Exceptions( mc::ASModel<T>::Exceptions::UNDEF );  }//{ return mc::log( x ); }
  static ASV myExp( const ASV& x ) { throw typename mc::ASModel<T>::Exceptions( mc::ASModel<T>::Exceptions::UNDEF );  }//{ return mc::exp( x ); }
  static ASV mySin( const ASV& x ) { throw typename mc::ASModel<T>::Exceptions( mc::ASModel<T>::Exceptions::UNDEF );  }//{ return mc::sin( x ); }
  static ASV myCos( const ASV& x ) { throw typename mc::ASModel<T>::Exceptions( mc::ASModel<T>::Exceptions::UNDEF ); }//{ return mc::cos( x ); }
  static ASV myTan( const ASV& x ) { throw typename mc::ASModel<T>::Exceptions( mc::ASModel<T>::Exceptions::UNDEF ); } //{ return mc::tan( x ); }
  static ASV myAsin( const ASV& x ) { throw typename mc::ASModel<T>::Exceptions( mc::ASModel<T>::Exceptions::UNDEF ); }
  static ASV myAcos( const ASV& x ) { throw typename mc::ASModel<T>::Exceptions( mc::ASModel<T>::Exceptions::UNDEF ); }
  static ASV myAtan( const ASV& x ) { throw typename mc::ASModel<T>::Exceptions( mc::ASModel<T>::Exceptions::UNDEF ); }
  static ASV mySinh( const ASV& x ) { throw typename mc::ASModel<T>::Exceptions( mc::ASModel<T>::Exceptions::UNDEF ); }
  static ASV myCosh( const ASV& x ) { throw typename mc::ASModel<T>::Exceptions( mc::ASModel<T>::Exceptions::UNDEF ); }
  static ASV myTanh( const ASV& x ) { throw typename mc::ASModel<T>::Exceptions( mc::ASModel<T>::Exceptions::UNDEF ); }
  static bool myEq( const ASV& x, const ASV& y ) { return mc::Op<T>::eq(const_cast<ASV*>(&x)->bound(),const_cast<ASV*>(&y)->bound()); } 
  static bool myNe( const ASV& x, const ASV& y ) { return mc::Op<T>::ne(const_cast<ASV*>(&x)->bound(),const_cast<ASV*>(&y)->bound()); }
  static bool myLt( const ASV& x, const ASV& y ) { return mc::Op<T>::lt(const_cast<ASV*>(&x)->bound(),const_cast<ASV*>(&y)->bound()); }
  static bool myLe( const ASV& x, const ASV& y ) { return mc::Op<T>::le(const_cast<ASV*>(&x)->bound(),const_cast<ASV*>(&y)->bound()); }
  static bool myGt( const ASV& x, const ASV& y ) { return mc::Op<T>::gt(const_cast<ASV*>(&x)->bound(),const_cast<ASV*>(&y)->bound()); }
  static bool myGe( const ASV& x, const ASV& y ) { return mc::Op<T>::ge(const_cast<ASV*>(&x)->bound(),const_cast<ASV*>(&y)->bound()); }
};

} // end namespace fadbad

namespace mc
{

//! @brief C++ structure for specialization of the mc::Op templated structure for use of mc::ASVar in DAG evaluation and as template parameter in other MC++ types
template< typename T > struct Op< mc::ASVar<T> >
{
  typedef mc::ASVar<T> ASV;
  static ASV point( const double c ) { return ASV(c); }
  static ASV zeroone() { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } //{ return ASV( mc::Op<T>::zeroone() ); }
  static void I(ASV& x, const ASV&y) { x = y; }
  static double l(const ASV& x) { return mc::Op<T>::l(x.B()); }
  static double u(const ASV& x) { return mc::Op<T>::u(x.B()); }
  static double abs (const ASV& x) { return mc::Op<T>::abs(x.B());  }
  static double mid (const ASV& x) { return mc::Op<T>::mid(x.B());  }
  static double diam(const ASV& x) { return mc::Op<T>::diam(x.B()); }
  //static ASV intsct(const ASV& x, double l, double u) {return mc::intersect(x,T(l,u));}
  static ASV inv (const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );  }
  static ASV sqr (const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );  }
  static ASV sqrt(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );}
  static ASV exp (const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );  }
  static ASV log (const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );  }
  static ASV xlog(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV lmtd(const ASV& x, const ASV& y) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV rlmtd(const ASV& x, const ASV& y) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV fabs(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV sin (const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );  }
  static ASV cos (const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );  }
  static ASV tan (const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } //{ return mc::tan(x);  }
  static ASV asin(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV acos(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV atan(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV sinh(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV cosh(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV tanh(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } //{ ASV y(std::move(x)); return mc::tanh(std::move(y));} //{ ASV y(std::move(x)); return mc::tanh(std::move(y));} // { return mc::tanh(x); } //
  static ASV tanh(ASV&& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } //{ ASV y(std::move(x)); return mc::tanh(std::move(y));} //{ ASV y(std::move(x)); return mc::tanh(std::move(y));} // { return mc::tanh(x); } //
  static ASV erf (const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV erfc(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV fstep(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV bstep(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV hull(const ASV& x, const ASV& y) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } //{ return mc::hull(x,y); }

  static ASV min (const ASV& x, const ASV& y) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  //static ASV max (const ASV& x, const double & y) { std::cout << "in a func" << std::endl;y == 0.?mc::relu(x):throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }  
  //static ASV max (const ASV& x, const ASV& y) {  std::cout << "in b func" << std::endl; throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV max (const ASV& x, const ASV& y) {  return mc::max(x,y); }
  static ASV arh (const ASV& x, const double k) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  template <typename X, typename Y> static ASV pow(const X& x, const Y& y) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV cheb(const ASV& x, const unsigned n) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV prod (const unsigned n, const ASV* x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static ASV monom (const unsigned n, const ASV* x, const unsigned* k) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  static bool inter(ASV& xIy, const ASV& x, const ASV& y) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } //{ return mc::inter(xIy,x,y); }
  static bool eq(const ASV& x, const ASV& y) { return mc::Op<T>::eq(x.B(),y.B()); }
  static bool ne(const ASV& x, const ASV& y) { return mc::Op<T>::ne(x.B(),y.B()); }
  static bool lt(const ASV& x, const ASV& y) { return mc::Op<T>::lt(x.B(),y.B()); }
  static bool le(const ASV& x, const ASV& y) { return mc::Op<T>::le(x.B(),y.B()); }
  static bool gt(const ASV& x, const ASV& y) { return mc::Op<T>::gt(x.B(),y.B()); }
  static bool ge(const ASV& x, const ASV& y) { return mc::Op<T>::ge(x.B(),y.B()); }

  // The following functions haven't been implemented in any version of ASModel.hpp by Jul-2022.
  // They are listed here merely for the compatibility with MAiNGO solver, to allow the relevent computation in LBP to be based on McCormick<ASV> instead of McCormick<Interval>.
  // McCormick<ASV> instances are McCormick relaxations whose supporting interval bounds are replaced with ASMs theryby got enhanced by ISA. 
  // static ASV fabsx_times_x(const ASV& x) {throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );} //{ return mc::fabsx_times_x(x); }
  // static ASV xexpax(const ASV& x, const double a) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } //{ return mc::xexpax(x,a); }
  // static ASV centerline_deficit(const ASV& x, const double xLim, const double type) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } //{ return mc::centerline_deficit(x,xLim,type); }
  // static ASV wake_profile(const ASV& x, const double type) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } //{ return mc::wake_profile(x,type); }
  // static ASV wake_deficit(const ASV& x, const ASV& r, const double a, const double alpha, const double rr, const double type1, const double type2) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::wake_deficit(x,r,a,alpha,rr,type1,type2); }
  // static ASV power_curve(const ASV& x, const double type) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::power_curve(x,type); }
  // static ASV mid(const ASV& x, const ASV& y, const double k) { return mc::mid(x, y, k); }
  // static ASV pinch(const ASV& Th, const ASV& Tc, const ASV& Tp) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } //{ return mc::pinch(Th, Tc, Tp); }
  // static ASV euclidean_norm_2d(const ASV& x,const ASV& y) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::euclidean_norm_2d(x,y); }
  // static ASV expx_times_y(const ASV& x,const ASV& y) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::expx_times_y(x,y); }
  // static ASV vapor_pressure(const ASV& x, const double type, const double p1, const double p2, const double p3, const double p4 = 0, const double p5 = 0, const double p6 = 0,
	// 						const double p7 = 0, const double p8 = 0, const double p9 = 0, const double p10 = 0) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } //{ return mc::vapor_pressure(x,type,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10);}
  // static ASV ideal_gas_enthalpy(const ASV& x, const double x0, const double type, const double p1, const double p2, const double p3, const double p4, const double p5, const double p6 = 0,
	// 						   const double p7 = 0) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } //{ return mc::ideal_gas_enthalpy(x,x0,type,p1,p2,p3,p4,p5,p6,p7);}
  // static ASV saturation_temperature(const ASV& x, const double type, const double p1, const double p2, const double p3, const double p4 = 0, const double p5 = 0, const double p6 = 0,
	// 								const double p7 = 0, const double p8 = 0, const double p9 = 0, const double p10 = 0) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::saturation_temperature(x,type,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10);}
  // static ASV enthalpy_of_vaporization(const ASV& x, const double type, const double p1, const double p2, const double p3, const double p4, const double p5, const double p6 = 0) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } //{ return mc::enthalpy_of_vaporization(x,type,p1,p2,p3,p4,p5,p6);}
  // static ASV cost_function(const ASV& x, const double type, const double p1, const double p2, const double p3) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::cost_function(x,type,p1,p2,p3);}
  // static ASV sum_div(const std::vector< ASV > &x, const std::vector<double> &coeff) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::sum_div(x,coeff);  }
  // static ASV xlog_sum(const std::vector< ASV > &x, const std::vector<double> &coeff) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::xlog_sum(x,coeff);  }  
  // static ASV affine_transform(const std::vector< ASV > &x, const std::vector<double> &weights, const double bias) { return mc::affine_transform(x,weights,bias); }  
  // static ASV affine_transform(std::vector< ASV > &&x, const std::vector<double> &weights, const double bias) { return mc::affine_transform(std::forward<std::vector< ASV > >(x),weights,bias); }  
  // static ASV affine_transform(std::vector< std::reference_wrapper<ASV>> x, const std::vector<double> &weights, const double bias) { return mc::affine_transform(x,weights,bias); }  
  // static ASV nrtl_tau(const ASV& x, const double a, const double b, const double e, const double f) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::nrtl_tau(x,a,b,e,f);}
  // static ASV nrtl_dtau(const ASV& x, const double b, const double e, const double f) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::nrtl_dtau(x,b,e,f);}
  // static ASV nrtl_G(const ASV& x, const double a, const double b, const double e, const double f, const double alpha ) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::nrtl_G(x,a,b,e,f,alpha);}
  // static ASV nrtl_Gtau(const ASV& x, const double a, const double b, const double e, const double f, const double alpha) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::nrtl_Gtau(x,a,b,e,f,alpha);}
  // static ASV nrtl_Gdtau(const ASV& x, const double a, const double b, const double e, const double f, const double alpha) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::nrtl_Gdtau(x,a,b,e,f,alpha);}
  // static ASV nrtl_dGtau(const ASV& x, const double a, const double b, const double e, const double f, const double alpha) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::nrtl_dGtau(x,a,b,e,f,alpha);}
  // static ASV iapws(const ASV& x, const double type) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } //{ return mc::iapws(x,type); }
  // static ASV iapws(const ASV& x, const ASV& y, const double type) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::iapws(x,y,type); }
  // static ASV p_sat_ethanol_schroeder(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::p_sat_ethanol_schroeder(x); }
  // static ASV rho_vap_sat_ethanol_schroeder(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::rho_vap_sat_ethanol_schroeder(x); }
  // static ASV rho_liq_sat_ethanol_schroeder(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::rho_liq_sat_ethanol_schroeder(x); }
  // static ASV covariance_function(const ASV& x, const double type) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::covariance_function(x,type); }
  // static ASV acquisition_function(const ASV& x, const ASV& y, const double type, const double fmin) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::acquisition_function(x,y,type,fmin); }
  // static ASV gaussian_probability_density_function(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::gaussian_probability_density_function(x); }
  // static ASV regnormal(const ASV& x, const double a, const double b) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::regnormal(x,a,b); }
  // static ASV squash_node (const ASV& x, const double lb, const double ub) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); } // { return mc::squash_node(x,lb,ub);  }

  // static ASV lb_func (const ASV& x, const double lb) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );  }
  // static ASV ub_func (const ASV& x, const double ub) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );  }
  // static ASV mc_print ( const ASV& x, const int number) {throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  // static ASV pos (const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );  }
  // static ASV neg (const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );  }
  // static ASV coth(const ASV& x) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF ); }
  // static ASV bounding_func (const ASV& x, const double lb, const double ub) { throw typename ASModel<T>::Exceptions( ASModel<T>::Exceptions::UNDEF );   }


};

} // namespace mc

#endif

