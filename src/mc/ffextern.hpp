// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__FFEXTERN_HPP
#define MC__FFEXTERN_HPP

// #define MC__FFEXTERN_DEBUG
#define MC__FFEXTERN_CHECK

#include <armadillo>
#include <bitset>

#include "ffdep.hpp"
#include "ffexpr.hpp"
#include "ffinv.hpp"
#include "ffunc.hpp"
#include "mccormick.hpp"
#include "mcfunc.hpp"
#include "polimage.hpp"
#include "pwcu.hpp"
#include "pwlu.hpp"
#include "scmodel.hpp"
#include "slift.hpp"
#include "specbnd.hpp"
#include "supmodel.hpp"

namespace mc
{

//! @brief C++ class defining external DAG operations in MC++ for template-type
//! objects
////////////////////////////////////////////////////////////////////////
//! mc::FFEXTERN is a C++ class for defining an external DAG operations
//! in MC++ for a generic template-type object O. The other template
//! parameter T specifies the type for interval arithmetic.
////////////////////////////////////////////////////////////////////////
template <typename T, typename O>
class FFEXTERN
    ////////////////////////////////////////////////////////////////////////
    : public FFOp
{
 protected:
  // pointer to the object
  O* _ptrObj;
  // whether this class owns the object
  bool _ownObj;
  // whether the internal arrays need resizing
  mutable bool _reset;

  // set tracking any exception thrown during evaluation
  mutable std::set<int> _excpEval;

  //! @brief Storage for floating-point evaluation
  mutable std::vector<double> _DVar;

  //! @brief Storage for Polyhedral relaxation
  mutable PolImg<T>* _POLEnv;
  mutable std::vector<T> _POLRng;
  mutable std::vector<PolVar<T>> _POLVar;
  mutable std::vector<PolVar<T>> _POLRes;
  mutable std::map<PolVar<T> const*, PolVar<T>, lt_PolVar<T>> _POLMap;

  //! @brief Storage for Interval bounds
  mutable std::vector<T> _IVar;
  mutable std::vector<T> _IRes;

  //! @brief Storage for McCormick relaxation
  mutable std::vector<T> _MCRng;
  mutable std::vector<McCormick<T>> _MCVar;
  mutable std::vector<McCormick<T>> _MCRes;

  //! @brief Storage for spectral bounds
  mutable std::vector<T> _SBRng;
  mutable std::vector<std::vector<size_t>> _SBBin;
  mutable std::vector<Specbnd<T>> _SBVar;
  mutable std::vector<Specbnd<T>> _SBRes;
  mutable std::vector<std::vector<double>> _SBCvURef;
  mutable std::vector<std::vector<double>> _SBCvUCoef;
  mutable std::vector<std::vector<double>> _SBCcORef;
  mutable std::vector<std::vector<double>> _SBCcOCoef;

  //! @brief Storage for sparse Chebyshev models
  mutable SCModel<T>* _SCEnv;
  mutable std::vector<T> _SCRng;
  mutable std::vector<SCVar<T>> _SCVar;
  mutable std::vector<SCVar<T>> _SCRes;
  // mutable std::map<t_mon,FFVar,lt_mon>           _SCMon;
  // mutable std::map<t_mon,PolVar<T>,lt_mon>       _SCPOLMon;
  // mutable std::vector<FFVar>                     _SCAux;
  // mutable std::vector<PolVar<T>>                 _SCPOLAux;
  mutable std::vector<PolVar<T>> _SCPOLVScaled;
  mutable std::vector<arma::mat> _SCSlope;
  mutable std::vector<arma::vec> _SCShift;
  mutable std::vector<arma::vec> _SCVert;

  //! @brief Storage for superposition models on fixed partition
  mutable SupModel<PWCU>* _PWCSEnv;
  mutable std::vector<T> _PWCSRng;
  mutable std::vector<SupVar<PWCU>> _PWCSVar;
  mutable std::vector<SupVar<PWCU>> _PWCSRes;
  mutable std::vector<std::vector<PolVar<T>>> _POLPWCSAux;
  mutable std::vector<std::vector<PolVar<T>>> _POLPWCSDis;

  //! @brief Storage for superposition models on adaptive partition
  mutable SupModel<PWLU>* _PWLSEnv;
  mutable std::vector<T> _PWLSRng;
  mutable std::vector<SupVar<PWLU>> _PWLSVar;
  mutable std::vector<SupVar<PWLU>> _PWLSRes;
  mutable std::vector<PolVar<T>> _POLPWLSAux;
  mutable std::vector<double> _DXPWLSAux;
  mutable std::vector<double> _DYPWLSAux;

  //! @brief Resize relaxation containers for the object
  void _resize_relax(PolImg<T>* img) const;

  //! @brief Propagate polyhedral image for the object
  void _propagate_relax(PolImg<T>* img, FFVar** pRes, PolVar<T>* vRes,
                        PolVar<T> const* vVar) const;

  //! @brief Append polyhedral image cuts for the object
  void _backpropagate_relax(PolImg<T>* img, FFOp* pOp, PolVar<T>* vRes,
                            PolVar<T>* vVar) const;

 private:
  //! @brief Initialize LHS sampling
  void _resize_lhs(size_t const nPts, size_t const nDim,
                   std::vector<std::vector<size_t>>& vBin, int seed) const;

  //! @brief Linearize spectral bound relaxation at selected sample points
  void _sample_SBcuts(std::vector<Specbnd<T>> const& vSBVar,
                      std::vector<Specbnd<T>> const& vSBRes,
                      std::vector<std::vector<size_t>> const& vSBBin,
                      std::vector<double>& DVar,
                      std::vector<std::vector<double>>& DCvURef,
                      std::vector<std::vector<double>>& DCvUCoef,
                      std::vector<std::vector<double>>& DCcORef,
                      std::vector<std::vector<double>>& DCcOCoef) const;

  //! @brief Append polyhedral image cuts for spectral bound relaxation
  void _append_SBcuts(PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res,
                      std::set<unsigned> const& depVar, PolVar<T> const* vVar,
                      double& DCvURef, double const* DCvUCoef, double& DCcORef,
                      double const* DCcOCoef) const;

  //! @brief Append polyhedral image cuts for McCormick relaxation
  template <typename U>
  void _append_MCcuts(PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res,
                      McCormick<U> const& MCRes, PolVar<T> const* vVar,
                      std::vector<McCormick<U>> const& vMCVar) const;

  //! @brief Append polyhedral image cuts for variables in piecewise-constant
  //! superposition relaxation
  void _append_PWCScuts(PolImg<T>* img, FFOp* pOp, PolVar<T>& Var,
                        std::vector<PolVar<T>>& vAux, PWCU const& uest,
                        PWCU const& oest) const;

  //! @brief Append polyhedral image cuts for disaggregated variables in
  //! piecewise-constant superposition relaxation
  void _append_PWCScuts(PolImg<T>* img, FFOp* pOp, PolVar<T>& Var,
                        std::vector<PolVar<T>>& vAux,
                        std::vector<PolVar<T>>& vDis, PWCU const& uest,
                        PWCU const& oest) const;

  //! @brief Append polyhedral image cuts for dependent in piecewise-constant
  //! superposition relaxation
  void _append_PWCScuts(PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res,
                        std::vector<std::vector<PolVar<T>>> const& vAux,
                        std::vector<PWCU> const& est,
                        std::set<unsigned int> const& dep,
                        bool const under) const;

  //! @brief Append polyhedral image cuts for dependent in piecewise-constant +
  //! slope superposition relaxation
  void _append_PWCScuts(PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res,
                        PolVar<T> const* vVar,
                        std::vector<std::vector<PolVar<T>>> const& vAux,
                        std::vector<std::vector<PolVar<T>>> const& vDis,
                        std::vector<PWCU> const& uest,
                        std::vector<PWCU> const& oest,
                        std::set<unsigned int> const& dep) const;

  //! @brief Append polyhedral image cuts for variables and dependent in
  //! piecewise-linear superposition relaxation
  void _append_PWLScuts(PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res,
                        PolVar<T>* vVar, std::vector<PolVar<T>>& vAux,
                        std::vector<double>& Xk, std::vector<double>& Yk,
                        std::vector<PWLU> const& est,
                        std::set<unsigned int> const& dep,
                        bool const under) const;

  //! @brief Generate under- nad overestimator hypeprlanes in [0,1] using a
  //! hybrid corner cutting method
  template <typename KEY, typename COMP>
  void _generate_corner_cuts(
      typename SCVar<T, KEY, COMP>::t_poly const& coefbern,
      std::map<KEY, unsigned, COMP> const& degmax,
      std::vector<arma::vec>& verts, arma::mat& slopes,
      arma::vec& shifts) const;

  //! @brief Append polyhedral image cuts for sparse Chebyshev model relaxation
  void _append_SCMcuts(PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res,
                       PolVar<T> const* vVar, std::vector<PolVar<T>>& vAux,
                       arma::mat const& A, arma::vec const& b, T const& R,
                       bool const first) const;

 public:
  //! @brief FFEXTERN options
  struct Options
  {
    //! @brief Relaxation type
    enum RELAX_TYPE
    {
      INT = 0,  //!< Interval bounds
      AUX,      //!< Auxiliary variable polyhedral relaxation
      MC,       //!< McCormick relaxation with interval bounds
      SB,       //!< Convexification using spectral bounds
      SCM,      //!< Sparse Chebyshev model relaxation
      PWCS,     //!< Piecewise-constant superposition bounds
      PWLS      //!< Piecewise-linear superposition bounds
    };

    //! @brief Default constructor
    Options() { reset(); }

    //! @brief Assignment operator
    template <typename OPT>
    Options&
    operator=(OPT const& opt)
    {
      RELAX.clear();
      for (const auto& r : opt.RELAX) RELAX.insert(static_cast<RELAX_TYPE>(r));
      // RELAX      = opt.RELAX;

      POLDEF = opt.POLDEF;
      POLIMG = opt.POLIMG;

      SBLIN  = opt.SBLIN;
      SBSEED = opt.SBSEED;

      SCMODORD  = opt.SCMODORD;
      SCBERNORD = opt.SCBERNORD;
      SCMODEL   = opt.SCMODEL;
      // CMODCUTS
      // CMODPROP
      // CMODDMAX
      // MONSCALE
      // SCQUAD.BASIS
      // MIG_ATOL
      // MIG_RTOL

      PWCDIV    = opt.PWCDIV;
      PWCREL    = opt.PWCREL;
      PWCSLOPE  = opt.PWCSLOPE;
      PWCSHADOW = opt.PWCSHADOW;
      PWCSUP    = opt.PWCSUP;

      PWLINI    = opt.PWLINI;
      PWLMAX    = opt.PWLMAX;
      PWLREL    = opt.PWLREL;
      PWLSHADOW = opt.PWLSHADOW;
      PWLSUP    = opt.PWLSUP;

      return *this;
    }

    //! @brief Reset options
    void
    reset()
    {
      RELAX = {RELAX_TYPE::MC};

      POLDEF = 1;
      POLIMG.reset();

      SBLIN  = 0;
      SBSEED = 42;

      SCMODORD  = 3;
      SCBERNORD = 0;
      SCMODEL.reset();

      PWCDIV    = 16;
      PWCREL    = 0;
      PWCSLOPE  = 1;
      PWCSHADOW = 1;
      PWCSUP.reset();

      PWLINI    = 16;
      PWLMAX    = 16;
      PWLREL    = 0;
      PWLSHADOW = 1;
      PWLSUP.reset();
    }

    //! @brief Type of relaxation
    std::set<RELAX_TYPE> RELAX;
    //! @brief Whether options default to the parent polyhedral relaxation
    //! environment
    bool POLDEF;
    //! @brief Options for polyhedral relaxation environment
    typename PolImg<T>::Options POLIMG;

    //! @brief Number of linearization points of spectral bound relaxation - 0:
    //! N+1 points (centre point + N LHS points)
    size_t SBLIN;
    //! @brief Random number generator seed in LHS sampler for linearization
    //! points of spectral bound relaxation
    int SBSEED;

    //! @brief Maximal order of sparse Chebyshef model
    size_t SCMODORD;
    //! @brief Degree of Bernstein basis conversion
    size_t SCBERNORD;
    //! @brief Options for sparse Chebyshev model environment
    typename SCModel<T>::Options SCMODEL;

    //! @brief Equipartition size in piecewise-constant superposition model
    size_t PWCDIV;
    //! @brief Representation of piecewise-constant univariates - 0: continuous
    //! relaxation; 1: binary encoding
    int PWCREL;
    //! @brief Whether to append cuts from slopes in piecewise-constant
    //! superposition model relaxation
    bool PWCSLOPE;
    //! @brief Whether to append cuts from shadow estimators in
    //! piecewise-constant superposition model relaxation
    bool PWCSHADOW;
    //! @brief Options for superposition model with piecewise-constant
    //! univariate estimators
    typename SupModel<PWCU>::Options PWCSUP;

    //! @brief Initial partition size in piecewise-linear superposition model
    size_t PWLINI;
    //! @brief Maximal partition size in piecewise-linear superposition model
    size_t PWLMAX;
    //! @brief Representation of piecewise-linear univariates - 0: continuous
    //! relaxation; 1: binary encoding; 2: SOS2 encoding
    int PWLREL;
    //! @brief Whether to append cuts from shadow estimators in piecewise-linear
    //! superposition model relaxation
    bool PWLSHADOW;
    //! @brief Options for superposition model with adaptive piecewise-linear
    //! univariate estimators
    typename SupModel<PWLU>::Options PWLSUP;
  } options;

  //! @brief Default Constructor
  FFEXTERN()
      : FFOp(EXTERN),
        _ptrObj(nullptr),
        _ownObj(false),
        _reset(true),
        _POLEnv(nullptr),
        _SCEnv(nullptr),
        _PWCSEnv(nullptr),
        _PWLSEnv(nullptr)
  {
  }

  // Destructor
  virtual ~FFEXTERN()
  {
    if (_ownObj && _ptrObj) delete _ptrObj;
    if (_POLEnv) delete _POLEnv;
    if (_SCEnv) delete _SCEnv;
    if (_PWCSEnv) delete _PWCSEnv;
    if (_PWLSEnv) delete _PWLSEnv;
  }

  // Copy constructor
  FFEXTERN(FFEXTERN<T, O> const& other)
      : FFOp(other),
        _reset(true),
        _POLEnv(nullptr),
        _SCEnv(nullptr),
        _PWCSEnv(nullptr),
        _PWLSEnv(nullptr),
        options(other.options)
  {
#ifdef MC__FFEXTERN_TRACE
    std::cout << "FFEXTERN::copy constructor\n";
#endif
    if (!other._ptrObj)
      throw std::runtime_error(
          "FFEXTERN::copy constructor ** Null pointer to external object\n");

    _ownObj = other._ownObj;
    if (_ownObj)
      _ptrObj = new O(*other._ptrObj);
    else
      _ptrObj = other._ptrObj;
  }
};

template <typename T, typename O>
inline void
FFEXTERN<T, O>::_resize_relax(PolImg<T>* img) const
{
  size_t const nin  = varin.size();
  size_t const nout = varout.size();
  size_t const npts = (!options.SBLIN ? 2 * nin + 1 : options.SBLIN);

  // Used as intermediate to intersect bounds computed with different
  // arithmetics
  // std::cout << "\nnout = " << nout << std::endl;
  _IRes.resize(nout);

  // Set relaxation environment and containers
  for (auto relax : options.RELAX)
  {
    switch (relax)
    {
      case Options::INT:
      default:
        _IVar.resize(nin);
        //_IRes.resize( nout );
        break;

      case Options::AUX:
        if (!_POLEnv) _POLEnv = new PolImg<T>;
        _POLEnv->options = (options.POLDEF ? img->options : options.POLIMG);
        // std::cout << "options.POLDEF = " << options.POLDEF << std::endl;
        // std::cout << "options.POLIMG.ALLOW_NLIN = " <<
        // options.POLIMG.ALLOW_NLIN.size() << std::endl; std::cout <<
        // "options.POLIMG.ALLOW_NLIN = " <<
        // options.POLIMG.ALLOW_NLIN.count(FFOp::TANH) << std::endl;
        _POLVar.resize(nin);
        _POLRes.resize(nout);
        _POLRng.clear();
        break;

      case Options::MC:
        _MCVar.resize(nin);
        _MCRes.resize(nout);
        _MCRng.clear();
        break;

      case Options::SB:
        _SBVar.resize(nin);
        _SBRes.resize(nout);
        _SBCvURef.resize(npts);
        _SBCvUCoef.resize(npts);
        _SBCcORef.resize(npts);
        _SBCcOCoef.resize(npts);
        for (size_t k = 0; k < npts; ++k)
        {
          _SBCvURef[k].resize(nout);
          _SBCvUCoef[k].resize(nin * nout);
          _SBCcORef[k].resize(nout);
          _SBCcOCoef[k].resize(nin * nout);
        }
        _resize_lhs(npts - 1, nin, _SBBin, options.SBSEED);
        _DVar.resize(nin);
        _SBRng.clear();
        break;

      case Options::SCM:
        if (!_SCEnv || _SCEnv->maxord() != options.SCMODORD)
        {
          delete _SCEnv;
          _SCEnv = new SCModel<T>(options.SCMODORD);
        }
        _SCEnv->options = options.SCMODEL;
        _SCVar.resize(nin);
        _SCRes.resize(nout);
        _SCRng.clear();
        break;

      case Options::PWCS:
        if (!_PWCSEnv || _PWCSEnv->nvar() != nin)
        {
          delete _PWCSEnv;
          _PWCSEnv = new SupModel<PWCU>(nin);
        }
        _PWCSEnv->options = options.PWCSUP;
        _PWCSVar.resize(nin);
        _PWCSRes.resize(nout);
        _POLPWCSAux.resize(nin);
        _POLPWCSDis.resize(nin);
        _PWCSRng.clear();
        break;

      case Options::PWLS:
        if (!_PWLSEnv || _PWLSEnv->nvar() != nin)
        {
          delete _PWLSEnv;
          _PWLSEnv = new SupModel<PWLU>(nin);
        }
        _PWLSEnv->options = options.PWLSUP;
        _PWLSVar.resize(nin);
        _PWLSRes.resize(nout);
        _POLPWLSAux.resize(nin);
        _PWLSRng.clear();
        break;
    }
  }
}

template <typename T, typename O>
inline void
FFEXTERN<T, O>::_propagate_relax(PolImg<T>* img, FFVar** pRes, PolVar<T>* vRes,
                                 PolVar<T> const* vVar) const
{
  size_t const nin  = varin.size();
  size_t const nout = varout.size();

  // track any exceptions thrown
  _excpEval.clear();

  // Forward-propagate relaxations through object
  bool firstEval = true, needEval = false;
  for (auto relax : options.RELAX)
  {
    try
    {
      switch (relax)
      {
        default:
          throw std::runtime_error("FFEXTERN::unavailable relaxation");

        // Interval bounds
        case Options::INT:
          for (unsigned i = 0; i < nin; ++i) _IVar[i] = vVar[i].range();
          _ptrObj->eval(_IVar.data(), _IRes.data());
#ifdef MC__FFEXTERN_DEBUG
          for (unsigned j = 0; j < nout; ++j)
            std::cerr << "_IRes[" << j << "]: " << _IRes[j] << std::endl;
#endif
          break;

        // Polyhedral relaxation with auxiliary variables
        case Options::AUX:
          needEval = false;
          for (unsigned i = 0; i < nin; ++i)
          {
            // Check against stored bounds for reevaluation
            if (!_POLRng.empty() && Op<T>::eq(_POLRng[i], vVar[i].range()))
              continue;
            needEval = true;
          }
          if (needEval)
          {
            _POLEnv->reset();
            //_POLEnv->options = (options.POLDEF? img->options: options.POLIMG);
            _POLMap.clear();
            for (unsigned i = 0; i < nin; ++i)
            {
              _POLVar[i].set(_POLEnv, _ptrObj->varin()[i], vVar[i].range(),
                             true);
              _POLMap[&_POLVar[i]] = vVar[i];
            }
            _ptrObj->eval(_POLVar.data(), _POLRes.data());
            // COULD DO CONSTRAINT PROPAGATION BASED ON vRes.range() - NEEDS A
            // PRIORI BoUNDS IN vRes? Update stored bounds
            _POLRng.resize(nin);
            for (unsigned i = 0; i < nin; ++i) _POLRng[i] = vVar[i].range();
            for (unsigned j = 0; j < nout; ++j)
            {
              if (firstEval)
                _IRes[j] = _POLRes[j].range();
              else if (!Op<T>::inter(_IRes[j], _IRes[j], _POLRes[j].range()))
                std::cerr << "Empty intersection of relaxation bounds: "
                          << _IRes[j] << "  " << _POLRes[j].range()
                          << std::endl;
#ifdef MC__FFEXTERN_DEBUG
              std::cerr << "_POLRes[" << j << "]: " << _POLRes[j] << std::endl;
#endif
            }
          }
#ifdef MC__FFEXTERN_DEBUG
          else
            std::cerr << "By-passing forward propagation of polyhedral image"
                      << std::endl;
#endif
          break;

        // McCormick relaxations at mid-point with subgradient in each direction
        case Options::MC:
          needEval = false;
          for (unsigned i = 0; i < nin; ++i)
          {
            // Check against stored bounds for reevaluation
            if (!_MCRng.empty() && Op<T>::eq(_MCRng[i], vVar[i].range()))
              continue;
            _MCVar[i] =
                McCormick<T>(vVar[i].range(), Op<T>::mid(vVar[i].range()))
                    .sub(nin, i);
            needEval = true;
          }
          if (needEval)
          {
            _ptrObj->eval(_MCVar.data(), _MCRes.data());
            // Update stored bounds
            _MCRng.resize(nin);
            for (unsigned i = 0; i < nin; ++i) _MCRng[i] = vVar[i].range();
            for (unsigned j = 0; j < nout; ++j)
            {
              if (firstEval)
                _IRes[j] = _MCRes[j].I();
              else if (!Op<T>::inter(_IRes[j], _IRes[j], _MCRes[j].I()))
                std::cerr << "Empty intersection of relaxation bounds: "
                          << _IRes[j] << "  " << _MCRes[j].I() << std::endl;
#ifdef MC__FFEXTERN_DEBUG
              std::cerr << "_MCRes[" << j << "]: " << _MCRes[j] << std::endl;
#endif
            }
          }
          break;

        // Spectral bounds
        case Options::SB:
          needEval = false;
          for (unsigned i = 0; i < nin; ++i)
          {
            // Check against stored bounds for reevaluation
            if (!_SBRng.empty() && Op<T>::eq(_SBRng[i], vVar[i].range()))
              continue;
            _SBVar[i].set(vVar[i].range(), i, nin);
            needEval = true;
          }
          if (needEval)
          {
            _ptrObj->eval(_SBVar.data(), _SBRes.data());
            // Update stored bounds
            _SBRng.resize(nin);
            for (unsigned i = 0; i < nin; ++i) _SBRng[i] = vVar[i].range();
            for (unsigned j = 0; j < nout; ++j)
            {
              if (firstEval)
                _IRes[j] = _SBRes[j].I();
              else if (!Op<T>::inter(_IRes[j], _IRes[j], _SBRes[j].I()))
                std::cerr << "Empty intersection of relaxation bounds: "
                          << _IRes[j] << "  " << _SBRes[j].I() << std::endl;
#ifdef MC__FFEXTERN_DEBUG
              std::cerr << "_SBRes[" << j << "]: " << _SBRes[j] << std::endl;
#endif
            }
          }
          break;

        // Sparse Chebyshev model
        case Options::SCM:
          _SCEnv->reset_aux();
          needEval = false;
          for (unsigned i = 0; i < nin; ++i)
          {
            // Check against stored bounds for reevaluation
            if (!_SCRng.empty() && Op<T>::eq(_SCRng[i], vVar[i].range()))
              continue;
            _SCVar[i].set(_SCEnv, i, vVar[i].range());
            needEval = true;
          }
          if (needEval)
          {
            _ptrObj->eval(_SCVar.data(), _SCRes.data());
            // Update stored bounds
            _SCRng.resize(nin);
            for (unsigned i = 0; i < nin; ++i) _SCRng[i] = vVar[i].range();
            // Get dependent bounds from sparse Chebyshev model
            for (unsigned j = 0; j < nout; ++j)
            {
              _SCRes[j].simplify();
              if (firstEval)
                _IRes[j] = _SCRes[j].B();
              else if (!Op<T>::inter(_IRes[j], _IRes[j], _SCRes[j].B()))
                std::cerr << "Empty intersection of relaxation bounds: "
                          << _IRes[j] << "  " << _SCRes[j].B() << std::endl;
#ifdef MC__FFEXTERN_DEBUG
              std::cerr << "_SCRes[" << j << "]: " << _SCRes[j] << std::endl;
#endif
            }
          }
          break;

        // Piecewise-constant superposition models
        case Options::PWCS:
          needEval = false;
          for (unsigned i = 0; i < nin; ++i)
          {
            // Check against stored bounds for reevaluation
            if (!_PWCSRng.empty() && Op<T>::eq(_PWCSRng[i], vVar[i].range()))
              continue;
            _PWCSVar[i].set(*_PWCSEnv, i, vVar[i].range(), options.PWCDIV);
            needEval = true;
          }
          if (needEval)
          {
            _ptrObj->eval(_PWCSVar.data(), _PWCSRes.data());
            // Update stored bounds
            _PWCSRng.resize(nin);
            for (unsigned i = 0; i < nin; ++i) _PWCSRng[i] = vVar[i].range();
            for (unsigned j = 0; j < nout; ++j)
            {
              if (firstEval)
                _IRes[j] = T(_PWCSRes[j].l(), _PWCSRes[j].u());
              else if (!Op<T>::inter(_IRes[j], _IRes[j],
                                     T(_PWCSRes[j].l(), _PWCSRes[j].u())))
                std::cerr << "Empty intersection of relaxation bounds: "
                          << _IRes[j] << "  "
                          << T(_PWCSRes[j].l(), _PWCSRes[j].u()) << std::endl;
#ifdef MC__FFEXTERN_DEBUG
              std::cerr << "_PWCSRes[" << j << "]: " << _PWCSRes[j]
                        << std::endl;
#endif
            }
          }
          break;

        // Piecewise-linear superposition models
        case Options::PWLS:
          needEval = false;
          for (unsigned i = 0; i < nin; ++i)
          {
            // Check against stored bounds for reevaluation
            if (!_PWLSRng.empty() && Op<T>::eq(_PWLSRng[i], vVar[i].range()))
              continue;
            _PWLSVar[i].set(*_PWLSEnv, i, vVar[i].range(), options.PWLINI);
            needEval = true;
          }
          if (needEval)
          {
            _ptrObj->eval(_PWLSVar.data(), _PWLSRes.data());
            // Update stored bounds
            _PWLSRng.resize(nin);
            for (unsigned i = 0; i < nin; ++i) _PWLSRng[i] = vVar[i].range();
            for (unsigned j = 0; j < nout; ++j)
            {
              if (firstEval)
                _IRes[j] = T(_PWLSRes[j].l(), _PWLSRes[j].u());
              else if (!Op<T>::inter(_IRes[j], _IRes[j],
                                     T(_PWLSRes[j].l(), _PWLSRes[j].u())))
                std::cerr << "Empty intersection of relaxation bounds: "
                          << _IRes[j] << "  "
                          << T(_PWLSRes[j].l(), _PWLSRes[j].u()) << std::endl;
#ifdef MC__FFEXTERN_DEBUG
              std::cerr << "_PWLSRes[" << j << "]: " << _PWLSRes[j]
                        << std::endl;
#endif
            }
          }
          break;
      }

      firstEval = false;
    }

    catch (...)
    {
      // #ifdef MC__FFEXTERN_DEBUG
      std::cerr << "Set-arithmetic propagation failed" << std::endl;
      // #endif
      _excpEval.insert(relax);
      continue;
    }
  }

  // Set intersected bounds
  for (unsigned j = 0; j < nout; ++j) vRes[j].set(img, *pRes[j], _IRes[j]);

  // Track polyhedral relaxation dependents
  if (options.RELAX.count(Options::AUX))
    for (unsigned j = 0; j < nout; ++j) _POLMap[&_POLRes[j]] = vRes[j];
}

template <typename T, typename O>
inline void
FFEXTERN<T, O>::_backpropagate_relax(PolImg<T>* img, FFOp* pOp, PolVar<T>* vRes,
                                     PolVar<T>* vVar) const
{
  size_t const nin  = varin.size();
  size_t const nout = varout.size();
  size_t const npts = (!options.SBLIN ? 2 * nin + 1 : options.SBLIN);

  // Back-propagate relaxations through object
  for (auto relax : options.RELAX)
  {
    if (_excpEval.count(relax)) break;
    switch (relax)
    {
      // Interval bounds
      case Options::INT:
      default:
        // No cuts apart from interval bounds on vRes
        break;

      // Polyhedral relaxation with auxiliary variables
      case Options::AUX:
        _POLEnv->generate_cuts(nout, _POLRes.data());
#ifdef MC__FFEXTERN_DEBUG
        std::cerr << "_POLEnv:" << *_POLEnv << std::endl;
#endif
        img->insert_cuts(_POLEnv, _POLMap);
        break;

      // McCormick relaxations at mid-point with subgradient in each direction
      case Options::MC:
        for (unsigned j = 0; j < nout; ++j)
        {
          // polyhedral cut generation for MC
          _append_MCcuts(img, pOp, vRes[j], _MCRes[j], vVar, _MCVar);
        }
        break;

      // Spectral bounds
      case Options::SB:
#ifdef MC__FFEXTERN_DEBUG
        for (unsigned j = 0; j < nout; ++j)
          std::cerr << "_SBRes[" << j << "]: " << _SBRes[j] << std::endl;
#endif
        // sample linearization points
        _sample_SBcuts(_SBVar, _SBRes, _SBBin, _DVar, _SBCvURef, _SBCvUCoef,
                       _SBCcORef, _SBCcOCoef);
        // add polyhedral cuts for each output at sampling points
        for (unsigned j = 0; j < nout; ++j)
        {
          for (unsigned k = 0; k < npts; ++k)
          {
            // polyhedral cut generation for SB
            _append_SBcuts(img, pOp, vRes[j], _SBRes[j].D(), vVar,
                           _SBCvURef[k][j], &_SBCvUCoef[k][j * nin],
                           _SBCcORef[k][j], &_SBCcOCoef[k][j * nin]);
            // Append a single SB cut in the case of a linear function
            if (_SBRes[j].N().empty()) break;
          }
        }
        break;

      // Sparse Chebyshev model
      case Options::SCM:
        _SCSlope.resize(nout);
        _SCShift.resize(nout);
        for (unsigned j = 0; j < nout; ++j)
        {
#ifdef MC__FFEXTERN_DEBUG
          std::cerr << "_SCRes[" << j << "]: " << _SCRes[j] << std::endl;
#endif
          unsigned ORD =
              (options.SCBERNORD > options.SCMODORD ? options.SCBERNORD
                                                    : options.SCMODORD);
          auto&& [coefbern, degmax] = _SCRes[j].to_bernstein(ORD);
          // Update polyhedral variable bounds
          if (!Op<T>::inter(_IRes[j], _IRes[j],
                            _SCRes[j].bound_bernstein(coefbern)))
            std::cerr << "Empty intersection of relaxation bounds: " << _IRes[j]
                      << "  " << _SCRes[j].bound_bernstein(coefbern)
                      << std::endl;
          vRes[j].update(_IRes[j]);
#ifdef MC__FFEXTERN_DEBUG
          std::cerr << "bndbern[" << j
                    << "]:" << _SCRes[j].bound_bernstein(coefbern);
          // std::cerr << "coefbern[" << j << "]:" << SCVar<T>::display(
          // coefbern );
          std::cerr << "degmax[" << j << "]:\n";
          for (auto const& [var, deg] : degmax)
            std::cerr << var << ": " << deg << std::endl;
#endif
          _generate_corner_cuts(coefbern, degmax, _SCVert, _SCSlope[j],
                                _SCShift[j]);
#ifdef MC__FFEXTERN_DEBUG
          std::cerr << "_SCVert[" << j << "] =\n";
          for (auto& vert : _SCVert) std::cerr << vert.t();
          std::cerr << "_SCSlope[" << j << "] =\n" << _SCSlope[j];
          std::cerr << "_SCShift[" << j << "] =\n" << _SCShift[j].t();
#endif
          // polyhedral cut generation from Bernstein corner cuts
          _append_SCMcuts(img, pOp, vRes[j], vVar, _SCPOLVScaled, _SCSlope[j],
                          _SCShift[j], _SCRes[j].R(), !j);
        }
        break;

      // Piecewise-constant superposition models
      case Options::PWCS:
        // define auxiliary variables
        for (unsigned i = 0; i < nin; ++i)
        {
          if (!options.PWCSLOPE)
            _append_PWCScuts(img, pOp, vVar[i], _POLPWCSAux[i],
                             _PWCSVar[i].uest()[i], _PWCSVar[i].oest()[i]);
          else
            _append_PWCScuts(img, pOp, vVar[i], _POLPWCSAux[i], _POLPWCSDis[i],
                             _PWCSVar[i].uest()[i], _PWCSVar[i].oest()[i]);
        }

        // add polyhedral cuts
        for (unsigned j = 0; j < nout; ++j)
        {
#ifdef MC__FFEXTERN_DEBUG
          std::cerr << "_PWCSRes[" << j << "] in " << _PWCSRes[j] << std::endl;
#endif
          // constant superposition relaxation
          if (_PWCSRes[j].sdep().empty())
          {
            *img->add_cut(pOp, PolCut<T>::EQ, _PWCSRes[j].cst(), vRes[j], 1.);
            continue;
          }

          // polyhedral cuts from superposition relaxation
          if (!options.PWCSLOPE)
          {
            _append_PWCScuts(img, pOp, vRes[j], _POLPWCSAux,
                             _PWCSRes[j].uest(0), _PWCSRes[j].sdep(), 1);
            _append_PWCScuts(img, pOp, vRes[j], _POLPWCSAux,
                             _PWCSRes[j].oest(0), _PWCSRes[j].sdep(), 0);
            if (options.PWCSHADOW)
            {
              _append_PWCScuts(img, pOp, vRes[j], _POLPWCSAux,
                               _PWCSRes[j].uest(1), _PWCSRes[j].sdep(), 1);
              _append_PWCScuts(img, pOp, vRes[j], _POLPWCSAux,
                               _PWCSRes[j].oest(1), _PWCSRes[j].sdep(), 0);
            }
          }
          else
          {
            _append_PWCScuts(img, pOp, vRes[j], vVar, _POLPWCSAux, _POLPWCSDis,
                             _PWCSRes[j].uest(0), _PWCSRes[j].oest(0),
                             _PWCSRes[j].sdep());
            if (options.PWCSHADOW)
              _append_PWCScuts(img, pOp, vRes[j], vVar, _POLPWCSAux,
                               _POLPWCSDis, _PWCSRes[j].uest(1),
                               _PWCSRes[j].oest(1), _PWCSRes[j].sdep());
          }
        }
        break;

      // Piecewise-linear superposition models
      case Options::PWLS:
        // add polyhedral cuts
        for (unsigned j = 0; j < nout; ++j)
        {
#ifdef MC__FFEXTERN_DEBUG
          std::cerr << "PWLSRes[" << j << "]: " << _PWLSRes[j] << std::endl;
#endif
          // constant superposition relaxation
          if (_PWLSRes[j].sdep().empty())
          {
            *img->add_cut(pOp, PolCut<T>::EQ, _PWLSRes[j].cst(), vRes[j], 1.);
            continue;
          }

          // polyhedral cuts for superposition relaxation
          for (auto& summand : _PWLSRes[j].uest(0))
            summand.reduce(true, options.PWLMAX);
          _append_PWLScuts(img, pOp, vRes[j], vVar, _POLPWLSAux, _DXPWLSAux,
                           _DYPWLSAux, _PWLSRes[j].uest(0), _PWLSRes[j].sdep(),
                           true);
          for (auto& summand : _PWLSRes[j].oest(0))
            summand.reduce(false, options.PWLMAX);
          _append_PWLScuts(img, pOp, vRes[j], vVar, _POLPWLSAux, _DXPWLSAux,
                           _DYPWLSAux, _PWLSRes[j].oest(0), _PWLSRes[j].sdep(),
                           false);
          if (options.PWLSHADOW)
          {
            for (auto& summand : _PWLSRes[j].uest(1))
              summand.reduce(true, options.PWLMAX);
            _append_PWLScuts(img, pOp, vRes[j], vVar, _POLPWLSAux, _DXPWLSAux,
                             _DYPWLSAux, _PWLSRes[j].uest(1),
                             _PWLSRes[j].sdep(), true);
            for (auto& summand : _PWLSRes[j].oest(1))
              summand.reduce(false, options.PWLMAX);
            _append_PWLScuts(img, pOp, vRes[j], vVar, _POLPWLSAux, _DXPWLSAux,
                             _DYPWLSAux, _PWLSRes[j].oest(1),
                             _PWLSRes[j].sdep(), false);
          }
        }
        break;
    }  // end switch
  }

#ifdef MC__FFEXTERN_DEBUG
  std::cerr << *img;
  {
    int dum;
    std::cout << "PAUSED, ENTER 1";
    std::cin >> dum;
  }
#endif
}

template <typename T, typename O>
template <typename U>
inline void
FFEXTERN<T, O>::_append_MCcuts(PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res,
                               McCormick<U> const& MCRes, PolVar<T> const* vVar,
                               std::vector<McCormick<U>> const& vMCVar) const
{
  double rhs1 = -MCRes.cv(), rhs2 = -MCRes.cc();
  for (unsigned i = 0; i < vMCVar.size(); ++i)
  {
    rhs1 += MCRes.cvsub(i) * vMCVar[i].cv();
    rhs2 += MCRes.ccsub(i) * vMCVar[i].cc();
  }
  img->add_cut(pOp, PolCut<T>::LE, rhs1, vMCVar.size(), vVar, MCRes.cvsub(),
               Res, -1.);
  img->add_cut(pOp, PolCut<T>::GE, rhs2, vMCVar.size(), vVar, MCRes.ccsub(),
               Res, -1.);
}

template <typename T, typename O>
inline void
FFEXTERN<T, O>::_resize_lhs(size_t const nPts, size_t const nDim,
                            std::vector<std::vector<size_t>>& vBin,
                            int seed) const
{
  if (!nPts || !nDim || (vBin.size() == nDim && vBin.back().size() == nPts))
    return;

  // Create a "bin" indices vector from 0 to npts-1
  vBin.resize(nDim);
  vBin.back().resize(nPts);
  std::iota(vBin.back().begin(), vBin.back().end(),
            0);  // Fill 0, 1, 2... nPts-1

  // Initialize the Random Engine ONCE outside the loop
  if (seed < 0)
  {
    // std::random_device provides a non-deterministic seed
    std::random_device rng;
    seed = rng();
  }
  // std::mt19937 is a standard mersenne_twister_engine
  std::mt19937 g(seed);

  // Shuffle the indices for each dimension to ensure randomness
  for (size_t i = 0; i < nDim; ++i)
  {
    if (i + 1 < nDim) vBin[i] = vBin.back();
    std::shuffle(vBin[i].begin(), vBin[i].end(), g);
  }
}

template <typename T, typename O>
inline void
FFEXTERN<T, O>::_sample_SBcuts(std::vector<Specbnd<T>> const& vSBVar,
                               std::vector<Specbnd<T>> const& vSBRes,
                               std::vector<std::vector<size_t>> const& vSBBin,
                               std::vector<double>& DVar,
                               std::vector<std::vector<double>>& DCvURef,
                               std::vector<std::vector<double>>& DCvUCoef,
                               std::vector<std::vector<double>>& DCcORef,
                               std::vector<std::vector<double>>& DCcOCoef) const
{
  // Get relaxation values and gradients at sampling points
  size_t const nin  = vSBVar.size();
  size_t const nout = vSBRes.size();
#ifdef MC__FFEXTERN_CHECK
  assert(vSBBin.size());
#endif
  size_t const npts = vSBBin.back().size() + 1;

  for (size_t k = 0; k < npts; ++k)
  {
    // Set linearization points
#ifdef MC__FFEXTERN_DEBUG
    std::cout << "Point " << k << ": [";
#endif
    for (size_t i = 0; i < nin; ++i)
    {
      // First point chosen as centre
      DVar[i] = (k == 0 ? Op<T>::mid(vSBVar[i].I())
                        : Op<T>::l(vSBVar[i].I()) +
                              Op<T>::diam(vSBVar[i].I()) *
                                  ((_SBBin[i][k - 1] + 0.5) / (npts - 1.)));
#ifdef MC__FFEXTERN_DEBUG
      std::cout << " " << DVar[i];
#endif
    }
#ifdef MC__FFEXTERN_DEBUG
    std::cout << " ]" << std::endl;
#endif

    // Evaluate function value and gradient at linearization point
    _ptrObj->grad(DVar.data(), DCvUCoef[k].data(), DCvURef[k].data());
    DCcORef[k]  = DCvURef[k];
    DCcOCoef[k] = DCvUCoef[k];

    // Add spectral bound to convex and concave estimators
    for (size_t j = 0; j < nout; ++j)
    {
      double alphal = -0.5 * Op<T>::l(vSBRes[j].SI()),
             alphau = -0.5 * Op<T>::u(vSBRes[j].SI());
      for (auto const& i : vSBRes[j].D())
        DCvURef[k][j] -= DCvUCoef[k][j * nin + i] * DVar[i];
      if (alphal > 0)
      {
        for (auto const& i : vSBRes[j].N())
        {
          double const corrSlope =
              alphal * (2. * DVar[i] - Op<T>::u(vSBVar[i].I()) -
                        Op<T>::l(vSBVar[i].I()));
          DCvUCoef[k][j * nin + i] += corrSlope;
          DCvURef[k][j] += alphal * (DVar[i] - Op<T>::u(vSBVar[i].I())) *
                               (DVar[i] - Op<T>::l(vSBVar[i].I())) -
                           corrSlope * DVar[i];
        }
      }
      for (auto const& i : vSBRes[j].D())
        DCcORef[k][j] -= DCcOCoef[k][j * nin + i] * DVar[i];
      if (alphau < 0)
      {
        for (auto const& i : vSBRes[j].N())
        {
          double const corrSlope =
              alphau * (2. * DVar[i] - Op<T>::u(vSBVar[i].I()) -
                        Op<T>::l(vSBVar[i].I()));
          DCcOCoef[k][j * nin + i] += corrSlope;
          DCcORef[k][j] += alphau * (DVar[i] - Op<T>::u(vSBVar[i].I())) *
                               (DVar[i] - Op<T>::l(vSBVar[i].I())) -
                           corrSlope * DVar[i];
        }
      }
    }
  }
}

template <typename T, typename O>
inline void
FFEXTERN<T, O>::_append_SBcuts(PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res,
                               std::set<unsigned> const& depVar,
                               PolVar<T> const* vVar, double& DCvURef,
                               double const* DCvUCoef, double& DCcORef,
                               double const* DCcOCoef) const
{
  img->add_cut(pOp, PolCut<T>::LE, -DCvURef, depVar, vVar, DCvUCoef, Res, -1.);
  img->add_cut(pOp, PolCut<T>::GE, -DCcORef, depVar, vVar, DCcOCoef, Res, -1.);
}

template <typename T, typename O>
inline void
FFEXTERN<T, O>::_append_PWCScuts(PolImg<T>* img, FFOp* pOp, PolVar<T>& Var,
                                 std::vector<PolVar<T>>& vAux, PWCU const& uest,
                                 PWCU const& oest) const
{
#ifdef MC__FFEXTERN_CHECK
  assert(uest.yL().size() == options.PWCDIV &&
         oest.yU().size() == options.PWCDIV);
#endif
  vAux.resize(options.PWCDIV);
  for (unsigned k = 0; k < options.PWCDIV; ++k)
    vAux[k].set(img, Op<T>::zeroone(), options.PWCREL ? false : true);
  // sum of auxiliaries equal to 1
  img->add_cut(pOp, PolCut<T>::EQ, 1., options.PWCDIV, vAux.data(), 1.);

  // link auxiliaries to independent variables
  img->add_cut(pOp, PolCut<T>::LE, 0., options.PWCDIV, vAux.data(),
               uest.yL().data(), Var, -1.);
  img->add_cut(pOp, PolCut<T>::GE, 0., options.PWCDIV, vAux.data(),
               oest.yU().data(), Var, -1.);
}

template <typename T, typename O>
inline void
FFEXTERN<T, O>::_append_PWCScuts(PolImg<T>* img, FFOp* pOp, PolVar<T>& Var,
                                 std::vector<PolVar<T>>& vAux,
                                 std::vector<PolVar<T>>& vDis, PWCU const& uest,
                                 PWCU const& oest) const
{
#ifdef MC__FFEXTERN_CHECK
  assert(uest.yL().size() == options.PWCDIV &&
         oest.yU().size() == options.PWCDIV);
#endif
  vAux.resize(options.PWCDIV);
  vDis.resize(options.PWCDIV);
  for (unsigned k = 0; k < options.PWCDIV; ++k)
  {
    vAux[k].set(img, Op<T>::zeroone(), options.PWCREL ? false : true);
    // vDis[k].set( img, T( -1e20, 1e20 ), true );
    vDis[k].set(img, Op<T>::hull(Var.range(), 0.), true);

    // bound disaggregated variable
    img->add_cut(pOp, PolCut<T>::LE, 0., vAux[k], -oest.yU().at(k), vDis[k],
                 1.);
    img->add_cut(pOp, PolCut<T>::GE, 0., vAux[k], -uest.yL().at(k), vDis[k],
                 1.);
  }

  // sum of auxiliaries equal to 1
  img->add_cut(pOp, PolCut<T>::EQ, 1., options.PWCDIV, vAux.data(), 1.);

  // link disaggregated variables to actual variable
  img->add_cut(pOp, PolCut<T>::EQ, 0., options.PWCDIV, vDis.data(), 1., Var,
               -1.);
}

template <typename T, typename O>
inline void
FFEXTERN<T, O>::_append_PWCScuts(
    PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res, PolVar<T> const* vVar,
    std::vector<std::vector<PolVar<T>>> const& vAux,
    std::vector<std::vector<PolVar<T>>> const& vDis,
    std::vector<PWCU> const& uest, std::vector<PWCU> const& oest,
    std::set<unsigned int> const& dep) const
{
  if (dep.empty()) return;

  // polyhedral cuts for superposition relaxation
  auto cutL = (!uest.empty() ? *img->add_cut(pOp, PolCut<T>::LE, 0., Res, -1.)
                             : nullptr);
  auto cutU = (!oest.empty() ? *img->add_cut(pOp, PolCut<T>::GE, 0., Res, -1.)
                             : nullptr);

  for (auto const& i : dep)
  {
    double const xiL = Op<T>::l(vVar[i].range()),
                 dxi = Op<T>::diam(vVar[i].range()) / (double)options.PWCDIV;
    double xiref     = xiL;

    // add cuts for disaggregated variables
    for (unsigned k = 0; k < options.PWCDIV; ++k, xiref += dxi)
    {
      // lower bounding cuts
      if (cutL)
      {
        // PolVar<T> ResikL( img, T( -1e20, 1e20 ), true );
        PolVar<T> ResikL(
            img, T(std::min(uest[i].l(), 0.), std::max(uest[i].u(), 0.)), true);
        cutL->append(ResikL, 1.);
        if (uest[i].yL().size() == options.PWCDIV)
        {
          img->add_cut(pOp, PolCut<T>::GE, 0., ResikL, 1., vAux[i][k],
                       -uest[i].yL().at(k));
        }
        if (uest[i].sLL().size() == options.PWCDIV &&
            uest[i].y0().size() == options.PWCDIV + 1)
        {
          double const siLL = uest[i].sLL().at(k);
          double const biLL = uest[i].y0().at(k) - siLL * xiref;
          img->add_cut(pOp, PolCut<T>::GE, 0., ResikL, 1., vDis[i][k], -siLL,
                       vAux[i][k], -biLL);
        }
        if (uest[i].sLU().size() == options.PWCDIV &&
            uest[i].y0().size() == options.PWCDIV + 1)
        {
          double const siLU = uest[i].sLU().at(k);
          double const biLU = uest[i].y0().at(k + 1) - siLU * (xiref + dxi);
          img->add_cut(pOp, PolCut<T>::GE, 0., ResikL, 1., vDis[i][k], -siLU,
                       vAux[i][k], -biLU);
        }
      }

      // upper bounding cuts
      if (cutU)
      {
        // PolVar<T> ResikU( img, T( -1e20, 1e20 ), true );
        PolVar<T> ResikU(
            img, T(std::min(oest[i].l(), 0.), std::max(oest[i].u(), 0.)), true);
        cutU->append(ResikU, 1.);
        if (oest[i].yU().size() == options.PWCDIV)
        {
          img->add_cut(pOp, PolCut<T>::LE, 0., ResikU, 1., vAux[i][k],
                       -oest[i].yU().at(k));
        }
        if (oest[i].sUL().size() == options.PWCDIV &&
            oest[i].y0().size() == options.PWCDIV + 1)
        {
          double const siUL = oest[i].sUL().at(k);
          double const biUL = oest[i].y0().at(k) - siUL * xiref;
          img->add_cut(pOp, PolCut<T>::LE, 0., ResikU, 1., vDis[i][k], -siUL,
                       vAux[i][k], -biUL);
        }
        if (oest[i].sUU().size() == options.PWCDIV &&
            oest[i].y0().size() == options.PWCDIV + 1)
        {
          double const siUU = oest[i].sUU().at(k);
          double const biUU = oest[i].y0().at(k + 1) - siUU * (xiref + dxi);
          img->add_cut(pOp, PolCut<T>::LE, 0., ResikU, 1., vDis[i][k], -siUU,
                       vAux[i][k], -biUU);
        }
      }
    }
  }
}

template <typename T, typename O>
inline void
FFEXTERN<T, O>::_append_PWCScuts(
    PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res,
    std::vector<std::vector<PolVar<T>>> const& vAux,
    std::vector<PWCU> const& est, std::set<unsigned int> const& dep,
    bool const under) const
{
  if (est.empty()) return;

  // polyhedral cuts for superposition relaxation
  if (under)
  {
    auto cut = *img->add_cut(pOp, PolCut<T>::LE, 0., Res, -1.);
    for (auto const& i : dep)
    {
#ifdef MC__FFEXTERN_CHECK
      assert(est[i].yL().size() == options.PWCDIV);
#endif
      cut->append(options.PWCDIV, vAux[i].data(), est[i].yL().data());
    }
  }
  else
  {
    auto cut = *img->add_cut(pOp, PolCut<T>::GE, 0., Res, -1.);
    for (auto const& i : dep)
    {
#ifdef MC__FFEXTERN_CHECK
      assert(est[i].yU().size() == options.PWCDIV);
#endif
      cut->append(options.PWCDIV, vAux[i].data(), est[i].yU().data());
    }
  }
}

template <typename T, typename O>
inline void
FFEXTERN<T, O>::_append_PWLScuts(
    PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res, PolVar<T>* vVar,
    std::vector<PolVar<T>>& vAux, std::vector<double>& Xk,
    std::vector<double>& Yk, std::vector<PWLU> const& est,
    std::set<unsigned int> const& dep, bool const under) const
{
  if (est.empty()) return;

  for (auto const& i : dep)
  {
    vAux[i].set(img, T(est.at(i).l(), est.at(i).u()),
                true);  // continuous auxiliary

    // Generate breakpoints from segment-based representation in mc::PWLU
    auto const& est_i = est.at(i);
    size_t const NK   = est_i.dx().size() + 1;
    Xk.resize(NK);
    Xk[0] = est_i.xL();
    Yk.resize(NK);
    Yk[0]    = est_i.yL();
    auto idx = est_i.dx().cbegin(), idy = est_i.dy().cbegin();
    for (unsigned k = 1; k < NK; ++k, ++idx, ++idy)
    {
      Xk[k] = Xk[k - 1] + *idx;
      Yk[k] = Yk[k - 1] + *idx * *idy;
    }

    img->add_semilinear_cuts(pOp, NK, vVar[i], Xk.data(), vAux[i], Yk.data(),
                             mc::PolCut<T>::EQ, options.PWLREL);
  }

  img->add_cut(pOp, (under ? PolCut<T>::LE : PolCut<T>::GE), 0., dep,
               vAux.data(), 1., Res, -1.);
}
/*
template< typename T, typename O >
inline FFVar const&
FFEXTERN<T,O>::_get_mon
( t_mon const& mon, int const BASIS, bool const DAGINSERT )
{
  assert( mon.tord );
  auto itXmapmon = _Xmon.find( mon );
  if( itXmapmon != _Xmon.end() ) return itXmapmon->second;

  if( mon.tord == 1 || !DAGINSERT ){
    FFVar Xmon( _dag );
    auto itXmon = _dag->Vars().find( &Xmon );
    _Xmon[mon] = **itXmon;
    return **itXmon;
  }

  _Xaux.resize( mon.expr.size() );
  unsigned nvar = 0;
  for( auto const& [ivar,ord] : mon.expr ){
    // This assumes that variables are already available in _Xmon
    FFVar const& var = _Xmon[t_mon(ivar)];
    switch( BASIS ){
      case t_poly::Options::MONOM: _Xaux[nvar++] = pow( var, (int)ord ); break;
      case t_poly::Options::CHEB:  _Xaux[nvar++] = cheb(var, ord );      break;
    }
  }
  FFVar Xmon;
  switch( BASIS ){
    case t_poly::Options::MONOM: Xmon = FFBase::prod( nvar, _Xaux.data() );
break; case t_poly::Options::CHEB:  Xmon = Op<FFVar>::prod( nvar, _Xaux.data()
); break;
  }
#ifdef MC__MINLPBND_DEBUG_MONDRL
  std::cout << "Subgraph of monomial " << mon.display(BASIS) << std::endl;
  _dag->output( _dag->subgraph( 1, &Xmon ) );
  { int dum; std::cin >> dum; }
#endif
  auto itXmon = _dag->Vars().find( &Xmon );
  _Xmon[mon] = **itXmon;
  return **itXmon;
}

*/
// ---------------------------------------------------------
// Corner Cutting (Graph Bounding / Envelope)
// ---------------------------------------------------------
// Generates a polyhedral relaxation of a polynomial in Bernstein basis with
// coefficients 'coefbern' Polyhedral cuts are in the form Ax <= b, with x =
// [u_0, ..., u_N-1, z]
//
// For each corner c:
// 1. Hyperplane H(u) = S'u
// 2. Max/Min deviations (errors) computed via 'verts'
// 3. Cuts:
//    Row 2c:   z <= S'u + max_dev  =>  -S.u + z <=  K + max_dev
//    Row 2c+1: z >= S'u + min_dev  =>   S.u - z <= -K - min_dev
// ---------------------------------------------------------
template <typename T, typename O>
template <typename KEY, typename COMP>
inline void
FFEXTERN<T, O>::_generate_corner_cuts(
    typename SCVar<T, KEY, COMP>::t_poly const& coefbern,
    std::map<KEY, unsigned, COMP> const& degmax, std::vector<arma::vec>& verts,
    arma::mat& slopes, arma::vec& shifts) const
{
  if (coefbern.empty() || degmax.empty())
  {
    slopes.clear();
    shifts.clear();
    return;
  }

  // Resize variables
  const size_t nVar = degmax.size();
  const size_t nVer = (1 << nVar);
  const size_t nCut = nVer * 2;
  slopes.zeros(nVar + 1, nCut);  // last column correspond to function value 'z'
  shifts.zeros(nCut);
  arma::vec hslope(nVar);

  // Prepare variable keys
  std::vector<KEY> Vars;
  Vars.reserve(nVar);
  for (const auto& [var, deg] : degmax) Vars.push_back(var);

  // Generate (N_points x (d+1)) matrix of control points (u_i, val) from sparse
  // Berstein coefficients
  verts.clear();
  verts.reserve(coefbern.size());
  for (auto const& [mon, coef] : coefbern)
  {
    arma::vec vert(nVar + 1);
    for (unsigned i = 0; i < nVar; ++i)
    {
      unsigned k_i = 0;
      auto it_deg  = mon.expr.find(Vars[i]);
      if (it_deg != mon.expr.end()) k_i = it_deg->second;
      vert(i) = (double)k_i / (double)degmax.at(Vars[i]);
    }
    vert.back() = coef;
    verts.push_back(vert);
  }

  // Helper: Bitmask to Index Map
  auto get_corner_index = [&](size_t mask) -> std::map<KEY, unsigned, COMP>
  {
    std::map<KEY, unsigned, COMP> ndx;
    for (size_t j = 0; j < nVar; ++j)
      ndx[Vars[j]] = ((mask >> j) & 1) ? degmax.at(Vars[j]) : 0;
    return ndx;
  };

  // Process All Corners
  for (size_t c = 0; c < nVer; ++c)
  {
    // 1. Identify Corner Value
    std::map<KEY, unsigned, COMP> cndx = get_corner_index(c);
    double cval                        = 0.0;
    auto it_c                          = coefbern.find(cndx);
    if (it_c != coefbern.end()) cval = it_c->second;

    // 2. Estimate Slopes (Normalized Domain [0,1])
    hslope.zeros();
    for (size_t i = 0; i < nVar; ++i)
    {
      std::map<KEY, unsigned, COMP> nndx = cndx;
      unsigned dmax                      = degmax.at(Vars[i]);
      unsigned dcur                      = cndx[Vars[i]];
      double nval                        = 0.0;

      if (dcur == 0)
      {
        nndx[Vars[i]] = 1;
        auto it_n     = coefbern.find(nndx);
        if (it_n != coefbern.end()) nval = it_n->second;
        hslope(i) = (nval - cval) * dmax;
      }
      else
      {
        nndx[Vars[i]] = dmax - 1;
        auto it_n     = coefbern.find(nndx);
        if (it_n != coefbern.end()) nval = it_n->second;
        hslope(i) = (cval - nval) * dmax;
      }
    }

    // 3. Compute Deviations (Point - Plane)
    // Hyperplane H(u) = cval + Sum( slope_i * (u_i - c_i) )
    double min_dev = 0.0;
    double max_dev = 0.0;

    arma::vec uval(nVar);
    for (size_t i = 0; i < nVar; ++i) uval(i) = ((c >> i) & 1) ? 1.0 : 0.0;

    for (const auto& vert : verts)
    {
      double vval = vert.back();  // The actual z value
      double hval = cval;
      for (size_t i = 0; i < nVar; ++i) hval += hslope(i) * (vert(i) - uval(i));

      double dev = vval - hval;
      if (dev < min_dev) min_dev = dev;
      if (dev > max_dev) max_dev = dev;
    }

    // 4. Fill Matrix A and Vector b
    size_t o = c * 2;      // overestimator column
    size_t u = c * 2 + 1;  // underestimator column

    // We compute the constant part of the plane equation:
    // hshift = cval - Sum(slope_i * uval_i)
    double hshift = cval;
    for (size_t i = 0; i < nVar; ++i) hshift -= hslope(i) * uval(i);

    // Overestimator
    // z <= Sum(s_i * u_i) + hshift + max_dev
    // => -Sum(s_i * u_i) + z <= hshift + max_dev
    for (size_t i = 0; i < nVar; ++i) slopes(i, o) = -hslope(i);
    slopes(nVar, o) = 1.0;  // Coeff for z
    shifts(o)       = hshift + max_dev;

    // Underestimator
    // z >= Sum(s_i * u_i) + K + min_dev
    // => Sum(s_i * u_i) - z <= -K - min_dev
    for (size_t i = 0; i < nVar; ++i) slopes(i, u) = hslope[i];
    slopes(nVar, u) = -1.0;  // Coeff for z
    shifts(u)       = -(hshift + min_dev);
  }
}

template <typename T, typename O>
inline void
FFEXTERN<T, O>::_append_SCMcuts(PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res,
                                PolVar<T> const* vVar,
                                std::vector<PolVar<T>>& vVSca,
                                arma::mat const& A, arma::vec const& b,
                                T const& R, bool const first) const
{
  size_t const nVar = A.n_rows - 1;
  size_t const nCut = A.n_cols;

  // Define auxiliaries for scaled independent variables
  if (first)
  {
    vVSca.resize(nVar);
    for (unsigned i = 0; i < nVar; ++i)
    {
      vVSca[i].set(img, Op<T>::zeroone(), true);
      double const d = Op<T>::diam(vVar[i].range()),
                   l = Op<T>::l(vVar[i].range());
      img->add_cut(pOp, PolCut<T>::EQ, l, vVSca[i], -d, vVar[i], 1.);
    }
  }

  // Add Berstein corner cuts
  for (unsigned k = 0; k < nCut; ++k)
  {
    double const* coefVar = A.colptr(k);
    if (coefVar[nVar] < 0.)
      img->add_cut(pOp, PolCut<T>::LE, b(k) - Op<T>::l(R), nVar, vVSca.data(),
                   coefVar, Res, coefVar[nVar]);
    else
      img->add_cut(pOp, PolCut<T>::LE, b(k) + Op<T>::u(R), nVar, vVSca.data(),
                   coefVar, Res, coefVar[nVar]);
  }
}

}  // end namespace mc

#endif
