// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__MCHSL_HPP
#define MC__MCHSL_HPP

#include "ffunc.hpp"
#include "ffdep.hpp"

//#define MC__MC13_DEBUG
//#define MC__MC33_DEBUG
//#define MC__MC13_DISABLE_MC21A

//#ifndef MC__USE_HSL

// HSL functions for block decomposition
extern "C" void mc13d_
  ( const int*, const int*, const int*, const int*, const int*, int*, int*, int*, int* );
extern "C" void mc21a_
  ( const int*, const int*, const int*, const int*, const int*, int*, int*, int* );
extern "C" void mc33ad_
  ( const int*, const int*, const int*, int*, const int*, double*, int*, int*, int*,
    int*, int*, int*, int*, int*, int* );

namespace mc
{

//! @brief C++ class for block decomposition of DAG using HSL
////////////////////////////////////////////////////////////////////////
//! mc::FFHSL is a C++ class to enable block decomposition algorithms
//! implemented in the Harwell scientific library (HSL) for DAGs defined
//! within MC++.
////////////////////////////////////////////////////////////////////////
class FFHSL
////////////////////////////////////////////////////////////////////////
{

public:

  //! @brief Default Constructor
  FFHSL
    ( FFGraph* dag=nullptr )
    : _dag( dag )
    {}
  
  //! @brief Destructor
  virtual ~FFHSL
    ()
    {}

  // Retreive pointer to DAG
  FFGraph* dag
    ()
    const
    { return _dag; };

  //! @brief Set DAG environment
  void set
    ( FFGraph* dag )
    { _dag = dag; }

  //! @brief Perform lower triangular block reaarangement of a square system - the output arguments are the same as in the <a href="http://www.hsl.rl.ac.uk/catalogue/mc13.html">documentation of MC13</a>
  bool MC13
    ( const unsigned nDep, const FFVar*pDep, const FFVar*pIndep,
      int*IPERM, int*IOR, int*IB, int&NB, const bool disp=false,
      std::ostream&os=std::cout );
  //! @brief Perform bordered-block triangular reaarangement of a possibly non-square system - the output arguments are the same as in the <a href="http://www.hsl.rl.ac.uk/catalogue/mc33.html">documentation of MC33</a>
  bool MC33
    ( const unsigned nDep, const FFVar*pDep, const unsigned nIndep,
      const FFVar*pIndep, int*IP, int*IQ, int*IPROF, int*IFLAG, const bool disp=false,
      std::ostream&os=std::cout );

  //! @brief FFHSL Exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for exception handling
    enum TYPE{
      DAG = 1,		//!< DAG not set
      UNDEF = -33 	//!< Undocumented error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr=UNDEF ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case DAG:
        return "mc::FFHSL::Exceptions\t DAG not set";
      case UNDEF:
      default:
        return "mc::FFHSL::Exceptions\t Undocumented error";
      }
    }
  private:
    TYPE _ierr;
  };

protected:

  //! @brief pointer to underlying dag
  FFGraph* _dag;

private:
  //! @brief Private methods to block default compiler methods
  FFHSL
    ( FFHSL const& ) = delete;
  FFHSL& operator=
    ( FFHSL const& ) = delete;
};

////////////////////////////////// FFNum ///////////////////////////////////////

inline bool
FFHSL::MC13
( unsigned const nDep, FFVar const* pDep, FFVar const* pIndep,
  int* IPERM, int* IOR, int* IB, int& NB, bool const disp, std::ostream& os )
{
  // Get list of operations
  std::vector<FFVar> pVar( pIndep, pIndep+nDep );
  std::vector<FFDep> vVar( nDep );
  for( unsigned i=0; i<nDep; i++ ) vVar[i].indep(i);
  auto&& sgDep = _dag->subgraph( nDep, pDep );
  for( auto const& op : sgDep.l_op ){
    // Operation not a variable
    if( op->type != FFOp::VAR ) continue;
    FFVar var = *op->varout.front();
    bool isParam = true;
    // Operation is a current independent
    for( unsigned i=0; isParam && i<nDep; i++ )
      if( var.id() == pIndep[i].id() ) isParam = false;
    if( !isParam ) continue;
    // Add dummy variable for parameter
    pVar.push_back( var );
    vVar.push_back( FFDep() );
  }
  const unsigned nVar = pVar.size();
  std::vector<FFDep> wkDep( sgDep.len_tap );
  std::vector<FFDep> vDep( nDep );
  //for( unsigned i=0; i<nVar; i++ )
  //  std::cout << pVar[i] << " = " << vVar[i] << std::endl;
  _dag->eval( sgDep, wkDep, nDep, pDep, vDep.data(), nVar, pVar.data(), vVar.data() );

  // Populate sparse arrays
  std::vector<int> IP(nDep), LENR(nDep), ICN;
  for( unsigned i=0; i<nDep; i++ ){
    IP[i] = ICN.size()+1;
    LENR[i] = vDep[i].dep().size();
    auto cit = vDep[i].dep().begin();
    for( ; cit != vDep[i].dep().end(); ++cit )
      ICN.push_back( (*cit).first+1 );
  }

  // Make a row permutation to remove nonzeros on diagonal: MC21A
  int N = IP.size(), LICN = ICN.size();
  std::vector<int> IW(4*nDep);
#ifdef MC__MC13_DISABLE_MC21A
  for( unsigned i=0; i<nDep; i++ ) IPERM[i] = i+1;
  bool singDep = false;
#else
  int NUMNZ; 
  mc21a_( &N, ICN.data(), &LICN, IP.data(), LENR.data(), IPERM, &NUMNZ, IW.data() );
  bool singDep = NUMNZ<N? true: false;
#ifdef MC__MC13_DEBUG
  std::cout << "Structural singularity: " << (singDep?'Y':'N') << std::endl;
#endif
  if( singDep ) return false;

  // Permute order of equation system using IPERM (!!Fortran style indices!!)
  ICN.clear();
  for( unsigned i=0; i<nDep; i++ ){
    IP[i] = ICN.size()+1;
    LENR[i] = vDep[IPERM[i]-1].dep().size();
    auto cit = vDep[IPERM[i]-1].dep().begin();
    for( ; cit != vDep[IPERM[i]-1].dep().end(); ++cit )
      ICN.push_back( (*cit).first+1 );
  }
#ifdef MC__MC13_DEBUG
  std::cout << "Row reordering: ";
  for( unsigned i=0; i<nDep; i++) std::cout << " " << IPERM[i];
  std::cout << std::endl;
#endif
#endif

  // Make a block lower-triangular decomposition: MC13D
  mc13d_( &N, ICN.data(), &LICN, IP.data(), LENR.data(), IOR, IB, &NB, IW.data() );
#ifdef MC__MC13_DEBUG
  std::cout << "Number of blocks in permuted matrix: " << NB << std::endl;
  std::cout << "Row/Column reordering: ";
  for( unsigned i=0; i<nDep; i++) std::cout << " " << IOR[i];
  std::cout << std::endl;
#endif

  // Display permuted system structure
  if( disp ){
    std::cout << std::endl << "Number of Blocks: " << NB << std::endl;
    os << "Lower-triangular block structure:" << std::endl
       << std::right << "     ";
    for( unsigned j=0; j<nDep; j++ )
    //  os << " " << std::setw(3) << IOR[j]-1;
      os << " " << std::setw(4) << pIndep[IOR[j]-1];
    os << std::endl;
    for( unsigned i=0; i<nDep; i++ ){
      //os << std::setw(3) << IPERM[IOR[i]-1]-1 << " ";
      os << std::setw(4) << pDep[IPERM[IOR[i]-1]-1] << " ";
      for( unsigned j=0; j<nDep; j++ )
        os << std::setw(3) << " "
           << (vDep[IPERM[IOR[i]-1]-1].dep(IOR[j]-1).first?"X ":"  ");
      os << std::endl;
    }
    os << std::endl;
  }
  return true;
}

inline bool
FFHSL::MC33
( unsigned const nDep, FFVar const* pDep, unsigned const nIndep,
  FFVar const* pIndep, int* IP, int* IQ, int* IPROF, int* IFLAG, bool const disp,
  std::ostream& os )
{
  // Get list of operations
  std::vector<FFVar> pVar( pIndep, pIndep+nIndep );
  std::vector<FFDep> vVar( nIndep );
  for( unsigned i=0; i<nIndep; i++ ) vVar[i].indep(i);
  auto&& sgDep = _dag->subgraph( nDep, pDep );
  for( auto const& op : sgDep.l_op ){
    // Operation not a variable
    if( op->type != FFOp::VAR ) continue;
    FFVar var = *op->varout.front();
    bool isParam = true;
    // Operation is a current independent
    for( unsigned i=0; isParam && i<nIndep; i++ )
      if( var.id() == pIndep[i].id() ) isParam = false;
    if( !isParam ) continue;
    // Add dummy variable for parameter
    pVar.push_back( var );
    vVar.push_back( FFDep() );
  }
  const unsigned nVar = pVar.size();
  std::vector<FFDep> wkDep( sgDep.len_tap );
  std::vector<FFDep> vDep( nDep );
  //for( unsigned i=0; i<nVar; i++ )
  //  std::cout << pVar[i] << " = " << vVar[i] << std::endl;
  _dag->eval( sgDep, wkDep, nDep, pDep, vDep.data(), nVar, pVar.data(), vVar.data() );

  // Populate sparse arrays
  std::vector<int> IRN, JCN;
  std::vector<double> A;
  for( unsigned i=0; i<nDep; i++ ){
    auto cit = vDep[i].dep().begin();
    for( ; cit != vDep[i].dep().end(); ++cit ){
      IRN.push_back( i+1 );
      JCN.push_back( (*cit).first+1 );
      A.push_back( IRN.size() );
    }
  }

  // Make a bordered-block lower-triangular decomposition: MC33A
  const int ITYPE = 5;
  const int M = nDep, N = nIndep, NZI = IRN.size();
  int NZO, IERR;
  std::vector<int> IW(M+N), IW1(9*N+3*M);
  mc33ad_( &M, &N, &NZI, &NZO, &ITYPE, A.data(), IRN.data(), JCN.data(),
           IP, IQ, IPROF, IFLAG, IW.data(), IW1.data(), &IERR );
  if( IERR ) return false;
#ifdef MC__MC33_DEBUG
  std::cout << "Width of bordered block: " << IFLAG[2] << std::endl;
  std::cout << "Row reordering: ";
  for( unsigned i=0; i<nDep; i++) std::cout << " " << IP[i];
  std::cout << "Column reordering: ";
  for( unsigned i=0; i<nIndep; i++) std::cout << " " << IQ[i];
  std::cout << std::endl;
#endif

  // Display permuted system structure
  if( disp ){
    os << std::endl << "Bordered-block Lower-triangular structure:" << std::endl
       << std::right << "   ";
    for( unsigned j=0; j<nIndep; j++ )
      os << (j+IFLAG[2]-IFLAG[1]?" ":" | ") << std::setw(4) << pIndep[IQ[j]-1];
    os << std::endl;
    for( unsigned i=0; i<nDep; i++ ){
      //os << std::setw(3) << IP[i]-1;
      os << std::setw(4) << pDep[IP[i]-1];
      for( unsigned j=0; j<nIndep; j++ ){
        os << (j+IFLAG[2]-IFLAG[1]?"   ":"|    ")
           << (vDep[IP[i]-1].dep(IQ[j]-1).first?"X ":"  ");
      }
      os << std::endl;
    }
    os /* << "  Border bandwidth:" << IFLAG[2] */ << std::endl;
  }

  return true;
}

} // namespace mc

#endif

