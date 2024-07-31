// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_SRED Detection of Reduction Constraints in Sparse Polynomial Constraints
\author Benoit Chachuat
\date 2023

The class mc::SRed defined in <tt>sred.hpp</tt> identifies <b>reduction constraints</b> from a given set of sparse polynomial equality constraints. The starting point is the reduced RLT approach introduced by <A href="https://doi.org/10.1111/j.1475-3995.2004.00438.x">Liberti (2004)</A> and <A href="https://doi.org/10.1007/s10898-006-9005-4">Liberti & Pantelides (2006)</A>, which uses graph theory to identify new quadratic constraints from existing linear equality constraints and participating bilinear terms, namely reduction constraints. The class mc::SRed generalizes this approach in three ways:
- It searches for new polynomial equality constraints based on existing polynomial equality constraints and participating monomial terms;
- It enables the formation of reduction constraint both through multiplication and division of existing polynomial equality constraints with (monomials of) the participating variables;
- It conducts the search via an integer-programming-based approach that supersedes the original graph-theoretic approach by <A href="https://doi.org/10.1007/s10898-006-9005-4">Liberti & Pantelides (2006)</A>.
.
In essence, mc::SRed identifies linearity embedded into polynomial programs and may be applied as a pre-processing step to simplify polynomial programs in practice. 


\section sec_SRED_idea How Does the Detection Algorithm Work?

The starting point is a set of participating monomials, \f$\mathcal{M} := \{ m_i({\bf v}) \}_{1\leq i\leq n_{\rm mon}}\f$ and polynomial equality constraints,
\f{align*}
  \sum_{k=1}^{M_\lambda} a_{\lambda,k}\cdot m_{\lambda,k}({\bf v}) =\ & 0,\ \ \lambda\in 1\ldots n_{\rm pol}
\f}
with monomials \f$m_{\lambda,k}({\bf v})\f$ in the variables \f$v_j\f$, \f$1\leq j\leq n_{\rm var}\f$.

The first step entails identifying of those polynomial constraints which, when multiplied or divided by a given variable (or more generally monomial), lead to the creation of fewer new monomials than extra constraints. <A href="https://doi.org/10.1007/s10898-006-9005-4">Liberti & Pantelides (2006)</A> proposed a graph-theoretical approach to identify a valid assignment. Instead, the class mc::SRed constructs and solves the following integer linear program (ILP) to identify an optimal assignment:
\f{align*}
\min &\ \sum_{j=1}^{n_{\rm var}} \left[ (1+\epsilon)\sum_{i=1}^{n_{\rm mon}} \left(h_{i,j}+l_{i,j}\right) - \sum_{\lambda=1}^{n_{\rm pol}} \left(d_{\lambda,j}+m_{\lambda,j}\right) \right]\\
\text{s.t.}\ &\ h_{i,j} \geq m_{\lambda,j}, \ \ \forall i\in\mathcal{M}_{\lambda,j}, \forall j,\lambda\\
             &\ l_{i,j} \geq d_{\lambda,j}, \ \ \forall i\in\mathcal{D}_{\lambda,j}, \forall j,\lambda\\
             &\ h_{i,j},l_{i,j},m_{\lambda,j},d_{\lambda,j}\in\{0,1\},\ \ \forall i,j,\lambda\,.
\f}

- The binary variables have the following interpretation:
\f{align*}
m_{\lambda,j} :=\ & \left\{\begin{array}{ll} 1 & \text{if polynomial constraint $\lambda\in 1\ldots n_{\rm pol}$ is multiplied by variable $v_j, j\in 1\ldots n_{\rm var}$}\\ 0 & \text{otherwise} \end{array}\right.\\
d_{\lambda,j} :=\ & \left\{\begin{array}{ll} 1 & \text{if polynomial constraint $\lambda\in 1\ldots n_{\rm pol}$ is divided by variable $v_j, j\in 1\ldots n_{\rm var}$}\\ 0 & \text{otherwise} \end{array}\right.\\
h_{i,j} :=\ & \left\{\begin{array}{ll} 1 & \text{if the higher-degree monomial $m_i({\bf v})\cdot v_j$ is created, that is, $m_i({\bf v})\cdot v_j\notin \mathcal{S}_{\rm mon}$}\\ 0 & \text{otherwise} \end{array}\right.\\
l_{i,j} :=\ & \left\{\begin{array}{ll} 1 & \text{if the lower-degree monomial $m_i({\bf v})/ v_j$ is created, that is, $m_i({\bf v})/ v_j\notin \mathcal{S}_{\rm mon}$}\\ 0 & \text{otherwise} \end{array}\right.
\f}

- The index sets \f$\mathcal{M}_{\lambda,j}\f$ and \f$\mathcal{D}_{\lambda,j}\f$ associated with each polynomial constraint \f$\lambda\in 1\ldots n_{\rm pol}\f$ and each multiplier or divider variable \f$v_j, j\in 1\ldots n_{\rm var}\f$, respectively, are given by:
\f{align*}
\mathcal{M}_{\lambda,j} :=\ & \left\{ i~\middle| a_{\lambda, i} \neq 0 \wedge m_i({\bf v})\cdot v_j\notin \mathcal{S}_{\rm mon} \right\}\\
\mathcal{D}_{\lambda,j} :=\ & \left\{ i~\middle| a_{\lambda, i} \neq 0 \wedge m_i({\bf v})/ v_j\notin \mathcal{S}_{\rm mon} \right\}
\f}
This way, given \f$i\in\mathcal{M}_{\lambda,j}\f$ for some \f$j\in 1\ldots n_{\rm var}\f$ and \f$\lambda\in 1\ldots n_{\rm pol}\f$, the cut \f$h_{i,j} \geq m_{\lambda,j}\f$ accounts for the creation of the new monomial \f$m_i({\bf v})\cdot v_j\f$ on multiplying constraint \f$\lambda\f$ with variable \f$v_j\f$; likewise, given \f$i\in\mathcal{D}_{\lambda,j}\f$ for some \f$j\in 1\ldots n_{\rm var}\f$ and \f$\lambda\in 1\ldots n_{\rm pol}\f$, the cut \f$l_{i,j} \geq d_{\lambda,j}\f$ accounts for the creation of the new fractional term \f$m_i({\bf v})/ v_j\f$ on dividing constraint \f$\lambda\f$ with variable \f$v_j\f$.

- Finally, the factor \f$(1-\epsilon)\f$ multiplying the binary variables \f$b_{i,j}\f$ and \f$f_{i,j}\f$ in the objective function penalizes assignments that create as many new constraints as bilinear terms -- which would otherwise add up to zero in the objective function and yield a large number of (globally optimal) solutions for the ILP. The objective in doing so is therefore to determine those assignments yielding a minimal assignment. 
.

Notice that taking all of the binary variables \f$h_{i,j},l_{i,j},m_{\lambda,j},d_{\lambda,j}\f$ equal to 0 gives a feasible solution to the ILP with an objective of 0. Of interest for us, of course, are optimal solution with a strictly negative objective, so that exactly \f$M-N>0\f$ reduction constraints could be identified, with \f$M:=\sum_{j=1}^{n_{\rm var}} \sum_{\lambda=1}^{n_{\rm pol}} (d_{\lambda,j}+m_{\lambda,j})\f$ and \f$N:=\sum_{j=1}^{n_{\rm var}} \sum_{i=1}^{n_{\rm pol}} (h_{i,j}+l_{i,j})\f$.

In order to recover their symbolic expressions, we proceed by first writing the complete set of new polynomial constraints obtained from the multiplication or division of existing polynomial constraints with variables,
\f{align*}
\forall (\lambda,j): m_{\lambda,j}=1, \quad & \sum_{k=1}^{M_\lambda} a_{\lambda,k} \cdot m_{\lambda,k}({\bf v})\cdot v_j = 0\\
\forall (\lambda,j): d_{\lambda,j}=1, \quad & \sum_{k=1}^{M_\lambda} a_{\lambda,k} \cdot m_{\lambda,k}({\bf v})/ v_j = 0
\f}
These constraints can be rewritten in matrix form as
\f{align*}
  {\bf A}\cdot{\bf n}({\bf v}) = & {\bf b}({\bf v})
\f}
where:
- \f${\bf b}({\bf v})\f$ is an \f$M\f$-dimensional vector, with entries
\f{align*}
  {\bf b}({\bf v}) & := \left( \begin{array}{ll} \forall (\lambda,j): m_{\lambda,j}=1, & \sum_{i: m_i({\bf v})\cdot v_j\notin \mathcal{S}_{\rm mon}} a_{\lambda,i}\cdot m_i({\bf v})\cdot v_j \\ \forall (\lambda,j): d_{\lambda,j}=1, & \sum_{i: m_i({\bf v})/ v_j\notin \mathcal{S}_{\rm mon}} a_{\lambda,i}\cdot m_i({\bf v})/ v_j \end{array} \right)
\f}
- \f${\bf n}({\bf v})\f$ is an \f$N\f$-dimensional vector, comprised of the newly created monomial terms
\f{align*}
  {\bf n}({\bf v}) & := \left( \begin{array}{cl} \forall(i,j):h_{i,j}=1, & m({\bf v})\cdot v_j \\ \forall(i,j):l_{i,j}=1, &  m({\bf v})/ v_j \end{array} \right)
\f}
- \f${\bf A}\f$ is an \f$M\f$-by-\f$N\f$ real matrix
.
Next, we perform an LU factorization of the coefficient matrix \f${\bf A}\f$, such that \f${\bf P}\cdot{\bf A}={\bf L}\cdot{\bf U}\f$. The bottom \f$M-N>0\f$ rows of the vector \f${\bf p}({\bf v})\f$ such that \f${\bf L}\cdot{\bf r}({\bf v})={\bf P}\cdot{\bf b}({\bf v})\f$ -- as obtained for instance from a simple forward elimination -- are equal to zero by construction, and yield valid reduction polynomial constraints,
\f{align*}
  r_{N+1}({\bf v}) & = 0\\
  &\ \vdots\\
  r_{M}({\bf v}) & = 0
\f}
For large-scale problem, one may use efficient sparse implementations of LU factorization, based for instance on a bordered block decomposition of the original coefficient matrix -- see, e.g., available routines within the <A href="http://www.hsl.rl.ac.uk/catalogue/">Harwell Scientific Library (HSL)</A>.


\section sec_SRED_use How Do I Search for Reduction Constraints?

For illustration, consider the following pair of linear and quadratic constraints in three variables \f$x_0,x_1,x_2\f$:
\f{align*}
  0 = & x_0 + x_1 - 1\\ 
  0 = & x_0^2 - x_1^2 - x_2
\f}

The search for reduction constraints requires the header file <tt>sred.hpp</tt> to be included:

\code
#include "sred.hpp"
\endcode

Sparse polynomials are created for each constraint:

\code
const unsigned NX = 3, NF = 2;
mc::SPoly X[NX];
for( unsigned i=0; i<NX; i++ ) X[i].var( i );
mc::SPoly F[NF];
F[0] = X[0] + X[1] - 1.;
F[1] = sqr(X[0]) - sqr(X[1]) - X[2];
\endcode

Next, an environment <a>mc::SRed</a> is defined and the methods <a>mc::SRed::set_monomials</a> and <a>mc::SRed::search_reductions</a> are invoked to define the list of existing monomials and to search for new reduction constraints, respectively:

\code
mc::SRed SRF;
SRF.set_monomials( NF, F, &t_SPoly::mapmon );
SRF.search_reductions( NF, F, &t_SPoly::mapmon );
cout << SRF;
\endcode

The following information is displayed in this instance:

\verbatim
  2 Monomial List: [ [0]^2 [1]^2 ]

  3 Multiplier/Divider List: [ [0] [1] [2] ]

  1 Reduction Constraint:
   0 = 1.00000e+00 [0] - 1.00000e+00 [1] - 1.00000e+00 [0]^2 + 1.00000e+00 [1]^2

  1 Linearizable Monomial: [ [1]^2 ]
\endverbatim

These results show that the algorithm has identified the following reduction constraint:
\f{align*}
  0 =\ & x_0 - x_1 + x_0^2 - x_1^2
\f}
It is easy to check that this new constraint is obtained upon multiplying the existing linear constraint \f$0 = x_0 + x_1 - 1\f$ with either variable \f$x_0\f$ or \f$x_1\f$ then subtracting them to eliminate the cross-term \f$x_0x_1\f$. Notice also that the solution set of the new reduction constraint is \f$\left\{x_0,x_1\middle| x_0+x_1=1 \vee x_0-x_1=0 \right\}\f$, which is indeed redundant with the original constraint \f$ x_0+x_1=1\f$.


\section sec_SRED_opt What are the options in mc::SRed and how do I set them?

The public static class member mc::SRed::options that can be used to set/modify the options; e.g.,

\code
      mc::SQuad::options.ORDER        = 2;
      mc::SQuad::options.MIPDISPLEVEL = 1;    
\endcode

The available options are the following:

<TABLE border="1">
 <TR><TH><b>Name</b>  <TD><b>Type</b> <TD><b>Default</b> <TD><b>Description</b>
 <TR><TD><tt>mc::SRed::Options::ORDER</tt> <TD><tt>unsigned int</tt> <TD>1 <TD>Order of reduction, i.e. maximal degree of multiplier/divider monomials
 <TR><TD><tt>mc::SRed::Options::NODIV</tt> <TD><tt>bool</tt> <TD>true <TD>Whether to disable divider variables in addition to multiplier variables
 <TR><TD><tt>mc::SRed::Options::REDTOL</tt> <TD><tt>double</tt> <TD>1e-3 <TD>Tolerance on new monomials in objective function to penalize operations that do not create redundancy
 <TR><TD><tt>mc::SRed::Options::ZEROTOL</tt> <TD><tt>double</tt> <TD>mc::machprec <TD>Tolerance on polynomial coefficients in reduction constraints
 <TR><TD><tt>mc::SRed::Options::RKATOL</tt> <TD><tt>double</tt> <TD>mc::machprec <TD>Absolute tolerance for rank computation using SVD
 <TR><TD><tt>mc::SRed::Options::RKRTOL</tt> <TD><tt>double</tt> <TD>mc::machprec <TD>Relative tolerance for rank computation using SVD
 <TR><TD><tt>mc::SRed::Options::LPALGO</tt> <TD><tt>int</tt> <TD>-1 <TD>LP algorithm used by MIP solver
 <TR><TD><tt>mc::SRed::Options::LPPRESOLVE</tt> <TD><tt>int</tt> <TD>-1 <TD>LP presolve strategy in MIP solver
 <TR><TD><tt>mc::SRed::Options::LPFEASTOL</tt> <TD><tt>double</tt> <TD>1e-7 <TD>Tolerance on LP feasibility in MIP solver
 <TR><TD><tt>mc::SRed::Options::LPOPTIMTOL</tt> <TD><tt>double</tt> <TD>1e-7 <TD>Tolerance on LP optimality in MIP solver
 <TR><TD><tt>mc::SRed::Options::MIPRELGAP</tt> <TD><tt>double</tt> <TD>1e-3 <TD>Tolerance on relative gap in MIP solver
 <TR><TD><tt>mc::SRed::Options::MIPABSGAP</tt> <TD><tt>double</tt> <TD>1e-3 <TD>Tolerance on absolute gap in MIP solver
 <TR><TD><tt>mc::SRed::Options::MIPTHREADS</tt> <TD><tt>int</tt> <TD>0 <TD>Number of threads used by MIP solver (default value of 0 enables all available threads)
 <TR><TD><tt>mc::SRed::Options::MIPCONCURRENT</tt> <TD><tt>int</tt> <TD>1 <TD>Number of independent MIP solves in parallel
 <TR><TD><tt>mc::SRed::Options::MIPFOCUS</tt> <TD><tt>int</tt> <TD>0 <TD>MIP high-level solution strategy (default value of 0 seeks to a balance between finding new feasible solutions and proving that the current solution is optimal)
 <TR><TD><tt>mc::SRed::Options::MIPHEURISTICS</tt> <TD><tt>double</tt> <TD>0.2 <TD>Fraction of time spent in MIP heuristics (default value of 0.2 aims to spend 20% of runtime on heuristics)
 <TR><TD><tt>mc::SRed::Options::MIPNUMFOCUS</tt> <TD><tt>int</tt> <TD>0 <TD>Degree to which the code attempts to detect and manage numerical issues (default setting of 0 makes an automatic choice, with a slight preference for speed)
 <TR><TD><tt>mc::SRed::Options::MIPDISPLEVEL</tt> <TD><tt>int</tt> <TD>0 <TD>Display level for MIP solver
 <TR><TD><tt>mc::SRed::Options::MIPOUTPUTFILE</tt> <TD><tt>std::string</tt> <TD>"" <TD>Name of output file for MIP model
 <TR><TD><tt>mc::SRed::Options::MIPTIMELIMIT</tt> <TD><tt>double</tt> <TD>60 <TD>Maximum MIP runtime (seconds)
 <TR><TD><tt>mc::SRed::Options::DISPLEVEL</tt> <TD><tt>int</tt> <TD>0 <TD>Display level for reduction procedure
 <TR><TD><tt>mc::SRed::Options::DISPLEN</tt> <TD><tt>unsigned int</tt> <TD>5 <TD>Number of digits in output stream
</TABLE>


\section sec_SRED_err What Errors Can Be Encountered during the Quadratization of a Multivariate Polynomial?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type mc::SQuad::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during a quadratization, and then make the appropriate changes. Should an exception be thrown and not caught, the program will stop.

Possible errors encountered during quadratization of a multivariate polynomial are:

<TABLE border="1">
 <TR><TH><b>Number</b> <TD><b>Description</b>
 <TR><TH><tt>1</tt> <TD>Maximal problem size reached
 <TR><TH><tt>2</tt> <TD>Mixed-integer programming solver disabled
 <TR><TH><tt>3</tt> <TD>LU decomposition failed
 <TR><TH><tt>4</tt> <TD>SV decomposition failed
 <TR><TH><tt>-33</tt> <TD>Internal error
</TABLE>


\section sec_SRED_refs References

- L Liberti, <A href="https://doi.org/10.1111/j.1475-3995.2004.00438.x">Reduction constraints for the global optimization of NLPs</A>, <i>International Transactions in Operation Research</i>, <b>11</b>(1):33-41, 2004.
- L Liberti, CC Pantelides, <A href="https://doi.org/10.1007/s10898-006-9005-4">An exact reformulation algorithm for large nonconvex NLPs involving bilinear terms</A>, <I>Journal of Global Optimization</I> <B>36</B>:161–189, 2006
- H Sherali, A Alameddine, <A href="https://doi.org/10.1007/BF00122429">A new reformulation-linearization technique for
bilinear programming problems</A>, <i>Journal of Global Optimization</i>, <b>2</b>:379-410.
- JP Ruiz, IE Grossmann, <A href="https://doi.org/10.1016/j.compchemeng.2011.01.035">Using redundancy to strengthen the relaxation for the global optimization of MINLP problems</A>, <I>Computers & Chemical Engineering</I> <B>35</B>:2729–2740, 2011.
.
*/

#ifndef MC__SRED_H
#define MC__SRED_H

#include <list>
#include <tuple>
#include <climits>
#include "spoly.hpp"
#include "mclapack.hpp"

#if defined(MC__USE_GUROBI)
 #include "gurobi_c++.h"
 extern "C"{
  #include <fenv.h>
  int fedisableexcept( int );
 }
#endif

#define MC__SRED_CHECK
#undef  MC__SRED_PROCESS_DEBUG
#undef  MC__SRED_CHECK_START

namespace mc
{

//! @brief C++ structure for ordering of contraint-monomial pair
template <typename COMP=std::less<unsigned>>
struct lt_CtrMon
{
  // Comparison operator
  template <typename KEY>
  bool operator
    ()
    ( std::pair< unsigned, SMon<KEY,COMP> const* > const& iMon,
      std::pair< unsigned, SMon<KEY,COMP> const* > const& jMon )
    const
    {
      // Order based on constraint index first
      if( iMon.first < jMon.first ) return true;
      if( iMon.first > jMon.first ) return false;
      // Order based on multiplier monomial next
      if( lt_SMon<COMP>()( *iMon.second, *jMon.second ) ) return true;
      if( lt_SMon<COMP>()( *jMon.second, *iMon.second ) ) return false;
      // Pairs are identical on reaching this point
      return false;
    }
};

//! @brief C++ class for reformulation of sparse polynomial models as a set of quadratic forms
////////////////////////////////////////////////////////////////////////
//! mc::SQuad is a C++ class for reformulation of sparse polynomial
//! models as a set of quadratic forms via the introduction of
//! auxiliary variables (lifting). Sparse monomial are defined using the
//! class mc::SPolyModel. 
////////////////////////////////////////////////////////////////////////
template <typename KEY=unsigned, typename COMP=std::less<unsigned>>
class SRed
////////////////////////////////////////////////////////////////////////
{
  template <typename K, typename C> friend std::ostream& operator<< ( std::ostream&, SRed<K,C> const& );
 
public:

  typedef unsigned long long t_size;
  typedef SMon<KEY,COMP> t_SMon;
  typedef SPoly<KEY,COMP> t_SPoly;
  typedef SMonExt<KEY,COMP> t_SMonExt;
  typedef std::set< KEY, COMP > set_SVar;
  typedef std::set< t_SMon, lt_SMon<COMP> > set_SMon;
  typedef std::set< t_SMonExt, lt_SMonExt<COMP> > set_SMonExt;
  typedef std::set< t_SMonExt const*, lt_pSMonExt<COMP> > set_pSMonExt;
  typedef std::map< t_SMon, double, lt_SMon<COMP> > map_SPoly;
  typedef std::map< std::pair< unsigned, t_SMon const* >, set_pSMonExt, lt_CtrMon<COMP> > map_CtrMon;

  //! @brief Options of mc::SRed
  struct Options
  {
    //! @brief Constructor
    Options():
      ORDER(1), NODIV(true), REDTOL(1e-3), ZEROTOL(machprec()),
      RKATOL(machprec()), RKRTOL(machprec()),
#if defined(MC__USE_GUROBI)
      LPALGO( LPALGO_DEFAULT ), LPPRESOLVE(-1),
      LPFEASTOL(1e-7), LPOPTIMTOL(1e-7),
      MIPRELGAP(1e-3), MIPABSGAP(1e-3), MIPTHREADS(0),
      MIPCONCURRENT(1), MIPFOCUS(0), MIPHEURISTICS(0.2),
      MIPNUMFOCUS(0), MIPDISPLEVEL(0), MIPOUTPUTFILE(""),
      MIPTIMELIMIT(60),
#endif
      DISPLEVEL(0), DISPLEN(5)
      {}
    //! @brief Assignment operator
    Options& operator= ( Options const& opt ){
        ORDER           = opt.ORDER;
        NODIV           = opt.NODIV;
        REDTOL          = opt.REDTOL;
        ZEROTOL         = opt.ZEROTOL;
        RKATOL          = opt.RKATOL;
        RKRTOL          = opt.RKRTOL;
#if defined(MC__USE_GUROBI)
        LPALGO          = opt.LPALGO;
        LPPRESOLVE      = opt.LPPRESOLVE;
        LPFEASTOL       = opt.LPFEASTOL;
        LPOPTIMTOL      = opt.LPOPTIMTOL;
        MIPRELGAP       = opt.MIPRELGAP;
        MIPABSGAP       = opt.MIPABSGAP;
        MIPTHREADS      = opt.MIPTHREADS;
        MIPCONCURRENT   = opt.MIPCONCURRENT;
        MIPFOCUS        = opt.MIPFOCUS;
        MIPHEURISTICS   = opt.MIPHEURISTICS;
        MIPNUMFOCUS     = opt.MIPNUMFOCUS;
        MIPDISPLEVEL    = opt.MIPDISPLEVEL;
        MIPOUTPUTFILE   = opt.MIPOUTPUTFILE;
        MIPTIMELIMIT    = opt.MIPTIMELIMIT;
#endif
        DISPLEVEL       = opt.DISPLEVEL;
        DISPLEN         = opt.DISPLEN;
        return *this ;
      }
    //! @brief Reduction order
    unsigned ORDER;
    //! @brief Whether to disable divider variables in addition to multiplier variables
    bool NODIV;
    //! @brief Tolerance on variable to penalize operations that do not create redundancy (default: 1e-3)
    double REDTOL;
    //! @brief Tolerance on polynomial coefficients in reduction constraints (default: mc::machprec())
    double ZEROTOL;
    //! @brief Absolute tolerance for rank computation using SVD (default: mc::machprec())
    double RKATOL;
    //! @brief Relative tolerance for rank computation using SVD (default: mc::machprec())
    double RKRTOL;
#if defined(MC__USE_GUROBI)
    //! @brief LP algorithm used by MIP solver
    int LPALGO;
    //! @brief LP presolve strategy in MIP solver
    int LPPRESOLVE;
    //! @brief Tolerance on LP feasibility in MIP solver
    double LPFEASTOL;
     //! @brief Tolerance on LP optimality in MIP solver
    double LPOPTIMTOL;
    //! @brief Tolerance on relative gap in MIP solver
    double MIPRELGAP;
    //! @brief Tolerance on absolute gap in MIP solver
    double MIPABSGAP;
    //! @brief Number of threads used by MIP solver - default value of 0 allows to use all available threads
    int MIPTHREADS;
    //! @brief Number of independent MIP solves in parallel
    int MIPCONCURRENT;
    //! @brief MIP high-level solution strategy - default value of 0 seeks to a balance between finding new feasible solutions and proving that the current solution is optimal
    int MIPFOCUS;
    //! @brief Fraction of time spent in MIP heuristics - default value of 0.2 aims to spend 20% of runtime on heuristics
    double MIPHEURISTICS;
    //! @brief Degree to which the code attempts to detect and manage numerical issues - default setting of 0 makes an automatic choice, with a slight preference for speed 
    int MIPNUMFOCUS;
    //! @brief Display level for MIP solver
    int MIPDISPLEVEL;
    //! @brief Name of output file for MIP model
    std::string MIPOUTPUTFILE;
    //! @brief Maximum MIP run time (seconds)
    double MIPTIMELIMIT;
    //! @brief Default option for LP solver
    static const int LPALGO_DEFAULT = -1;
#endif
    //! @brief Display level for reduction procedure
    int DISPLEVEL;
    //! @brief Number of digits in output stream for sparse polynomial coefficients
    unsigned DISPLEN;
  } options;
    
  //! @brief Exceptions of mc::SRed
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SQuad exception handling
    enum TYPE{
      MAXSIZE  = 1,   //!< Maximal problem size reached
      MIPERR,         //!< Call to MIP solver disabled
      LUERR,          //!< LU decomposition failed
      SVDERR,         //!< SV decomposition of coefficient matrix failed
      INTERNAL = -33  //!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
       case MAXSIZE:
        return "mc::SRed\t Maximal problem size reached";
      case MIPERR:
        return "mc::SRed\t Mixed-integer programming solver disabled";
      case LUERR:
        return "mc::SRed\t LU decomposition failed";
      case SVDERR:
        return "mc::SRed\t SV decomposition failed";
       case INTERNAL:
       default:
        return "mc::SRed\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

protected:

  //! @brief Set of existing monomials in reduction problem
  set_SMon _SetMon;

  //! @brief Set of multiplier/divider monomials in reduction problem
  set_SMon _SetMonOp;

  //! @brief Set of new candidate monomials in reduction problem
  set_SMonExt _SetMonNew;

  //! @brief Set of particpating variables in reduction problem
  set_SVar _SetVar;

  //! @brief Map of new monomials created upon multiplying a polynomial constraint with a multiplier monomial
  map_CtrMon _NewMul;

  //! @brief Map of new monomials created upon dividing a polynomial constraint with a divider monomial
  map_CtrMon _NewDiv;

  //! @brief Set of original polynomial constraints participating in the reduction
  std::set<unsigned> _ndxCtrUse;

#if defined(MC__USE_GUROBI)
  //! @brief whether the MIP solver has sent an exception
  bool _MIPexcpt;
  
  //! @brief Set of monomials in MIP optimization model
  set_SMon _MIP_SetMon;
  
  //! @brief Dummy monomial storing partipating variables and their higest degrees in reduction problem
  t_SMon _VarDeg;
  
  //! @brief MIP environment
  GRBEnv* _GRBenv;
  
  //! @brief MIP model
  GRBModel* _GRBmodel;
  
  //! @brief map of binary variables for candidate multiplier monomials for each constraint
  std::map< std::pair< unsigned, t_SMon const* >, GRBVar, lt_CtrMon<COMP> > _MIP_CtrMul;
  
  //! @brief map of binary variables for candidate divider monomials for each constraint
  std::map< std::pair< unsigned, t_SMon const* >, GRBVar, lt_CtrMon<COMP> > _MIP_CtrDiv;
  
  //! @brief map of binary variables for new monomials from multiplication
  std::map< t_SMonExt const*, GRBVar, lt_pSMonExt<COMP> > _MIP_NewMul;
  
  //! @brief set of candidate divider monomials for each constraint
  std::map< t_SMonExt const*, GRBVar, lt_pSMonExt<COMP> > _MIP_NewDiv;
#endif

  //! @brief dense oefficient matrix for constraint reduction
  CPPL::dgematrix _mat;

  //! @brief sparse coefficient matrix for constraint reduction
  std::vector< std::map< t_SMonExt const*, double, lt_pSMonExt<COMP> > > _coef;
  
  //! @brief coefficient matrix for constraint reduction
  std::map< t_SMonExt const*, unsigned, lt_pSMonExt<COMP> > _new;

  //! @brief right-hand side vector for constraint reduction
  std::vector< t_SPoly > _rhs;

  //! @brief Reduction polynomial constraints
  std::vector< t_SPoly > _RedCtr;

  //! @brief sparse coefficient matrix for monomial elimination
  std::vector< std::map< t_SMon const*, double, lt_pSMon<COMP> > > _lhs;

  //! @brief Map of candidate monomials to variables for elimination
  std::map< t_SMon const*, unsigned, lt_pSMon<COMP> > _mapMonVar;

  //! @brief Map of variables for elimination to candidate monomials
  std::map< unsigned, t_SMon const* > _mapVarMon;

  //! @brief Map of variables with associated weights/priorities
  std::map<unsigned,double> _mapVarWeight;

  //! @brief Map of weights/priorities with associated variables
  std::multimap<double,unsigned> _mapWeightVar;

  //! @brief Set of eliminated monomials through exact linearization
  set_SMon _setMonElim;

#if defined(MC__USE_GUROBI)
  //! @brief Optimize for maximal reduction constraints 
  void _MIP_optimize
    ( std::set<unsigned> const& ndxSPol );

  //! @brief Solve MIP optimization model for reduction constraint search
  void _MIP_solve
    ();

  //! @brief Encode MIP optimization model for reduction constraint search
  void _MIP_encode
    ( std::set<unsigned> const& ndxSPol );

  //! @brief Display MIP optimal solution for reduction constraint search
  void _MIP_display
    ();

  //! @brief Reset MIP optimization model for reduction constraint search
  void _MIP_reset
    ();

  //! @brief Set options in MIP optimization model for reduction constraint search
  void _MIP_options
    ();
#endif

  //! @brief Reset the reduction constraints
  void _Reduc_reset
    ();

  //! @brief Reset the elimination variables
  void _Elim_reset
    ();

  //! @brief Search for exact linearization for monomials participating in the sparse polynomials in array <a>pSPol</a>
  template <typename POL>
  unsigned _search_linearizations
    ( std::set<unsigned> const& ndxSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const );

public:

  //! @brief Default Constructor
  SRed
    ()
    : _nord(0), _nvar(0), _nmon(0),
      _posord(nullptr), _expmon(nullptr), _binom(nullptr)
    {
#if defined(MC__USE_GUROBI)
      _GRBenv   = new GRBEnv();
      _GRBmodel = nullptr;
#endif
    }

  //! @brief Destructor
  virtual ~SRed
    ()
    {
      _cleanmon();
#if defined(MC__USE_GUROBI)
      delete _GRBmodel;
      delete _GRBenv;
#endif
    }

  //! @brief Process the sparse polynomials in array <a>pSPol</a> indexed by <a>ndxSPol</a> to extract their participating monomials
  template <typename POL>
  void set_monomials
    ( std::set<unsigned> const& ndxSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
      bool const APPEND=false );

  //! @brief Process the <a>nSPol</a> sparse polynomials in array <a>pSPol</a> to extract their participating monomials
  template <typename POL>
  void set_monomials
    ( unsigned const nSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
      bool const APPEND=false );

  //! @brief Process the sparse polynomials <a>SPol</a> to extract their participating monomials
  template <typename POL>
  void set_monomials
    ( POL const& SPol, map_SPoly const& (POL::*mapmon)() const,
      bool const APPEND=false );

  //! @brief Search for reduction constraints based on the sparse polynomials in array <a>pSPol</a> indexed by <a>ndxSPol</a>
  template <typename POL>
  unsigned search_reductions
    ( std::set<unsigned> const& ndxSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
      bool const LINEARIZE=true );

  //! @brief Search for reduction constraints based on the <a>nSPol</a >sparse polynomials in array <a>pSPol</a>
  template <typename POL>
  unsigned search_reductions
    ( unsigned const nSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
      bool const LINEARIZE=true );

  //! @brief Return reduction polynomial constraints
  std::vector<t_SPoly> const& RedCtr
    ()
    const
    { return _RedCtr; }

  //! @brief Return set of eliminated monomials through exact linearization
  set_SMon const& ElimMon
    ()
    { return _setMonElim; }

private:

  //! @brief Order of reduction
  unsigned _nord;

  //! @brief Number of variables
  unsigned _nvar;

  //! @brief Total number of monomial candidates
  t_size _nmon;

  //! @brief Positions of first monomial term of order <tt>iord=1,...,_nord</tt>
  t_size *_posord;

  //! @brief Variable exponents of monomial terms. The exponent for variable <tt>ivar</tt> in monomial term <tt>imon</tt> is at location <tt>imon*_nvar+ivar</tt>.
  unsigned *_expmon;

  //! @brief Array of <tt>(_nvar+_nord-1)*(_nord+1)</tt> contining binomial coefficients
  t_size *_binom;

  //! @brief Maximum binomial coefficients in array _binom
  std::pair<unsigned, unsigned> _binom_size;

  //! @brief Set reduction order and internal arrays
  void _sizemon
    ();

  //! @brief Clean up internal arrays
  void _cleanmon
    ();

  //! @brief Populate array <tt>_posord</tt> up to order <tt>nord</tt>
  void _set_posord
    ();

  //! @brief Populate array <tt>_expmon</tt> up to order <tt>nord</tt>
  void _set_expmon
    ();

  //! @brief Generate variable exponents <tt>iexp</tt> for subsequent monomial order <tt>iord</tt>
  void _next_expmon
    ( unsigned *iexp, const unsigned iord )
    const;

  //! @brief Get index of monomial term with variable exponents <tt>iexp</tt> in <tt>1,...,_nmon</tt>
  t_size _loc_expmon
    ( const unsigned *iexp )
    const;
    
  //! @brief Populate array <tt>_binom</tt> with binomial coefficients up to order <tt>nord</tt>
  void _set_binom
    ();

  //! @brief Get binomial coefficient \f$\left(\stackrel{n}{k}\right)\f$
  t_size _get_binom
    ( const unsigned n, const unsigned k )
    const;

  //! @brief Get matrix rank from diagonal array
  unsigned _rank
    ( unsigned const dim, double const* arr )
    const;
};

//template <typename KEY, typename COMP>
//inline typename SRed<KEY,COMP>::Options SRed<KEY,COMP>::options;

#if defined(MC__USE_GUROBI)
template <typename KEY, typename COMP>
inline int const SRed<KEY,COMP>::Options::LPALGO_DEFAULT;
#endif

////////////////////////////////////////////////////////////////////////


template <typename KEY, typename COMP>
template <typename POL >
inline 
void
SRed<KEY,COMP>::set_monomials
( std::set<unsigned> const& ndxSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
  bool const APPEND )
{
  if( !APPEND ) _SetMon.clear();
  for( unsigned const& i : ndxSPol )
    set_monomials( pSPol[i], mapmon, true );
}

template <typename KEY, typename COMP>
template <typename POL >
inline 
void
SRed<KEY,COMP>::set_monomials
( unsigned const nSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
  bool const APPEND )
{
  if( !APPEND ) _SetMon.clear();
  for( unsigned i=0; i<nSPol; ++i )
    set_monomials( pSPol[i], mapmon, true );
}

template <typename KEY, typename COMP>
template <typename POL >
inline 
void
SRed<KEY,COMP>::set_monomials
( POL const& SPol, map_SPoly const& (POL::*mapmon)() const,
  bool const APPEND )
{
  if( !APPEND ) _SetMon.clear();
  for( auto const& [mon,coef] : (SPol.*mapmon)() )
    if( mon.tord > 1 ) _SetMon.insert( mon );
}

template <typename KEY, typename COMP>
template <typename POL >
inline 
unsigned
SRed<KEY,COMP>::search_reductions
( unsigned const nSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
  bool const LINEARIZE )
{
  std::set<unsigned> ndxSPol;
  for( unsigned i=0; i<nSPol; ++i ) ndxSPol.insert(i);
  return search_reductions( ndxSPol, pSPol, mapmon, LINEARIZE );
}

template <typename KEY, typename COMP>
template <typename POL >
inline 
unsigned
SRed<KEY,COMP>::search_reductions
( std::set<unsigned> const& ndxSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
  bool const LINEARIZE )
{
  _Reduc_reset();
  if( ndxSPol.empty() ) return 0;
#ifdef MC__SRED_DEBUG
  for( auto const& i : ndxSPol )
    std::cout << "pSPol[" << i << "]: " << pSPol[i] << std::endl;
#endif
  
  // Keep track of participating variables in reduction problem
  _VarDeg.expr.clear();
  _VarDeg.tord = 0;
  for( auto const& mon : _SetMon )
    _VarDeg.hull( mon ); 
  for( auto const& i : ndxSPol )
    for( auto const& [mon,coef] : (pSPol[i].*mapmon)() ){
      if( coef == 0. ) continue;
      _VarDeg.hull( mon );
    }
  _SetVar.clear();
  for( auto const& [var,ord] : _VarDeg.expr )
    _SetVar.insert( var );

  // Generate set of multiplier/divider monomials
  _sizemon();

  // Create sets of new monomials for each polynomial constraint and each multiplier/divider monomial
  for( t_size k=0; k<_nmon; ++k ){
    // Create candidate multiplier/divider monomial
    t_SMon kmon;
    unsigned j = 0;
    for( auto const& var : _SetVar ){
      unsigned const& exp = _expmon[ k*_nvar+(j++) ];
      if( !exp ) continue;
      kmon.expr[var] = exp;
      kmon.tord += exp;
    }
    if( !kmon.tord ) continue;
    auto [ktmon,ins] = _SetMonOp.insert( kmon );
#ifdef MC__SRED_DEBUG
    std::cout << "kmon: " << kmon.display(0) << std::endl;
#endif
    // Test if multiplier monomial exists in product with monomials in constraints
    for( auto const& i : ndxSPol ){
      for( auto const& [mon,coef] : (pSPol[i].*mapmon)() ){
        if( coef == 0. ) continue;
        auto jtmon = (mon.tord>1? _SetMon.find( mon ): _SetMon.begin());
        if( jtmon == _SetMon.end() ) continue;
        // New multiplier term added if not present and has order greater than one
        t_SMon monnew( kmon );
        if( mon.tord ) monnew += mon;
        if( monnew.tord > 1 && !_SetMon.count( monnew ) ){
          auto [iktmon,ins] = _SetMonNew.insert( monnew );
          _NewMul[std::make_pair(i,&*ktmon)].insert( &*iktmon );
#ifdef MC__SRED_DEBUG
          std::cout << "  ctr " << i << ", mon " << jtmon->display(0) << std::endl;
#endif
        }
        if( options.NODIV ) continue;
        // New division term added if not present and has order greater than one
        if( !kmon.subseteq( mon ) || ( mon.tord > kmon.tord+1 && !_SetMon.count( mon - kmon ) ) ){
#ifdef MC__SRED_DEBUG
          std::cout << "  imon " << mon.display(0) << " - kmon " << kmon.display(0) << std::endl;
#endif
          t_SMonExt monnew( mon ); monnew -= kmon;
          auto [iktmon,ins] = _SetMonNew.insert( monnew );
          _NewDiv[std::make_pair(i,&*ktmon)].insert( &*iktmon );
        }
      }
    }
  }
  
  // Search reduction constraints
#if defined(MC__USE_GUROBI)
  _MIP_optimize( ndxSPol );
#else
  throw Exceptions( Exceptions::MIPERR );
#endif

  // Coefficient matrix and right-hand-side vector holding the reduction constraints
  for( auto const& [keyctr,grbvar] : _MIP_CtrMul ){
    if( grbvar.get(GRB_DoubleAttr_X) < 0.9 ) continue;
    auto const& [ctr,monop] = keyctr;
    _ndxCtrUse.insert( ctr );
    _rhs.push_back( 0. );
    _coef.push_back( std::map<t_SMonExt const*,double,lt_pSMonExt<COMP>>() );
    for( auto const& [mon,coef] : (pSPol[ctr].*mapmon)() ){
      if( coef == 0. ) continue;
      t_SMonExt monnew( mon ); monnew += *monop;
      auto itmonnew = _SetMonNew.find( monnew );
      if( itmonnew == _SetMonNew.end() )
        _rhs.back() += std::make_pair( mon + *monop, coef );
      else{
#ifdef MC__SRED_CHECK
        auto itmip = _MIP_NewMul.find( &monnew );
#ifdef MC__SRED_DEBUG
        std::cout << "ctr:" << ctr << " monop:" << monop->display(0) << " monnew:" << monnew.display(0)
                  << " " << itmip->second.get(GRB_DoubleAttr_X) << std::endl; 
#endif
        assert( itmip != _MIP_NewMul.end() && itmip->second.get(GRB_DoubleAttr_X) > 0.9 );
#endif
        _new.insert( std::make_pair( &*itmonnew, _new.size() ) );
        _coef.back().insert( std::make_pair( &*itmonnew, coef ) );
      }
    }
  }

  for( auto const& [keyctr,grbvar] : _MIP_CtrDiv ){
    if( grbvar.get(GRB_DoubleAttr_X) < 0.9 ) continue;
    auto const& [ctr,monop] = keyctr;
    _ndxCtrUse.insert( ctr );
    _rhs.push_back( 0. );
    _coef.push_back( std::map<t_SMonExt const*,double,lt_pSMonExt<COMP>>() );
    for( auto const& [mon,coef] : (pSPol[ctr].*mapmon)() ){
      if( coef == 0. ) continue;
      t_SMonExt monnew( mon ); monnew -= *monop;
      auto itmonnew = _SetMonNew.find( monnew );
      if( itmonnew == _SetMonNew.end() )
        _rhs.back() += std::make_pair( mon - *monop, coef );
      else{
#ifdef MC__SRED_CHECK
        auto itmip = _MIP_NewDiv.find( &monnew );
        assert( itmip != _MIP_NewDiv.end() && itmip->second.get(GRB_DoubleAttr_X) > 0.9 );
#endif
        _new.insert( std::make_pair( &*itmonnew, _new.size() ) );
        _coef.back().insert( std::make_pair( &*itmonnew, coef ) );
      }
    }
  }

#ifdef MC__SRED_DEBUG
  unsigned irhs = 0;
  for( auto const& pol : _rhs )
    std::cout << "_rhs[" << irhs++ << "]: " << pol << std::endl;
  for( auto const& [monnew,pos] : _new )
    std::cout << "_new[" << monnew->display(0) << "]: " << pos << std::endl;
#endif

  // Perform LU decomposition of coefficient matrix
  _mat.resize( _rhs.size(), _new.size() );
  _mat.zero();
  for( unsigned irow=0; irow<_coef.size(); ++irow )
    for( auto const& [monnew,coef] : _coef[irow] )
      _mat( irow, _new[monnew] ) = coef;

  CPPL::dgematrix L, U;
  std::vector<int> PIV;
  int INFO = CPPL::dgetrf( _mat, U, L, PIV );
  if( INFO < 0 ) throw Exceptions( Exceptions::LUERR );
#ifdef MC__SRED_DEBUG
  std::cout << "\nLU DECOMPOSITION OF COEFFICIENT MATRIX: "
            << INFO << std::endl
            << "\nMatrix U:\n" << U
            << "\nMatrix L:\n" << L
            << "\nVector PIV:\n";
  for( auto&& ip : PIV ) std::cout << ip << std::endl;
#endif

  // Recover reduction constraints 
  for( unsigned irow=0; irow<_new.size(); ++irow ){
    for( unsigned krow=0; krow<irow; ++krow )
      _rhs[PIV[irow]] -= L(irow,krow) * _rhs[PIV[krow]];
    _rhs[PIV[irow]] /= L(irow,irow);
  }
  for( unsigned irow=_new.size(); irow<_rhs.size(); ++irow ){
    for( unsigned krow=0; krow<_new.size(); ++krow )
      _rhs[PIV[irow]] -= L(irow,krow) * _rhs[PIV[krow]];
  }
  _RedCtr.clear();
  for( unsigned irow=0; irow<_rhs.size()-_new.size(); ++irow ){
    _rhs[PIV[_new.size()+irow]].clean( options.ZEROTOL );
    if( _rhs[PIV[_new.size()+irow]].mapmon().empty() ) continue;
    _RedCtr.push_back( _rhs[PIV[_new.size()+irow]] );
  }
  
  // Search linearization variables
  return( LINEARIZE? _search_linearizations( ndxSPol, pSPol, mapmon ): _RedCtr.size() );
}

template <typename KEY, typename COMP>
template <typename POL >
inline 
unsigned
SRed<KEY,COMP>::_search_linearizations
( std::set<unsigned> const& ndxSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const )
{
  _Elim_reset();
  //if( _ndxCtrUse.empty() && _RedCtr.empty() ) return 0;
  if( _RedCtr.empty() ) return 0;

  // Coefficient matrix and right-hand-side vector holding the linearized constraints
  // Process original polynomial constraints
  unsigned nvar = 0, nrhs = 0;
/*
  for( auto const& ictr : ndxSPol ){
  //for( auto const& ictr : _ndxCtrUse ){
    if( (pSPol[ictr].*mapmon)().rbegin()->first.tord < 2 ) continue;
    nrhs++;
    _lhs.push_back( std::map<t_SMon const*,double,lt_pSMon<COMP>>() );
    for( auto const& [mon,coef] : (pSPol[ictr].*mapmon)() ){
      if( coef == 0. || mon.tord < 2 )
        continue; // exclude zero, constant and linear terms from elimination
      auto const& [itmon,ins] = _mapMonVar.insert( std::make_pair( &mon, nvar ) );
      if( ins ){
        _mapVarMon.insert( std::make_pair( nvar, &mon ) );
        _mapVarWeight[nvar] = 1;
        nvar++;
      }
      else{
        _mapVarWeight[itmon->second]++;
      }
      _lhs.back().insert( std::make_pair( itmon->first, coef ) );
    }
  }
*/
  // Process reduction polynomial constraints
  for( unsigned ired=0; ired<_RedCtr.size(); ++ired ){
    if( _RedCtr[ired].mapmon().rbegin()->first.tord < 2 ) continue;
    nrhs++;
    _lhs.push_back( std::map<t_SMon const*,double,lt_pSMon<COMP>>() );
    for( auto const& [mon,coef] : _RedCtr[ired].mapmon() ){
      if( coef == 0. || mon.tord < 2 )
        continue; // exclude zero, constant and linear terms from elimination
      auto const& [itmon,ins] = _mapMonVar.insert( std::make_pair( &mon, nvar ) );
      if( ins ){
        _mapVarMon.insert( std::make_pair( nvar, &mon ) );
        _mapVarWeight[nvar] = 1;
        nvar++;
      }
      else{
        _mapVarWeight[itmon->second]++;
      }
      _lhs.back().insert( std::make_pair( itmon->first, coef ) );
    }
  }
#ifdef MC__SRED_DEBUG
  for( auto const& [mon,ndx] : _mapMonVar )
    std::cout << "_mapMonVar[" << mon->display(0) << "]: " << ndx << std::endl;
#endif

  // Terminate if size zero
  if( !nrhs || !nvar ) return 0;

  // Map of variable weights
  for( auto const& [var,wei] : _mapVarWeight ){
    _mapWeightVar.insert( std::make_pair( wei, var ) );
#ifdef MC__SRED_DEBUG
    std::cout << "_mapWeightVar[" << var << "]: " << wei << std::endl;
#endif
  }

  // Perform LU decomposition of coefficient matrix
  _mat.resize( nrhs, nvar );
  _mat.zero();
  for( unsigned irow=0; irow<_lhs.size(); ++irow )
    for( auto const& [mon,coef] : _lhs[irow] )
      _mat( irow, _mapMonVar[mon] ) = coef;
#ifdef MC__SRED_DEBUG
  std::cout << "\nMatrix A = \n" << _mat;
#endif

  CPPL::dgematrix cpmat( _mat );
  CPPL::dcovector S;
  CPPL::dgematrix U, VT;
  int INFO = cpmat.dgesvd( S, U, VT );
  if( INFO < 0 ) throw Exceptions( Exceptions::SVDERR );
#ifdef MC__SRED_DEBUG
  std::cout << "\nSV DECOMPOSITION OF COEFFICIENT MATRIX: "
            << INFO << std::endl
            << "\nMatrix U:\n" << U
            << "\nVector S:\n" << S
            << "\nMatrix VT:\n" << VT;
#endif

  unsigned rkmat = _rank( S.l, S.array );
#ifdef MC__SRED_DEBUG
  std::cout << "\n Rank: " << rkmat << std::endl;
#endif

  // Can all monomials be eliminated?
  if( rkmat == _mapMonVar.size() ){
    for( auto const& [mon,ndx] : _mapMonVar )
      _setMonElim.insert( *mon );
    return rkmat;
  }

  // Else identify maximal elimination
  CPPL::dgematrix U1T( rkmat, nrhs );
  for( unsigned i=0; i<nrhs; ++i )
    for( unsigned j=0; j<rkmat; ++j )
      U1T(j,i) = U(i,j);
  CPPL::dgematrix U1TA = U1T * _mat;
#ifdef MC__SRED_DEBUG
  std::cout << "\nMatrix U1T * A = \n" << U1TA;
  CPPL::dgematrix Smat( nrhs, nvar );
  Smat.zero();
  for( unsigned i=0; i<rkmat; ++i )
    Smat(i,i) = S(i);
  std::cout << "\nMatrix S * VT = \n" << Smat * VT;
#endif

  CPPL::dgematrix SVTbasis( rkmat, rkmat );
  CPPL::dgbmatrix Sband;
  SVTbasis.zero();
  unsigned irk = 0;
  for( auto itwei=_mapWeightVar.crbegin(); irk<rkmat && itwei!=_mapWeightVar.crend(); ++itwei ){
    auto const& [wei,var] = *itwei;
#ifdef MC__SRED_DEBUG
    std::cout << "var = " << var << "  wei = " << wei << "  irk = " << irk << std::endl;
#endif
    for( unsigned i=0; i<rkmat; ++i )
      SVTbasis(i,irk) = U1TA(i,var);
    if( !irk ){
      _setMonElim.insert( *_mapVarMon[var] );
      irk++;
      continue;
    }
    CPPL::dgematrix cpbasis( SVTbasis );
    cpbasis.dgesvd( Sband );
    if( _rank( Sband.m, Sband.array ) == irk+1 ){
      _setMonElim.insert( *_mapVarMon[var] );
      irk++;
      continue;
    }
  }
  
  return rkmat;
}

template <typename KEY, typename COMP>
inline
unsigned
SRed<KEY,COMP>::_rank
( unsigned const dim, double const* arr )
const
{
  unsigned rank = dim;
  for( unsigned i=0; i<dim; ++i, --rank ){
#ifdef MC__SRED_DEBUG
    std::cout << "S(" << dim-i-1 << ") = " << arr[dim-i-1] << std::endl;
#endif
    if( arr[dim-i-1] > options.RKATOL
     && arr[dim-i-1] > options.RKRTOL*arr[0] ) break;
  }
  return rank;
}

template <typename KEY, typename COMP>
inline
void
SRed<KEY,COMP>::_Elim_reset
()
{
  _mapMonVar.clear();
  _mapVarMon.clear();
  _lhs.clear();
  _mapVarWeight.clear();
  _mapWeightVar.clear();
  _setMonElim.clear();
}

template <typename KEY, typename COMP>
inline
void
SRed<KEY,COMP>::_Reduc_reset
()
{
  _cleanmon();
  
  _SetMonOp.clear();
  _SetMonNew.clear();
  _SetVar.clear();
  _NewMul.clear();
  _NewDiv.clear();
  _ndxCtrUse.clear();

  _rhs.clear();
  _coef.clear();
  _new.clear();
  _RedCtr.clear();
}

#if defined(MC__USE_GUROBI)
template <typename KEY, typename COMP>
inline void
SRed<KEY,COMP>::_MIP_optimize
( std::set<unsigned> const& ndxSPol )
{
  _MIPexcpt = false;

  try{
    // Run MIP optimization for a minimal representation
    _MIP_reset();
    _MIP_encode( ndxSPol );
    _MIP_options();
    _MIP_solve();
    if( options.DISPLEVEL ) _MIP_display();
  }
  
  catch(GRBException& e){
    if( options.MIPDISPLEVEL )
      std::cout << "Error code = " << e.getErrorCode() << std::endl
                << e.getMessage() << std::endl;
    _MIPexcpt = true;
  }
}

template <typename KEY, typename COMP>
inline void
SRed<KEY,COMP>::_MIP_solve
()
{
  _GRBmodel->update();
  if( options.MIPOUTPUTFILE != "" )
    _GRBmodel->write( options.MIPOUTPUTFILE );
  fedisableexcept(FE_ALL_EXCEPT);
  _GRBmodel->set( GRB_IntAttr_ModelSense, 1 );
  _GRBmodel->optimize();
}

template <typename KEY, typename COMP>

inline void
SRed<KEY,COMP>::_MIP_encode
( std::set<unsigned> const& ndxSPol )
{
  double const coef = 1. + options.REDTOL;//1./(_SetMon.size()+1.);

  for( auto const& monop : _SetMonOp ){
    for( auto const& ctr : ndxSPol ){
      std::pair< unsigned, t_SMon const* > keyctr = std::make_pair( ctr, &monop );

      // Define Gurobi binary variables and constraints for multiplication operations
      auto [mipctrmul,insm] = _MIP_CtrMul.insert( std::make_pair( keyctr, _GRBmodel->addVar( 0., 1., -1., GRB_BINARY ) ) );
      auto itnewmul = _NewMul.find( std::make_pair( ctr, &monop ) );
      if( itnewmul !=  _NewMul.end() ){
        for( auto const& keymon : itnewmul->second ){
          auto mipnewmul = _MIP_NewMul.find( keymon );
          if( mipnewmul == _MIP_NewMul.end() )
            mipnewmul = _MIP_NewMul.insert( std::make_pair( keymon, _GRBmodel->addVar( 0., 1., coef, GRB_CONTINUOUS ) ) ).first;
          _GRBmodel->addConstr( mipnewmul->second, GRB_GREATER_EQUAL, mipctrmul->second );
        }
      }

      // Define Gurobi binary variables and constraints for division operations
      if( options.NODIV ) continue;
      auto [mipctrdiv,insd] = _MIP_CtrDiv.insert( std::make_pair( keyctr, _GRBmodel->addVar( 0., 1., -1., GRB_BINARY ) ) );
      auto itnewdiv = _NewDiv.find( std::make_pair( ctr, &monop ) );
      if( itnewdiv !=  _NewDiv.end() ){
        for( auto const& keymon : itnewdiv->second ){
          auto mipnewdiv = _MIP_NewDiv.find( keymon );
          if( mipnewdiv == _MIP_NewDiv.end() )
            mipnewdiv = _MIP_NewDiv.insert( std::make_pair( keymon, _GRBmodel->addVar( 0., 1., coef, GRB_CONTINUOUS ) ) ).first;
          _GRBmodel->addConstr( mipnewdiv->second, GRB_GREATER_EQUAL, mipctrdiv->second );
        }
      }
    }
  }
}

template <typename KEY, typename COMP>
inline void
SRed<KEY,COMP>::_MIP_display
()
{
  // Selected constraints and multiplier variables
  for( auto const& [keyctr,grbvar] : _MIP_CtrMul ){
    if( grbvar.get(GRB_DoubleAttr_X) < 0.9 ) continue;
    std::cout << "Multiply Ctr " << keyctr.first << " with Mon " << keyctr.second->display(0) << std::endl;
  }

  // New multiplier monomials
  for( auto const& [keymon,grbvar] : _MIP_NewMul ){
    if( grbvar.get(GRB_DoubleAttr_X) < 0.9 ) continue;
    std::cout << "New Mon " << keymon->display(0) << std::endl;
  }

  // Selected constraints and divider variables
  for( auto const& [keyctr,grbvar] : _MIP_CtrDiv ){
    if( grbvar.get(GRB_DoubleAttr_X) < 0.9 ) continue;
    std::cout << "Divide Ctr " << keyctr.first << " with Mon " << keyctr.second->display(0) << std::endl;
  }

  // New divider monomials
  for( auto const& [keymon,grbvar] : _MIP_NewDiv ){
    if( grbvar.get(GRB_DoubleAttr_X) < 0.9 ) continue;
    std::cout << "New Mon " << keymon->display(0) << std::endl;
  }
}

template <typename KEY, typename COMP>
inline void
SRed<KEY,COMP>::_MIP_options
()
{
  _GRBmodel->getEnv().set( GRB_IntParam_OutputFlag,        options.MIPDISPLEVEL );
  _GRBmodel->getEnv().set( GRB_IntParam_Method,            options.LPALGO );
  _GRBmodel->getEnv().set( GRB_IntParam_Presolve,          options.LPPRESOLVE );
  _GRBmodel->getEnv().set( GRB_IntParam_Threads,           options.MIPTHREADS );
  _GRBmodel->getEnv().set( GRB_IntParam_ConcurrentMIP,     options.MIPCONCURRENT );
  _GRBmodel->getEnv().set( GRB_IntParam_MIPFocus,          options.MIPFOCUS );
  _GRBmodel->getEnv().set( GRB_IntParam_NumericFocus,      options.MIPNUMFOCUS );
  _GRBmodel->getEnv().set( GRB_DoubleParam_Heuristics,     options.MIPHEURISTICS );
  _GRBmodel->getEnv().set( GRB_DoubleParam_FeasibilityTol, options.LPFEASTOL );
  _GRBmodel->getEnv().set( GRB_DoubleParam_OptimalityTol,  options.LPOPTIMTOL );
  _GRBmodel->getEnv().set( GRB_DoubleParam_MIPGap,         options.MIPRELGAP );
  _GRBmodel->getEnv().set( GRB_DoubleParam_MIPGapAbs,      options.MIPABSGAP );
  _GRBmodel->getEnv().set( GRB_DoubleParam_TimeLimit,      options.MIPTIMELIMIT );
}

template <typename KEY, typename COMP>
inline void
SRed<KEY,COMP>::_MIP_reset
()
{
  _MIP_CtrMul.clear();
  _MIP_CtrDiv.clear();
  _MIP_NewMul.clear();
  _MIP_NewDiv.clear();
  delete _GRBmodel;
  _GRBmodel = new GRBModel( *_GRBenv );
}
#endif // #if defined(MC__USE_GUROBI)

template <typename KEY, typename COMP>
inline 
void
SRed<KEY,COMP>::_sizemon
()
{
  //if( _SetVar.empty() ) throw Exceptions( Exceptions::SIZE );

  if( _nvar == _SetVar.size() && _nord == options.ORDER ) return;
  _nvar = _SetVar.size();
  _nord = options.ORDER;

  _cleanmon();
  _binom = new t_size[_nord?(_nvar+_nord-1)*(_nord+1):_nvar];
  _binom_size = std::make_pair( _nvar+_nord-1, _nord+1 );
  _set_binom();
  _posord = new t_size[_nord+2];
  _set_posord();
  _nmon = _posord[_nord+1];
  _expmon = new unsigned[_nmon*_nvar];
  _set_expmon();
}

template <typename KEY, typename COMP>
inline 
void
SRed<KEY,COMP>::_cleanmon()
{
  delete[] _expmon; _expmon = nullptr;
  delete[] _posord; _posord = nullptr;
  delete[] _binom;  _binom  = nullptr;
}

template <typename KEY, typename COMP>
inline 
void
SRed<KEY,COMP>::_set_posord
()
{
  _posord[0] = 0;
  _posord[1] = 1;
  for( unsigned i=1; i<=_nord; i++ ){
    t_size _posord_next = _posord[i] + _get_binom( _nvar+i-1, i );
    if( _posord_next > UINT_MAX )
      throw Exceptions( Exceptions::MAXSIZE );
    _posord[i+1] = _posord_next;
  }

#ifdef MC__SRED_DEBUG
  mc::display( 1, _nord+2, _posord, 1, "_posord", std::cout );
#endif
}

template <typename KEY, typename COMP>
inline 
void
SRed<KEY,COMP>::_set_expmon
()
{
  unsigned *iexp = new unsigned[_nvar] ;
  for( unsigned k=0; k<_nvar; k++ )
    _expmon[k] = iexp[k] = 0;
  for( unsigned i=1; i<=_nord; i++ ){
    for( t_size j=_posord[i]; j<_posord[i+1]; j++ ){
      _next_expmon( iexp, _nvar );
      for( unsigned k=0; k<_nvar; k++ )
        _expmon[j*_nvar+k] = iexp[k];
    }
  }
  delete[] iexp;

#ifdef MC__SRED_DEBUG
  mc::display( _nvar, _nmon, _expmon, _nvar, "_expmon", std::cout );
#endif
}

// Based on a code by John Burkardt
// http://people.sc.fsu.edu/~jburkardt/cpp_src/monomial/monomial.html
template <typename KEY, typename COMP>
inline 
void
SRed<KEY,COMP>::_next_expmon
( unsigned *x, const unsigned m )
const
{
  // Ensure that M>=1
  assert( m );

  // Find I, the index of the rightmost nonzero entry of X.
  unsigned i=0, im1, t;
  for( unsigned j=m; j>=1; j-- ){
    if ( x[j-1] ){
      i = j;
      break;
    }
  }

  // set T = X(I)
  // set X(I) to zero,
  // increase X(I-1) by 1,
  // increment X(D) by T-1.
  if( !i ){
    x[m-1] = 1;
    return;
  }
  else if( i == 1 ){
    t = x[0] + 1;
    im1 = m;
  }
  else{
    t = x[i-1];
    im1 = i - 1;
  }

  x[i-1] = 0;
  x[im1-1] = x[im1-1] + 1;
  x[m-1] = x[m-1] + t - 1;
}

template <typename KEY, typename COMP>
inline 
typename SRed<KEY,COMP>::t_size
SRed<KEY,COMP>::_loc_expmon
( const unsigned *iexp )
const
{
  unsigned ord = 0;
  for( unsigned i=0; i<_nvar; i++ ) ord += iexp[i];
#ifdef MC__SRED_CHECK
  assert( ord<_nord+2 );
#endif
  t_size pos = _posord[ord];
  
  unsigned p = _nvar ; 
  for( unsigned i=0; i<_nvar-1; i++ ){
    p--;
    for( unsigned j=0; j<iexp[i]; j++ )
      pos += _get_binom( p-1+ord-j, ord-j );
    ord -= iexp[i];
  }

  return pos;    
}

template <typename KEY, typename COMP>
inline 
void
SRed<KEY,COMP>::_set_binom
()
{
  t_size *p;
  unsigned k;
  for( unsigned i=0; i<_nvar+_nord-1; i++ ){
    p = &_binom[i*(_nord+1)];
    *p = 1;
    p++;
    *p = i+1;
    p++;
    k = ( i+1<_nord? i+1: _nord );
    for( unsigned j=2; j<=k; j++, p++ ) *p = *(p-1) * (i+2-j)/j;
    for( unsigned j=k+1; j<=_nord; j++, p++ ) *p = 0.;
  }
#ifdef MC__SRED_DEBUG
  std::cout << "nvar: " << _nvar << "  nord: " << _nord << std::endl;
  mc::display( _binom_size.second, _binom_size.first, _binom,
    _binom_size.second, "_binom", std::cout );
#endif
}

template <typename KEY, typename COMP>
inline 
typename SRed<KEY,COMP>::t_size
SRed<KEY,COMP>::_get_binom
( const unsigned n, const unsigned k )
const
{
#ifdef MC__SRED_CHECK
  assert( n<=_binom_size.first );
  assert( k<=_binom_size.second );
  assert( k<=n );
#endif
  return( n? _binom[(n-1)*_binom_size.second+k]: 1 );
}

template <typename KEY, typename COMP>
inline std::ostream&
operator<<
( std::ostream& out, SRed<KEY,COMP> const& red )
{
  const int BASIS = 0;
  const unsigned DISPLEN = red.options.DISPLEN;
  out << std::scientific << std::setprecision(DISPLEN)
      << std::right;

  // Output current monomials
  out << std::endl << "  " << red._SetMon.size() << " Monomial List: [ ";
  for( auto const& mon : red._SetMon )
    out << mon.display( BASIS ) << " ";
  out << "]" << std::endl;

  // Output multiplier/divider monomials
  out << std::endl << "  " << red._SetMonOp.size() << " Multiplier/Divider List: [ ";
  for( auto const& mon : red._SetMonOp )
    out << mon.display( BASIS ) << " ";
  out << "]" << std::endl;

  // Output reduction polynomial constraints
  if( red._RedCtr.size() ){
    out << std::endl << "  " << red._RedCtr.size() << " Reduction Constraint"
        << (red._RedCtr.size()==1? ":": "s:") << std::endl;
    for( auto const& ctr : red._RedCtr )
      out << "   0 = " << ctr.display( ctr.mapmon(), 0, red.options.DISPLEN, 1 ) << std::endl;
  }
  else
    out << std::endl << "  No Reduction Constraint" << std::endl;
  
  // Output eliminated monomials
  if( red._setMonElim.size() ){
    out << std::endl << "  " << red._setMonElim.size() << " Linearizable Monomial"
        << (red._RedCtr.size()==1? ": [ ": "s: [ ");
    for( auto const& mon : red._setMonElim )
      out << mon.display( BASIS ) << " ";
    out << "]" << std::endl;
  }
  else
    out << std::endl << "  No Linearizable Monomial" << std::endl;

  return out;
}

} // namespace mc

#endif
