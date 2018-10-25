// Copyright (C) 2015-2018 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_RLTRED Reduced RLT Generation in Factorable Expressions
\author Benoit Chachuat & OMEGA Research Group (http://www3.imperial.ac.uk/environmentenergyoptimisation)
\date 2018
\bug No known bugs.

The class mc::RLTRED in <tt>rltred.hpp</tt> identifies <b>reduction constraints</b> from a given set of factorable equality constraints; namely, constraints that are redundant because they do not cut off solutions from the original equation system, yet strengthen relaxations of the solution set. The starting point is the reduced RLT approach introduced by Liberti (2004) and Liberti & Pantelides (2006) based on the RLT approach by Sherali & Alameddine (1992), which identifies reduction constraints between existing nonlinearities, i.e. does not introduce any additional nonlinear term. This approach essentially identifies linearity embedded in a set of nonlinear constraints and may be applied as a pre-processing step in practice. 

Herein, we implement an extension of the original reduced RLT approach by: (i) enabling the creation of reduction constraint via both multiplication and division of existing linear expression by participating variables; and (ii) identifying reduction constraints via an integer-programming-based approach that supersedes the original graph-theoretical approach by Liberti & Pantelides (2006). We also discuss the case with linear inequality constraints using additional bound factor information.


\section sec_RLTRED_idea How Does the Reduced RLT Algorithm Work?

The algorithm starts by reformulating the general nonlinear constraints \f${\bf f}({\bf x}) = {\bf 0}\f$ in the following form (Smith & Pantelides, 1999):
\f{align*}
  {\bf A}\cdot{\bf z} =\ & {\bf b}\\
  z_i =\ & z_j\cdot z_k, \ \forall (i,j,k)\in\mathcal{B}\\
  z_i =\ & g_i(z_j), \ \forall (i,j)\in\mathcal{N}
\f}
where \f$z=(z_1,\ldots,z_p)\in\mathbb{R}^{n_z}\f$ is a vector of variables, comprising the original variables \f${\bf x}\f$ alongside auxiliary variables; \f${\bf A}\in\mathbb{R}^{n_\ell\times n_z}\f$ is a full row-rank matrix, and \f${\bf b}\in\mathbb{R}^{n_\ell}\f$; \f$\mathcal{N}\f$ is a set of index pairs \f$1\leq i,j\leq n_z\f$, and \f$g_i\f$ are nonlinear continuous univariate functions; and \f$\mathcal{B}\f$ is a set of index triplets \f$1\leq i,j,k\leq n_z\f$ with \f$j\leq k\f$. Notice that, in addition to bilinear terms, the following nonlinearities may be accounted for as part of the set \f$\mathcal{B}\f$:
- fractional terms, \f$z_i = \frac{z_j}{z_k}\f$ so that \f$(j,i,k)\in\mathcal{B}\f$;
- square terms, \f$z_i = z_j^2\f$ so that \f$(i,j,j)\in\mathcal{B}\f$;
- square-root terms, \f$z_i = \sqrt{z_j}\f$ so that \f$(j,i,i)\in\mathcal{B}\f$.
.
The next step entails the identification of those subsets of linear constraints which, when multiplied or divided by a given variable subset, lead to the creation of fewer new nonlinear terms than extra constraints. Liberti & Pantelides (2006) proposed a graph-theoretical approach to making this assignment. An equivalent approach consists in constructing and solving the following integer linear program (ILP):
\f{align*}  \min_{{\bf b},{\bf f},{\bf d},{\bf m}} &\ \sum_{j=1}^{n_z} \left[ (1+\epsilon)\sum_{i=1}^{n_z} \left(b_{i,j}+f_{i,j}\right) - \sum_{\lambda=1}^{n_\ell} \left(d_{\lambda,j}+m_{\lambda,j}\right) \right]\\
\text{s.t.}\ &\ b_{i,j} \geq m_{\lambda,j}, \ \ \forall i\in\mathcal{M}_{\lambda,j}, \forall j,\lambda\\
             &\ f_{i,j} \geq d_{\lambda,j}, \ \ \forall i\in\mathcal{D}_{\lambda,j}, \forall j,\lambda\\
             &\ b_{i,j},f_{i,j},m_{\lambda,j},d_{\lambda,j}\in\{0,1\},\ \ \forall i,j,\lambda\,.
\f}

- The binary variables have the following interpretation:
\f{align*}
m_{\lambda,j} :=\ & \left\{\begin{array}{ll} 1 & \text{if linear constraint $\lambda\in 1\ldots n_\ell$ is multiplied by variable $z_j, j\in 1\ldots n_z$}\\ 0 & \text{otherwise} \end{array}\right.\\
d_{\lambda,j} :=\ & \left\{\begin{array}{ll} 1 & \text{if linear constraint $\lambda\in 1\ldots n_\ell$ is divided by variable $z_j, j\in 1\ldots n_z$}\\ 0 & \text{otherwise} \end{array}\right.\\
b_{i,j} :=\ & \left\{\begin{array}{ll} 1 & \text{if the bilinear term $z_i\cdot z_j$ is created, i.e. $\not\exists k: (k,i,j)\in\mathcal{B}\vee (k,j,i)\in\mathcal{B}$}\\ 0 & \text{otherwise} \end{array}\right.\\
f_{i,j} :=\ & \left\{\begin{array}{ll} 1 & \text{if the fractional term $\frac{z_i}{z_j}$ is created, i.e. $\not\exists k: (j,i,k)\in\mathcal{B}\vee (j,k,i)\in\mathcal{B}$}\\ 0 & \text{otherwise} \end{array}\right.
\f}

- The index sets \f$\mathcal{M}_{\lambda,j}\f$ and \f$\mathcal{D}_{\lambda,j}\f$ associated with each linear constraint \f$\lambda\in 1\ldots n_\ell\f$ and each multiplier or divider variable \f$z_j, j\in 1\ldots n_z\f$, respectively, are given by:
\f{align*}
\mathcal{M}_{\lambda,j} :=\ & \left\{ i \middle| a_{\lambda, i} \neq 0, \not\exists k: (k,i,j)\in\mathcal{B}\vee (k,j,i)\in\mathcal{B} \right\}\\
\mathcal{D}_{\lambda,j} :=\ & \left\{ i \middle| a_{\lambda, i} \neq 0, \not\exists k: (i,j,k)\in\mathcal{B}\vee (i,k,j)\in\mathcal{B} \right\}
\f}
This way, given \f$i\in\mathcal{M}_{\lambda,j}\f$ for some \f$j\in 1\ldots n_z\f$ and \f$\lambda\in 1\ldots n_\ell\f$, the cut \f$b_{i,j} \geq m_{\lambda,j}\f$ accounts for the creation of the new bilinear term \f$z_i\cdot z_j\f$ on multiplying constraint \f$\lambda\f$ with variable \f$z_j\f$; likewise, given \f$i\in\mathcal{D}_{\lambda,j}\f$ for some \f$j\in 1\ldots n_z\f$ and \f$\lambda\in 1\ldots n_\ell\f$, the cut \f$f_{i,j} \geq d_{\lambda,j}\f$ accounts for the creation of the new fractional term \f$\frac{z_i}{z_j}\f$ on dividing constraint \f$\lambda\f$ with variable \f$z_j\f$.

- Finally, the role of the factor \f$(1-\epsilon)\f$ multiplying the binary variables \f$b_{i,j}\f$ and \f$f_{i,j}\f$ in the objective function is to penalize/exclude assignments introducing as many new constraints as bilinear terms -- which would otherwise add up to zero in the objective function and could give rise to a large number of (globally optimal) solutions for the ILP. The objective in doing so is therefore to determine those assignments yielding a valid reduced RLT constraint set (Liberti & Pantelides, 2006). 
.

Notice that taking all of the binary variables \f$b_{i,j},f_{i,j},m_{\lambda,j},d_{\lambda,j}\f$ equal to 0 gives a feasible solution to the ILP with an objective of 0. Of interest here are of course optimal solution with a strictly negative objective, so that exactly \f$M-N>0\f$ reduction constraints could be identified, with \f$M:=\sum_{j=1}^{n_z} \sum_{\lambda=1}^{n_\ell} (d_{\lambda,j}+m_{\lambda,j})\f$ and \f$N:=\sum_{j=1}^{n_z} \sum_{i=1}^{n_z} (b_{i,j}+f_{i,j})\f$. In order to recover their symbolic expressions, we proceed by first writing the complete set of new constraints obtained from the multiplication or division of existing linear constraints with variables,
\f{align*}
\forall (\lambda,j): m_{\lambda,j}=1, \quad & \sum_{i=1}^{n_z} a_{\lambda,i}\cdot z_i\cdot z_j = b_{\lambda}\cdot z_j\\
\forall (\lambda,j): d_{\lambda,j}=1, \quad & \sum_{i=1}^{n_z} a_{\lambda,i}\cdot \frac{z_i}{z_j} = b_{\lambda}\cdot \frac{1}{z_j}
\f}
These constraints can be rewritten in matrix form as
\f{align*}
  {\bf R}\cdot{\bf n}({\bf z}) =\ & {\bf m}({\bf z})
\f}
where:
- \f${\bf m}({\bf z})\f$ is an \f$M\f$-dimensional vector, with entries
\f{align*}
  {\bf m}({\bf z}) &\ := \left( \begin{array}{cl} \vdots & \\ b_{\lambda}\cdot z_j - \sum_{\{i|\exists k:(k,i,j)\in\mathcal{B}\}} a_{\lambda,i}\cdot z_i\cdot z_j & \forall (\lambda,j): m_{\lambda,j}=1 \\ \vdots & \\ b_{\lambda}\cdot \frac{1}{z_j} - \sum_{\{i|\exists k:(i,k,j)\in\mathcal{B}\}} a_{\lambda,i}\cdot \frac{z_i}{z_j} & \forall (\lambda,j): d_{\lambda,j}=1 \\ \vdots & \end{array} \right)
\f}
- \f${\bf n}({\bf z})\f$ is an \f$N\f$-dimensional vector, comprised of the newly created bilinear and fractional terms
\f{align*}
  {\bf n}({\bf z}) &\ := \left( \begin{array}{cl} \vdots & \\ z_i\cdot z_j & \forall(i,j):b_{i,j}=1 \\ \vdots & \\ \frac{z_i}{z_j} & \forall(i,j):f_{i,j}=1 \\ \vdots & \end{array} \right)
\f}
- \f${\bf R}\f$ is an \f$M\f$-by-\f$N\f$ real matrix
.
Finally, we perform an LU factorization of the coefficient matrix \f${\bf R}\f$, such that \f${\bf P}\cdot{\bf R}={\bf L}\cdot{\bf U}\f$. The bottom \f$M-N>0\f$ rows of the vector \f${\bf s}({\bf z})\f$ such that \f${\bf L}\cdot{\bf s}({\bf z})={\bf P}\cdot{\bf m}({\bf z})\f$ -- as obtained for instance from a simple forward elimination -- are equal to zero by construction, and yield valid reduced RLT constraints,
\f{align*}
  s_{N+1}({\bf z}) &\ = 0\\
  &\ \vdots\\
  s_{M}({\bf z}) &\ = 0
\f}
For large-scale problem, one may use efficient sparse implementations of LU factorization, based for instance on a bordered block decomposition of the original coefficient matrix -- see, e.g., available routines within the <A href="http://www.hsl.rl.ac.uk/catalogue/">Harwell Scientific Library (HSL)</A>.


\section sec_RLTRED_ineq What If My Problem Has Linear Inequality Constraints?

Consider the case that linear inequality constraints are present in the form of
\f{align*}
  {\bf A}\cdot{\bf z} \leq\ & {\bf b}
\f}
and the variables are bounded, \f$z^{\rm L}\leq z\leq z^{\rm U}\f$. A very similar procedure could be used to identify pairs of inequality reduction constraints upon multiplying the linear inequality constraints with the non-negative factors \f$z^{\rm U}-z\f$ and \f$z-z^{\rm L}\f$. Division with such factors to generate extra constraints no longer seems possible in this case, however.

Having determined optimal values for the binary variables \f$b_{i,j},m_{\lambda,j}\f$ like previously, new inequality constraints obtained from the multiplication of existing linear inequality constraints with the non-negative factors \f$z^{\rm U}-z\f$ and \f$z-z^{\rm L}\f$ read
\f{align*}
& \forall (\lambda,j): m_{\lambda,j}=1,\\ & b_{\lambda}\cdot (z_j-z_j^{\rm U}) + \sum_{i=1}^{n_z} a_{\lambda,i}\cdot z_i\cdot z_j^{\rm U}~\leq~\sum_{i=1}^{n_z} a_{\lambda,i}\cdot z_i\cdot z_j~\leq~b_{\lambda}\cdot (z_j-z_j^{\rm L}) + \sum_{i=1}^{n_z} a_{\lambda,i}\cdot z_i\cdot z_j^{\rm L}
\f}
Finally, one can apply a similar LU factorization strategy as with linear equality constraints in order to extract \f$M-N\f$ valid reduction inequality constraint pairs in terms of the existing linear and nonlinear terms only.


\section sec_RLTRED_use How Do I Search for Redundant Constraints in a Factorable Expression

For illustration, suppose we want to determine reduction constraints for the following set of equality constraints:
\f{align*}
  {\bf 0} = {\bf f}({\bf x}) = \left(\begin{array}{c} x_0 + x_1 - 1\\ x_0^2 - x_1^2 - x_2\end{array}\right)
\f}

The constructions require the header file <tt>rltred.hpp</tt> to be included:

\code
      #include "rltred.hpp"
\endcode

A DAG of the factorable function is first created:

\code
      mc::FFGraph DAG;
      const unsigned NX = 3, NF = 2;
      mc::FFVar X[NX];
      for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
      mc::FFVar F[NF];
      F[0] = X[0] + X[1] - 1.;
      F[1] = sqr(X[0]) - sqr(X[1]) - X[2];
      std::cout << DAG;
\endcode

The last line displays the following information about the factorable constraint DAG:

\verbatim
    DAG VARIABLES:
      X0	 => { Z0 Z4 }
      X1	 => { Z0 Z3 }
      X2	 => { Z6 }

    DAG INTERMEDIATES:
      Z0	<=  X0 + X1		 => { Z2 }
      Z2	<=  Z0 + Z1		 => { }
      Z3	<=  SQR( X1 )	 => { Z5 }
      Z4	<=  SQR( X0 )	 => { Z5 }
      Z5	<=  Z4 - Z3		 => { Z6 }
      Z6	<=  Z5 - X2		 => { }
      Z1	<=  -1(D)		 => { Z2 }
\endverbatim

Reduction constraint identification starts by creating an instance of the class mc::RLTRED for the DAG at hand. Various options may be specified at this point using the public member mc::RLTRED::options; e.g.:
- mc::RLTRED::Options::METHOD enables the selection of ILP-based search, instead of the default graph-based approach
- mc::RLTRED::Options::LEVEL enables simultaneous search in all the participating variables, instead of the default sequential, variable-by-variable strategy
- mc::RLTRED::Options::NODIV enables multiplying constraints with variables only, as opposed to both multiplication and division
- mc::RLTRED::Options::MIPDISPLAY enables display of the MIP solver
.

\code
      mc::RLTRED RRLT( &DAG );
      RRLT.options.METHOD     = mc::RLTRED::Options::ILP;
      RRLT.options.LEVEL      = mc::RLTRED::Options::FULLSIM;
      RRLT.options.NODIV      = true;
      RRLT.options.MIPDISPLAY = 1;
\endcode

The reduction constraints algorithm is called via the method mc::RLTRED::search. A vector of valid reduction constraint found can be retreived using the method mc::RLTRED::constraints, and then displayed as shown below.

\code
      RRLT.search( NF, F );
      auto FRED = RRLT.constraints();
      for( auto it=FRED.begin(); it!=FRED.end(); ++it ){
        std::ostringstream ostr; ostr << " OF REDUCTION CONSTRAINT " << **it;
        DAG.output( DAG.subgraph( 1, *it ), ostr.str() );
      }
\endcode

The following output is generated in this case:

\verbatim
      Found incumbent of value 0.000000 after 0.00 sec. (0.01 ticks)
      Tried aggregator 4 times.
      MIP Presolve eliminated 28 rows and 39 columns.
      Aggregator did 42 substitutions.
      Reduced MIP has 16 rows, 15 columns, and 32 nonzeros.
      Reduced MIP has 15 binaries, 0 generals, 0 SOSs, and 0 indicators.
      Presolve time = 0.00 sec. (0.16 ticks)
      Probing time = 0.00 sec. (0.01 ticks)
      Tried aggregator 1 time.
      Reduced MIP has 16 rows, 15 columns, and 32 nonzeros.
      Reduced MIP has 15 binaries, 0 generals, 0 SOSs, and 0 indicators.
      Presolve time = 0.00 sec. (0.02 ticks)
      Probing time = 0.00 sec. (0.01 ticks)
      Clique table members: 16.
      MIP emphasis: balance optimality and feasibility.
      MIP search method: dynamic search.
      Parallel mode: deterministic, using up to 4 threads.
      Root relaxation solution time = 0.00 sec. (0.02 ticks)

              Nodes                                         Cuts/
         Node  Left     Objective  IInf  Best Integer    Best Bound    ItCnt     Gap

      *     0+    0                           -0.9997       -7.9994           700.18%
            0     0        cutoff             -0.9997                      9    0.00%

      Root node processing (before b&c):
        Real time             =    0.01 sec. (0.27 ticks)
      Parallel b&c, 4 threads:
        Real time             =    0.00 sec. (0.00 ticks)
        Sync time (average)   =    0.00 sec.
        Wait time (average)   =    0.00 sec.
                          ------------
      Total (root+branch&cut) =    0.01 sec. (0.27 ticks)
        #reduction constraints: 1

      FACTORS IN SUBGRAPH OF REDUCTION CONSTRAINT Z12:
        X1	<=  VARIABLE
        Z3	<=  SQR( X1 )	
        X0	<=  VARIABLE
        Z8	<=  - X1	
        Z9	<=  X0 + Z8	
        Z10	<=  Z3 + Z9	
        Z4	<=  SQR( X0 )	
        Z11	<=  - Z4	
        Z12	<=  Z10 + Z11	
\endverbatim

The algorithm has therefore identified the following reduction constraint:
\f{align*}
  x_0 - x_1 + x_0^2 - x_1^2 = 0
\f}
It is not hard to check that the solution set of the former constraint is \f$\left\{x_0,x_1\middle| x_0+x_1=1 \vee x_0-x_1=0 \right\}\f$, which is indeed redundant with the original constraint \f$ x_0+x_1=1\f$.

\section sec_FRLTRED_refs References

- Liberti, L. <A href="https://doi.org/10.1111/j.1475-3995.2004.00438.x">Reduction constraints for the global optimization of NLPs</A>, <i>International Transactions in Operation Research</i>, <b>11</b>(1):33-41, 2004.
- Liberti, L., Pantelides, C. C. <A href="https://doi.org/10.1007/s10898-006-9005-4">An exact reformulation algorithm for large
nonconvex NLPs involving bilinear terms</A>, <i>Journal of Global Optimization</i>, <b>36</b>(2):161-189, 2006.
- Sherali, H., Alameddine, A. <A href="https://doi.org/10.1007/BF00122429">A new reformulation-linearization technique for
bilinear programming problems</A>, <i>Journal of Global Optimization</i>, <b>2</b>:379-410.
- Smith, E. M. B., Pantelides, C. C., <A href="http://dx.doi.org/10.1016/S0098-1354(98)00286-5">A symbolic reformulation/spatial branch-and-bound algorithm for the global optimisation of nonconvex MINLPs</A>, <i>Computers & Chemical Engineering</i>, <b>23</b>(4/5):457-478, 1999.
.
*/

//TODO: 
//- [TO DO] Documentation
//- [TO DO] Handle product (multilinear) terms
//- [TO DO] Detect that products of the same variable is a power term
//- [TO DO] Enable sparse LU factorization, e.g. using MA33 from HSL, instead of dgetrf from LAPACK
//- [TO DO] Implement elimination of nonlinear terms using the reduction constraints
//BUGS: 
//- [DONE] For division reduction, constants need to be taken into account
//- [DONE] Limit FULLSEQ / FULLSIM to nonlinear terms to avoid trivial / duplicate constraints

#ifndef MC__RLTRED_HPP
#define MC__RLTRED_HPP

#include "ffunc.hpp"
#include "spolyexpr.hpp"
#include "mclapack.hpp"
#if defined(MC__USE_CPLEX)
  #include "ilcplex/ilocplex.h"
#elif defined(MC__USE_GUROBI)
  #include "gurobi_c++.h"
#endif

#undef MC__RLTRED_DEBUG

#if defined(MC__USE_GUROBI)
extern "C"{
  #include <fenv.h>
  int fedisableexcept( int );
}
#endif

namespace mc
{

//! @brief C++ base class for ordering new potential bilinear/ fractional terms in RLT procedure
////////////////////////////////////////////////////////////////////////
//! mc::lt_RLTNew is a C++ structure for ordering new potential bilinear
//! / fractional terms in an RLT procedure.
////////////////////////////////////////////////////////////////////////
struct lt_RLTNew
////////////////////////////////////////////////////////////////////////
{
  typedef std::pair< const FFVar*, typename FFOp::TYPE > t_RLTVar;
  typedef std::pair< const t_RLTVar, const FFVar* > t_New;

  bool operator()
    ( const t_New& Term1, const t_New& Term2 ) const
    {
      // Order terms w.r.t. their RLT operation types first
      if( Term1.first.second < Term2.first.second ) return true;
      if( Term1.first.second > Term2.first.second ) return false;

      switch( Term1.first.second ){
       // Case of bilinear terms
       case FFOp::TIMES:
        // Order terms w.r.t. their lower variable next
        if( lt_FFVar()( Term1.first.first, Term1.second ) ){
          if( lt_FFVar()( Term2.first.first, Term2.second ) ){
            if( lt_FFVar()( Term1.first.first, Term2.first.first ) ) return true;
            if( lt_FFVar()( Term2.first.first, Term1.first.first ) ) return false;
            return lt_FFVar()( Term1.second, Term2.second );
          }
          if( lt_FFVar()( Term1.first.first, Term2.second ) ) return true;
          if( lt_FFVar()( Term2.second, Term1.first.first ) ) return false;
          return lt_FFVar()( Term1.second, Term2.first.first );
        }
        if( lt_FFVar()( Term2.first.first, Term2.second ) ){
          if( lt_FFVar()( Term1.second, Term2.first.first ) ) return true;
          if( lt_FFVar()( Term2.first.first, Term1.second ) ) return false;
          return lt_FFVar()( Term1.first.first, Term2.second );
        }
        if( lt_FFVar()( Term1.second, Term2.second ) ) return true;
        if( lt_FFVar()( Term2.second, Term1.second ) ) return false;
        return lt_FFVar()( Term1.first.first, Term2.first.first );

       // Case of fractional terms
       case FFOp::DIV:
        // Order terms w.r.t. their RLT variable next
        if( lt_FFVar()( Term1.first.first, Term2.first.first ) ) return true;
        if( lt_FFVar()( Term2.first.first, Term1.first.first ) ) return false;
        // Order terms w.r.t. their linear variable at last
        return lt_FFVar()( Term1.second, Term2.second );

       default:
        throw std::runtime_error( "Internal error" );
      }
    }
};

//! @brief C++ base class for ordering edges in RLT bi-graph
////////////////////////////////////////////////////////////////////////
//! mc::lt_RLTEdge is a C++ structure for ordering reduction constraints
//! in an RLT set.
////////////////////////////////////////////////////////////////////////
struct lt_RLTEdge
////////////////////////////////////////////////////////////////////////
{
  typedef std::pair< const FFVar*, typename FFOp::TYPE > t_RLTVar;
  typedef std::pair< const t_RLTVar, const FFOp* > t_Edge;

  bool operator()
    ( const t_Edge& Edge1, const t_Edge& Edge2 ) const
    {
      // Order RLT edges w.r.t. their RLT operation types first
      if( Edge1.second < Edge2.second ) return true;
      if( Edge1.second > Edge2.second ) return false;
      // Order RLT edges w.r.t. their RLT variable next
      if( lt_FFVar()( Edge1.first.first, Edge2.first.first ) ) return true;
      if( lt_FFVar()( Edge2.first.first, Edge1.first.first ) ) return false;
      // Order RLT edges w.r.t. their RLT terms at last
      return( Edge1.first.second < Edge2.first.second );
    }
};

//! @brief C++ base class for reduced RLT generation in a DAG
////////////////////////////////////////////////////////////////////////
//! mc::RLTRED is a C++ base class implementing the reduced RLT
//! approach by Liberti et al. to identifying redundant constraints
//! in a DAG based on linear equality constraints.
////////////////////////////////////////////////////////////////////////
class RLTRED
{
public:
  typedef typename FFVar::pt_idVar pt_idVar;
  typedef std::set< const FFVar*, lt_FFVar > t_Vars;
  typedef std::pair< const FFOp*, unsigned int > pt_Op;
  typedef std::list< const FFOp* > l_Ops;
  typedef std::set< const FFOp*,  lt_FFOp > t_Ops;
  typedef std::multimap< const FFVar*, const FFOp*, lt_FFVar > t_Edges;
  typedef std::pair< const FFVar*, typename FFOp::TYPE > t_RLTVar;
  typedef std::set< std::pair< const t_RLTVar, const FFOp* >, lt_RLTEdge > t_RLTRed;
#if defined(MC__USE_CPLEX)
  typedef std::map< std::pair< const t_RLTVar, const FFOp* >, IloNumVar, lt_RLTEdge > t_ILPCtr;
  typedef std::map< std::pair< const t_RLTVar, const FFVar* >, IloNumVar, lt_RLTNew > t_ILPVar;
#elif defined(MC__USE_GUROBI)
  typedef std::map< std::pair< const t_RLTVar, const FFOp* >, GRBVar, lt_RLTEdge > t_ILPCtr;
  typedef std::map< std::pair< const t_RLTVar, const FFVar* >, GRBVar, lt_RLTNew > t_ILPVar;
#endif

private:
  //! @brief Pointer to DAG
  FFGraph* _dag;

  //! @brief Subset of linear operations
  t_Ops _linTerms;
  //! @brief Subset of variables participating in linear operations
  t_Vars _linVars;
  //! @brief Map of linear operations and participating variables
  t_Edges _linEdges;
  //! @brief Subset of variables for RLT constraint search
  t_Vars _RLTVars;
  //! @brief Subset of edges for RLT constraint search
  t_Edges _RLTEdges;

  // subset of bilinear terms
  t_Ops _bilTerms;
  // subset of fractional terms
  t_Ops _divTerms;
  // subset of inverse terms
  t_Ops _invTerms;
  // subset of square terms
  t_Ops _sqrTerms;
  // subset of square root terms
  t_Ops _sqrtTerms;
  // subset of integer pow terms
  t_Ops _ipowTerms;

  //! @brief Map of RLT constraints and corresponding multiplier/divider variable
  t_RLTRed _RLTRed;
  //! @brief Coefficient matrix of "hidden" linear (reduction) expressions
  CPPL::dgematrix _coefRed;
  //! @brief New variables/operations participating in "hidden" linear (reduction) expressions
  std::map< const FFOp*, unsigned, lt_FFOp > _varRed;
  //! @brief Existing (linear combinations of) variables participating in "hidden" linear (reduction) expressions
  std::vector<FFVar*> _rhsRed;
  //! @brief DAG variables representing "hidden" linear (reduction) expressions in DAG
  std::vector<const FFVar*> _ctrRed;

  //! @brief Whether a given RLT variable has been assigned to a linear term and which term
  std::vector<pt_Op> _VarAssigned;
  //! @brief Whether a given RLT variable has been visited
  std::vector<unsigned> _VarVisited;
  //! @brief Whether a given linear term has been visited
  std::vector<unsigned> _TermVisited;
  
#if defined(MC__USE_CPLEX) || defined(MC__USE_GUROBI)
  //! @brief map of binary variables in variable operation in ILP
  t_ILPVar _ILPvar;
  //! @brief map of binary variables for constraint operation in ILP
  t_ILPCtr _ILPctr;
  //! @brief whether the ILP solver has sent an exception
  bool _ILPexcpt;
#endif
#if defined(MC__USE_CPLEX)
  //! @brief CPLEX environment for ILP
  IloEnv* _ILOenv;
  //! @brief CPLEX model for ILP
  IloModel* _ILOmodel;
  //! @brief CPLEX object for ILP
  IloCplex* _ILOcplex;
  //! @brief CPLEX objective for ILP
  IloExpr* _ILPobj;
#elif defined(MC__USE_GUROBI)
  //! @brief GUROBI environment for ILP
  GRBEnv* _GRBenv;
  //! @brief GUROBI model for ILP
  GRBModel* _GRBmodel;
#endif

public:
  
  //! @brief Default Constructor
  RLTRED
    ( FFGraph* dag )
    : _dag( dag )
    {
#if defined(MC__USE_CPLEX)
      _ILOenv   = new IloEnv;
      _ILOmodel = 0;
      _ILOcplex = 0;
      _ILPobj = 0;
#elif defined(MC__USE_GUROBI)
      _GRBenv   = new GRBEnv();
      _GRBmodel = 0;
#endif
    }
  //! @brief Destructor
  ~RLTRED()
    {
      for( auto itv=_varRed.begin(); itv!=_varRed.end(); ++itv ) delete itv->first;
      for( auto itv=_rhsRed.begin(); itv!=_rhsRed.end(); ++itv ) delete *itv;
#if defined(MC__USE_CPLEX)
      delete _ILPobj;
      delete _ILOmodel;
      delete _ILOcplex;
      _ILOenv->end();
      delete _ILOenv;
#elif defined(MC__USE_GUROBI)
      delete _GRBmodel;
      delete _GRBenv;
#endif
    }

  //! @brief Return a map of reduction constraints for the expressions in lhs
  t_RLTRed& search
    ( const std::vector<const FFVar*>& vlhs, const bool add2dag=true );

  //! @brief Return a map of reduction constraints for the expressions in lhs
  t_RLTRed& search
    ( const unsigned nlhs, const FFVar* lhs, const bool add2dag=true );

  //! @brief Return a vector of pointers to the reduction "hidden" constraints found
  const std::vector<const FFVar*>& constraints
    ()
    const
    { return _ctrRed; }

  //! @brief Options of mc::RLTRED
  struct Options
  {
    //! @brief Constructor
    Options():
      LEVEL(PRIMSEQ), METHOD(GRAPH), NODIV(false), VARTOL(1e-4), 
      LPALGO( RLTRED::LPALGO_DEFAULT ), LPPRESOLVE(-1),
      LPFEASTOL(1e-9), LPOPTIMTOL(1e-9), MIPRELGAP(1e-7), MIPABSGAP(1e-7),
      DISPLAY(0), MIPDISPLAY(0), MIPFILE(""), MAXCPU(7200)
      {}
    //! @brief Enumeration type for reduced RLT strategy
    enum RLTVARS{
      PRIMSEQ=0,//!< Reduced RLT with primary variables as operands only treated sequentially
      FULLSEQ,	//!< Reduced RLT with primary and auxiliary variables as operands treated sequentially
      PRIMSIM,	//!< Reduced RLT with primary variables as operands only treated simultaneously
      FULLSIM	//!< Reduced RLT with primary and auxiliary variables as operands treated simultaneously
    };
    //! @brief Enumeration type for reduced RLT strategy
    enum RLTMETH{
      GRAPH=0,  //!< Reduced RLT using graph theory
#if defined(MC__USE_CPLEX) || defined(MC__USE_GUROBI)
      ILP 	    //!< Reduced RLT using integer programming (ILP)
#endif
    };
    //! @brief Reduced RLT variable selection
    RLTVARS LEVEL;
    //! @brief Reduced RLT method selection
    RLTMETH METHOD;
    //! @brief Whether or not to consider divider variables in addition to mutliplier variables
    bool NODIV;
    //! @brief Tolerance on variable to penalize operations that do not create redundancy
    double VARTOL;
    //! @brief LP algorithm
    int LPALGO;
    //! @brief LP presolve
    int LPPRESOLVE;
    //! @brief Tolerance on LP feasibility
    double LPFEASTOL;
     //! @brief Tolerance on LP optimality
    double LPOPTIMTOL;
    //! @brief Tolerance on relative MIP gap
    double MIPRELGAP;
    //! @brief Tolerance on absolute MIP gap
    double MIPABSGAP;
    //! @brief Display level
    unsigned DISPLAY;
    //! @brief Display level for MIP
    int MIPDISPLAY;
    //! @brief Name of output file for optimization model
    std::string MIPFILE;
    //! @brief Maximum run time (seconds)
    double MAXCPU;
  } options;

  //! @brief Structure holding current statistics
  struct Stats{
    void reset()
      { tGRAPH = tILPSET = tCTRGEN = tCTRDAG = tTOTAL = 0.; } // = tILPAUX
    void display
      ( std::ostream&os=std::cout )
      { os << std::fixed << std::setprecision(2)
           << "#  TOTAL:  " << tTOTAL  << " CPU SEC" << std::endl
           << "#  GRAPH:  " << tGRAPH  << " CPU SEC" << std::endl
           //<< "#  ILPAUX:  " << tILPAUX  << " CPU SEC" << std::endl
           << "#  ILPSET: " << tILPSET << " CPU SEC" << std::endl
           << "#  CTRGEN: " << tCTRGEN << " CPU SEC" << std::endl
           << "#  CTRDAG: " << tCTRDAG << " CPU SEC" << std::endl;
        os << std::endl; }
    double tTOTAL;
    double tGRAPH;
    //double tILPAUX;
    double tILPSET;
    double tCTRGEN;
    double tCTRDAG;
  } stats;

  //! @brief Default option for LP solver
#if defined(MC__USE_CPLEX)
  static const int LPALGO_DEFAULT = 0;
#elif defined(MC__USE_GUROBI)
  static const int LPALGO_DEFAULT = -1;
#endif

private:
  //! @brief Define subsets and map of linear operations and participating variables
  void _extract_linearity
    ( const l_Ops& Ops, const std::vector<const FFVar*>& vlhs );
  //! @brief Define subsets of bilinear, fractional, square and square root terms
  void _target_nonlinearity
    ( const l_Ops& Ops );
  //! @brief Search for redution constraints using a bi-partite graph for each candidate variable sequentially
  void _search_sequential
    ( const l_Ops& Ops );
  //! @brief Search for redution constraints using a bi-partite graph for all candidate variable simultaneously
  void _search_simultaneous
    ( const l_Ops& Ops );
  //! @brief Identify variables not yet participating in any nonlinear terms with variable <a>idMult</a>, and create corresponding bigraph edges
  void _bigraph_RLT
    ( const t_RLTVar& VarRed, t_Vars& RLTVars, t_Edges& RLTEdges )
    const;
  //! @brief Identify variables not yet participating in any nonlinear terms multiplied by variable <a>idMult</a>, and create corresponding bigraph edges
  void _bigraph_RLT_mult
    ( const FFVar* VarMult, t_Vars &linVars, t_Edges &linEdges )
    const;
  //! @brief Identify variables not yet participating in any nonlinear terms divided by variable <a>idDiv</a>, and create corresponding bigraph edges
  void _bigraph_RLT_div
    ( const FFVar* VarDiv, t_Vars &linVars, t_Edges &linEdges )
    const;
  //! @brief Identify valid reduction constraints with the candidate reduction variable <a>VarRed</a>
  void _reduction_RLT
    ( const t_RLTVar& VarRed );
  //! @brief Determine whether an augmented path can be found that emanates from the linear term <a>pOp<\a>
  bool _augpath_RLT
    ( const pt_Op& pOp );

  //! @brief Add reduced constraint expressions to DAG and display
  void _add_to_dag
    ( const std::vector<const FFVar*>& vlhs );
  //! @brief Process all terms in reduced constraint
  void _process_constraint
    ( const unsigned iRow, const std::vector<const FFVar*>& vlhs, const t_RLTVar& RLTVar,
      const FFOp*RLTOp );
  //! @brief Process one term in reduced constraint
  void _process_term
    ( const unsigned iRow, const FFOp::TYPE& RLTOp, const FFVar& RLTVar,
      const FFVar& linVar, const double coef );

  //! @brief Create subset of operations of given type
  t_Ops _subset_op
    ( const l_Ops& Ops, const unsigned int nOp, const typename FFOp::TYPE*typeOp )
    const;
  //! @brief Create subset of operations of given type
  t_Ops _subset_op
    ( const l_Ops& Ops, const typename FFOp::TYPE&typeOp )
    const;
  //! @brief Create subset of variables participating in operations Ops and sunmap with their defining operations 
  void _subset_var_edge
    ( const std::vector<const FFVar*>& vlhs, const t_Ops& Ops,
      t_Vars& Vars, t_Edges& Edges )
    const;
  //! @brief Test if a variable is corresponds to a dependent
  bool _is_rhs
    ( const std::vector<const FFVar*>& vlhs, const FFVar& varRes )
    const;

#if defined(MC__USE_CPLEX)
  //! @brief Add new term binary variable to ILP model
  t_ILPVar::iterator _set_ILPvar
    ( const t_RLTVar& VarRed, const FFVar* pVarLin );
  //! @brief Add new constraint binary variable to ILP model
  t_ILPCtr::iterator _set_ILPctr
    ( const t_RLTVar& VarRed, const FFOp* pOpLin );
  //! @brief Add new edge constraint to ILP model
  void _add_ILPedge
    ( const IloNumVar& varbin, const IloNumVar& ctrbin );
#elif defined(MC__USE_GUROBI)
  //! @brief Add new term binary variable to ILP model
  t_ILPVar::iterator _set_ILPvar
    ( const t_RLTVar& VarRed, const FFVar* pVarLin );
  //! @brief Add new constraint binary variable to ILP model
  t_ILPCtr::iterator _set_ILPctr
    ( const t_RLTVar& VarRed, const FFOp* pOpLin );
  //! @brief Add new edge constraint to ILP model
  void _add_ILPedge
    ( const GRBVar& varbin, const GRBVar& ctrbin );
#endif
#if defined(MC__USE_CPLEX) || defined(MC__USE_GUROBI)
  //! @brief Append binary variables, cost contributions and constraints to ILP model
  void _append_ILPmodel
    ( const t_RLTVar& VarRed, t_Vars& RLTVars, t_Edges& RLTEdges );
  //! @brief Reset ILP model
  void _reset_ILPmodel
    ();
  //! @brief Set options of ILP model
  void _set_ILPoptions
    ();
  //! @brief Solve ILP model
  void _solve_ILPmodel
    ();
  //! @brief Append valid reduction constraints from the ILP model solution
  void _extract_RLT
    ();
#endif
};

inline typename RLTRED::t_RLTRed&
RLTRED::search
( const unsigned nlhs, const FFVar* lhs, const bool add2dag )
{
  std::vector<const FFVar*> vlhs;
  for( unsigned i=0; i<nlhs; i++ ) vlhs.push_back( lhs+i );
  return search( vlhs, add2dag );
}

inline typename RLTRED::t_RLTRed&
RLTRED::search
( const std::vector<const FFVar*>& vlhs, const bool add2dag )
{
  stats.reset();
  stats.tTOTAL -= cpuclock();

  // Reset reduction constraints
  _RLTRed.clear();  _ctrRed.clear();

  // Generate list of operation for constraint LHS and particpating variables
  stats.tGRAPH -= cpuclock();
  l_Ops Ops = _dag->subgraph( vlhs ).l_op;

  // Define subsets and map of linear operations and participating variables
  _extract_linearity( Ops, vlhs );

  // Define subsets of bilinear, fractional, square and square root terms
  _target_nonlinearity( Ops );
  stats.tGRAPH += cpuclock();

  // Bi-partite graph approach to identifying reduction constraints
  switch( options.LEVEL ){
   case Options::PRIMSEQ: case Options::FULLSEQ: _search_sequential( Ops );   break;
   case Options::PRIMSIM: case Options::FULLSIM: _search_simultaneous( Ops ); break;
  }

  // Add reduced constraints to DAG
  if( add2dag ) _add_to_dag( vlhs );

  stats.tTOTAL += cpuclock();
  return _RLTRed;
}

inline void
RLTRED::_add_to_dag
( const std::vector<const FFVar*>& vlhs )
{
  // RLT did not give any reduction constraints
  if( _RLTRed.begin() == _RLTRed.end() ){
    if( options.DISPLAY >= 1 )
      std::cout << "\nNO RLT REDUCTION CONSTRAINTS FOUND\n";
    return;
  }

  // Display "raw" reduction constraints
  if( options.DISPLAY >= 2 ){
    std::cout << "\nRLT REDUCTION CONSTRAINTS:\n";
    for( auto it = _RLTRed.begin(); it != _RLTRed.end(); ++it ){
      std::cout << "  " << *(*it).first.first << "  <";
      switch( (*it).first.second ){
       case FFOp::TIMES:
        std::cout << "MULT"; break;
       case FFOp::DIV:
        std::cout << "DIV";  break;
       default:
        throw std::runtime_error( "Internal error" );
      }
      std::cout << ">  " << *(*it).second->pres << " := " << *(*it).second << std::endl;
    }
  }

  // Start timing of RLT constraint generation
  stats.tCTRDAG -= cpuclock();

  // Coefficient matrix and right-hand-side vector holding the "hidden" linear relationships
  const unsigned nRLTRed = _RLTRed.size();
  _coefRed.resize( nRLTRed, nRLTRed ); _coefRed.zero(); // Oversize and initialize to zero
  for( auto itv=_rhsRed.begin(); itv!=_rhsRed.end(); ++itv ) delete *itv;
  _rhsRed.assign( nRLTRed, 0 );
  for( auto itv=_varRed.begin(); itv!=_varRed.end(); ++itv ) delete itv->first;
  _varRed.clear(); 

  // Fill-in coefficient matrix and right-hand-side vector
  auto itRed = _RLTRed.begin();
  for( unsigned iRow = 0; itRed != _RLTRed.end(); ++itRed, iRow++ )
    _process_constraint( iRow, vlhs, itRed->first, itRed->second ); 
  //std::cout << *_dag;
  if( options.DISPLAY >= 2 ){
    std::cout << "\nRIGHT-HAND-SIDE DEPENDENTS:\n";
    for( auto itv=_rhsRed.begin(); itv!=_rhsRed.end(); ++itv )
      std::cout << "  " << *(*itv) << std::endl;
    std::cout << "\nLEFT-HAND-SIDE NEW VARIABLES:\n";
    for( auto itv=_varRed.begin(); itv!=_varRed.end(); ++itv )
      std::cout << "  " << itv->second << ": " << *itv->first << std::endl;
    std::cout << "\nCOEFFICIENT MATRIX:\n" << _coefRed;
  }

  // Perform LU decomposition of coefficient matrix
  CPPL::dgematrix L, U;
  std::vector<int> PIV;
  const int INFO = CPPL::dgetrf( _coefRed, U, L, PIV );
  if( options.DISPLAY >= 2 ){
    std::cout << "\nLU DECOMPOSITION OF COEFFICIENT MATRIX: "
              << INFO << std::endl
              << "\nMatrix U:\n" << U
              << "\nMatrix L:\n" << L
              << "\nVector PIV:\n";
    for( auto&& ip : PIV ) std::cout << ip << std::endl;
  }

  // Error during LU decomposition
  if( INFO < 0 ){
    stats.tCTRDAG += cpuclock();
    throw std::runtime_error( "Error during LU decomposition" );
  }

  // Build "hidden" relationships using forward substitution
  else if( INFO > 0 ){
    std::vector<SPolyExpr> vSCtr;
    for( auto&& ip : PIV ){
      assert(_rhsRed[ip]);
      const unsigned k = vSCtr.size();
      SPolyExpr ctr = 0;
      if( _rhsRed[ip]->cst() )
        ctr += _rhsRed[ip]->num().val();
      else if( _rhsRed[ip]->ops().first->type != FFOp::NEG )
        ctr += SPolyExpr( *_rhsRed[ip] );
      else
        ctr -= SPolyExpr( *_rhsRed[ip]->ops().first->pops[0] );
      for( unsigned j=0; j<k; j++ )
        if( !isequal( L(k,j), 0 ) ) ctr += (-L(k,j)) * vSCtr[j];
      ctr /= L(k,k);
      if( options.DISPLAY >= 2 )
        std::cout << "New constraint #" << k << ":" << ctr;
      vSCtr.push_back( ctr );
    }

    
#if defined(MC__USE_CPLEX)
    try{
      vSCtr.erase( vSCtr.begin(), vSCtr.begin()+vSCtr.size()-std::ceil(-_ILOcplex->getObjValue()) );
    }
    catch(IloException& e){
      if( options.MIPDISPLAY )
        std::cout << "Error code = " << e.getMessage() << std::endl;
      _ILPexcpt = true;
    }
#elif defined(MC__USE_GUROBI)
    try{
      vSCtr.erase( vSCtr.begin(), vSCtr.begin()+vSCtr.size()-std::ceil(-_GRBmodel->get(GRB_DoubleAttr_ObjVal)) );
    }
    catch(GRBException& e){
      if( options.MIPDISPLAY )
        std::cout << "Error code = " << e.getErrorCode() << std::endl
                  << e.getMessage() << std::endl;
      _ILPexcpt = true;
    }
#else
    vSCtr.erase( vSCtr.begin(), vSCtr.begin()+INFO-1 );
#endif

    for( auto&& sctr : vSCtr ){
      FFVar var = sctr.insert(_dag);
      auto itvar = _dag->Vars().find( &var );
      if( itvar == _dag->Vars().end() ){
        if( !isequal( var.num().val(), 0 ) ){
          std::ostringstream msg; 
          msg << "Inconsistent constraint: 0 != " << var << std::endl;
          throw std::runtime_error( msg.str() );
        }
        continue;
      }
      _ctrRed.push_back( *itvar );
      if( options.DISPLAY >= 2 ){
        std::ostringstream ostr; ostr << " OF " << var;
        _dag->output( _dag->subgraph( 1, *itvar ), ostr.str() );
      }
    }
  }

  // End timing of RLT constraint generation
  stats.tCTRDAG += cpuclock();

  // Display "hidden" relationships
  if( options.DISPLAY == 1 )
    std::cout << std::endl << _ctrRed.size() << " REDUCED RLT CONSTRAINTS\n";
  if( options.DISPLAY >= 2 ){
    std::cout << std::endl << _ctrRed.size() << " REDUCED RLT CONSTRAINTS FOUND:";
    for( unsigned i=0; i<_ctrRed.size(); i++ ) std::cout << "  " << *_ctrRed[i];
    std::cout << std::endl;
    // Generate list of operation for reduction constraints
    _dag->output( _dag->subgraph( _ctrRed ) );
  }
}

inline void
RLTRED::_process_constraint
( const unsigned iRow, const std::vector<const FFVar*>& vlhs, const t_RLTVar& RLTVar, const FFOp* linOp )
{
  _rhsRed[iRow] = new FFVar(0.);

  switch( linOp->type ){
   case FFOp::SHIFT:
   case FFOp::PLUS:
    _process_term( iRow, RLTVar.second, *RLTVar.first, *linOp->pops[0], 1. );
    _process_term( iRow, RLTVar.second, *RLTVar.first, *linOp->pops[1], 1. );
    break;

  case FFOp::NEG:
    _process_term( iRow, RLTVar.second, *RLTVar.first, *linOp->pops[0], -1. );
    break;

   case FFOp::MINUS:
    _process_term( iRow, RLTVar.second, *RLTVar.first, *linOp->pops[0], 1. );
    _process_term( iRow, RLTVar.second, *RLTVar.first, *linOp->pops[1], -1. );
    break;

   case FFOp::SCALE:
    if( linOp->pops[1]->cst() )
      _process_term( iRow, RLTVar.second, *RLTVar.first, *linOp->pops[0], linOp->pops[1]->num().val() );
    else
      _process_term( iRow, RLTVar.second, *RLTVar.first, *linOp->pops[1], linOp->pops[0]->num().val() );
    break;

   default:
     throw std::runtime_error( "Internal error" );
  }

  if( !_is_rhs( vlhs, *linOp->pres ) )
    _process_term( iRow, RLTVar.second, *RLTVar.first, *linOp->pres, -1. );
}

inline void
RLTRED::_process_term
( const unsigned iRow, const FFOp::TYPE& RLTOp, const FFVar& RLTVar,
  const FFVar& linVar, const double coef )
{
  bool foundTerm = false;
  switch( RLTOp ){
  // Multiplier variable
  case FFOp::TIMES:
    // Multiplying a constant
    if( linVar.cst() ){
      *_rhsRed[iRow] += coef * linVar.num().val() * RLTVar;
      foundTerm = true;
    }
    // In a bilinear term
    for( auto itt=_bilTerms.begin(); !foundTerm && itt!=_bilTerms.end(); ++itt ){
      if( ( *(*itt)->pops[0] == linVar && *(*itt)->pops[1] == RLTVar )
       || ( *(*itt)->pops[1] == linVar && *(*itt)->pops[0] == RLTVar ) ){
        *_rhsRed[iRow] += coef * *(*itt)->pres;
        foundTerm = true;
      }
    }
    // In a fractional term
    for( auto itt=_divTerms.begin(); !foundTerm && itt!=_divTerms.end(); ++itt ){
      if( ( *(*itt)->pres == linVar && *(*itt)->pops[1] == RLTVar )
       || ( *(*itt)->pops[1] == linVar && *(*itt)->pres == RLTVar ) ){
        *_rhsRed[iRow] += coef * *(*itt)->pops[0];
        foundTerm = true;
      }
    }
    // In an inverse term
    for( auto itt=_invTerms.begin(); !foundTerm && itt!=_invTerms.end(); ++itt ){
      if( ( *(*itt)->pres == linVar && *(*itt)->pops[1] == RLTVar )
       || ( *(*itt)->pops[1] == linVar && *(*itt)->pres == RLTVar ) ){
        *_rhsRed[iRow] += coef * (*itt)->pops[0]->num().val();
        foundTerm = true;
      }
    }
    // In a square term
    for( auto itt=_sqrTerms.begin(); !foundTerm && itt!=_sqrTerms.end(); ++itt ){
      if( *(*itt)->pops[0] == linVar && *(*itt)->pops[0] == RLTVar ){
        *_rhsRed[iRow] += coef * *(*itt)->pres;
        foundTerm = true;
      }
    }
    // In a square-root term
    for( auto itt=_sqrtTerms.begin(); !foundTerm && itt!=_sqrtTerms.end(); ++itt ){
      if( *(*itt)->pres == linVar && *(*itt)->pres == RLTVar ){
        // Add term to rhs
        *_rhsRed[iRow] += coef * *(*itt)->pops[0];
        foundTerm = true;
      }
    }
    break;

  // Divider variable
  case FFOp::DIV:
    // Dividing itself
    if( linVar == RLTVar ){
      *_rhsRed[iRow] += coef;
      foundTerm = true;
    }
    // Dividing a constant
    if( linVar.cst() ){
      for( auto itt=_invTerms.begin(); !foundTerm && itt!=_invTerms.end(); ++itt ){
        if( *(*itt)->pops[1] == RLTVar ){
          *_rhsRed[iRow] += coef * linVar.num().val() * *(*itt)->pres / (*itt)->pops[0]->num().val();;
          foundTerm = true;
        }
      }
    }
    // In a bilinear term
    for( auto itt=_bilTerms.begin(); !foundTerm && itt!=_bilTerms.end(); ++itt ){
      if( *(*itt)->pres == linVar && *(*itt)->pops[0] == RLTVar ){
        *_rhsRed[iRow] += coef * *(*itt)->pops[1];
        foundTerm = true;
      }
      else if( *(*itt)->pres == linVar && *(*itt)->pops[1] == RLTVar ){
        *_rhsRed[iRow] += coef * *(*itt)->pops[0];
        foundTerm = true;
      }
    }
    // In a fractional term
    for( auto itt=_divTerms.begin(); !foundTerm && itt!=_divTerms.end(); ++itt ){
      if( *(*itt)->pops[0] == linVar && *(*itt)->pops[1] == RLTVar ){
        *_rhsRed[iRow] += coef * *(*itt)->pres;
        foundTerm = true;
      }
      else if( *(*itt)->pops[0] == linVar && *(*itt)->pres == RLTVar ){
        *_rhsRed[iRow] += coef * *(*itt)->pops[1];
        foundTerm = true;
      }
    }
    // In an inverse term
    for( auto itt=_invTerms.begin(); !foundTerm && itt!=_invTerms.end(); ++itt ){
      if( *(*itt)->pops[0] == linVar && *(*itt)->pops[1] == RLTVar ){
        *_rhsRed[iRow] += coef * *(*itt)->pres;
        foundTerm = true;
      }
      else if( *(*itt)->pops[0] == linVar && *(*itt)->pres == RLTVar ){
        *_rhsRed[iRow] += coef * *(*itt)->pres;
        foundTerm = true;
      }
    }
    // In a square term
    for( auto itt=_sqrTerms.begin(); !foundTerm && itt!=_sqrTerms.end(); ++itt ){
      if( *(*itt)->pres == linVar && *(*itt)->pops[0] == RLTVar ){
        *_rhsRed[iRow] += coef * *(*itt)->pops[0];
        foundTerm = true;
      }
    }
    // In a square-root term
    for( auto itt=_sqrtTerms.begin(); !foundTerm && itt!=_sqrtTerms.end(); ++itt ){
      if( *(*itt)->pops[0] == linVar && *(*itt)->pres == RLTVar ){
        *_rhsRed[iRow] += coef * *(*itt)->pres;
        foundTerm = true;
      }
    }
    break;

   default:
     throw std::runtime_error( "Internal error" );
  }
  if( foundTerm ) return;

  // Insert or retreive new operation
  FFOp*newOp = new FFOp( RLTOp, const_cast<FFVar*>(&linVar), const_cast<FFVar*>(&RLTVar), 0 );
  auto itv = _varRed.insert( std::make_pair( newOp, _varRed.size() ) );
  if( !itv.second ) delete newOp;
  const unsigned jCol = itv.first->second;
  // Assign corresponding entry in coefficent matrix
  _coefRed( iRow, jCol ) = -coef;
}

inline void
RLTRED::_extract_linearity
( const typename RLTRED::l_Ops& Ops, const std::vector<const FFVar*>& vlhs )
{
  static const unsigned int nlinOp = 5;
  static const typename FFOp::TYPE linOp[nlinOp] = { FFOp::PLUS, FFOp::NEG, FFOp::MINUS, FFOp::SHIFT, FFOp::SCALE };

  // Create subset of linear operations
  _linTerms  = _subset_op( Ops, nlinOp, linOp );
  if( options.DISPLAY >= 2 ){
    std::cout << "\nLINEAR TERMS ACCOUNTED FOR IN RLT:\n";
    for( auto it = _linTerms.begin(); it != _linTerms.end(); ++it )
      std::cout << *(*it) << std::endl;
  }

  // Create subset of variables participating in linear operations
  // and map between these variables and their defining operations
  _subset_var_edge( vlhs, _linTerms, _linVars, _linEdges );
  if( options.DISPLAY >= 2 ){
    std::cout << "\nVARIABLES PARTICIPATING IN RLT TERMS:\n";
    for( auto it = _linVars.begin(); it != _linVars.end(); ++it )
      std::cout << *(*it) << std::endl;
    std::cout << "\nOPERATION <-> VARIABLE:\n";
    for( auto it = _linEdges.begin(); it != _linEdges.end(); ++it )
      std::cout << *(*it).second << "  <->  " << *(*it).first << std::endl;
  }
}

inline void
RLTRED::_target_nonlinearity
( const typename RLTRED::l_Ops& Ops )
{
  // subset of bilinear terms
  _bilTerms  = _subset_op( Ops, FFOp::TIMES );

  // subset of fractional terms
  _divTerms  = _subset_op( Ops, FFOp::DIV );

  // subset of inverse terms
  _invTerms  = _subset_op( Ops, FFOp::INV );

  // subset of square terms
  _sqrTerms  = _subset_op( Ops, FFOp::SQR );

  // subset of square root terms
  _sqrtTerms = _subset_op( Ops, FFOp::SQRT );

  // subset of square root terms
  _ipowTerms = _subset_op( Ops, FFOp::IPOW );
}

inline void
RLTRED::_search_sequential
( const typename RLTRED::l_Ops& Ops )
{
  // MAIN LOOP: Valid RRLT cut for each candidate multiplier/divider variables
  for( auto ito = Ops.rbegin(); ito != Ops.rend(); ++ito ){
    assert( (*ito)->pres );
    const pt_idVar& idVar = (*ito)->pres->id();

    // Multiplying/dividing constraints by a constant is pointless
    if( idVar.first == FFVar::CINT || idVar.first == FFVar::CREAL || (*ito)->pres->cst() ) continue;

    // Multiplier/divider variables limited to primary variables (no auxiliary)
    if( options.LEVEL == Options::PRIMSEQ && idVar.first == FFVar::AUX ) continue;

    // Multiplier/divider variables limited to nonlinear operands
    static const unsigned nlinOp = 5;
    static const typename FFOp::TYPE linOp[nlinOp] = { FFOp::PLUS, FFOp::NEG, FFOp::MINUS, FFOp::SHIFT, FFOp::SCALE };
    bool is_lin = false;
    for( auto&& op : linOp ) if( (*ito)->type == op ){ is_lin = true; break; }
    if( is_lin ) continue;

    // Drop any linear variable participating in a nonlinear terms multiplied by variable <a>idVar</a> and the corresponding edges
    auto VarMult = std::make_pair( (*ito)->pres, FFOp::TIMES );
    stats.tGRAPH -= cpuclock();
    _RLTVars = _linVars;
    _RLTEdges = _linEdges;
    _bigraph_RLT( VarMult, _RLTVars, _RLTEdges );
    stats.tGRAPH += cpuclock();

    // Identify valid reduction constraints with the candidate multiplier variable
    switch( options.METHOD ){
#if defined(MC__USE_CPLEX) || defined(MC__USE_GUROBI)
     case Options::ILP:
      stats.tILPSET -= cpuclock();
      _reset_ILPmodel();
      _append_ILPmodel( VarMult, _RLTVars, _RLTEdges );
      stats.tILPSET += cpuclock();
      stats.tCTRGEN -= cpuclock();
      _solve_ILPmodel();
      _extract_RLT();
      stats.tCTRGEN += cpuclock();
      break;
#endif
     case Options::GRAPH:
      stats.tCTRGEN -= cpuclock();
      _reduction_RLT( VarMult );
      stats.tCTRGEN += cpuclock();

      break;
    }

     if( options.NODIV ) continue;
 
   // Drop any linear variable participating in a nonlinear terms divided by variable <a>idVar</a> and the corresponding edges
    auto VarDiv = std::make_pair( (*ito)->pres, FFOp::DIV );
    _RLTVars = _linVars;
    _RLTEdges = _linEdges;
    _bigraph_RLT( VarDiv, _RLTVars, _RLTEdges );

    // Identify valid reduction constraints with the candidate divider variable
    switch( options.METHOD ){
#if defined(MC__USE_CPLEX) || defined(MC__USE_GUROBI)
     case Options::ILP:
      stats.tILPSET -= cpuclock();
      _reset_ILPmodel();
      _append_ILPmodel( VarDiv, _RLTVars, _RLTEdges );
      stats.tILPSET += cpuclock();
      stats.tCTRGEN -= cpuclock();
      _solve_ILPmodel();
      _extract_RLT();
      stats.tCTRGEN += cpuclock();
      break;
#endif
     case Options::GRAPH:
      stats.tCTRGEN -= cpuclock();
      _reduction_RLT( VarDiv );
      stats.tCTRGEN += cpuclock();
      break;
    }

  }
}

inline void
RLTRED::_search_simultaneous
( const typename RLTRED::l_Ops& Ops )
{
  if( options.METHOD == Options::GRAPH )
    throw std::runtime_error( "Simultaneous RLT reduction using graph is unavailable.\n" );

#if defined(MC__USE_CPLEX) || defined(MC__USE_GUROBI)
  // Reset ILP model of RLT reduction
  stats.tILPSET -= cpuclock();
  _reset_ILPmodel();
  stats.tILPSET += cpuclock();

  // MAIN LOOP: Valid RRLT cut for all candidate multiplier/divider variables
  for( auto ito = Ops.rbegin(); ito != Ops.rend(); ++ito ){
    assert( (*ito)->pres );
    const pt_idVar& idVar = (*ito)->pres->id();

    // Multiplying/dividing constraints by a constant is pointless
    if( idVar.first == FFVar::CINT || idVar.first == FFVar::CREAL || (*ito)->pres->cst() ) continue;

    // Multiplier/divider variables limited to primary variables (no auxiliary)
    if( options.LEVEL == Options::PRIMSIM && idVar.first == FFVar::AUX ) continue;

    // Multiplier/divider variables limited to nonlinear operands
    static const unsigned nlinOp = 5;
    static const typename FFOp::TYPE linOp[nlinOp] = { FFOp::PLUS, FFOp::NEG, FFOp::MINUS, FFOp::SHIFT, FFOp::SCALE };
    bool is_lin = false;
    for( auto&& op : linOp ) if( (*ito)->type == op ){ is_lin = true; break; }
    if( is_lin ) continue;

    // Drop any linear variable participating in a nonlinear terms multiplied by variable <a>idVar</a> and the corresponding edges
    auto VarMult = std::make_pair( (*ito)->pres, FFOp::TIMES );
    stats.tGRAPH -= cpuclock();
    _RLTVars = _linVars;
    _RLTEdges = _linEdges;
    _bigraph_RLT( VarMult, _RLTVars, _RLTEdges );
    stats.tGRAPH += cpuclock();

    // Append variables and constraints in ILP model
    stats.tILPSET -= cpuclock();
    _append_ILPmodel( VarMult, _RLTVars, _RLTEdges );
    stats.tILPSET += cpuclock();

     if( options.NODIV ) continue;
 
   // Drop any linear variable participating in a nonlinear terms divided by variable <a>idVar</a> and the corresponding edges
    auto VarDiv = std::make_pair( (*ito)->pres, FFOp::DIV );
    stats.tGRAPH -= cpuclock();
    _RLTVars = _linVars;
    _RLTEdges = _linEdges;
    _bigraph_RLT( VarDiv, _RLTVars, _RLTEdges );
    stats.tGRAPH += cpuclock();

    // Append variables and constraints in ILP model
    stats.tILPSET -= cpuclock();
    _append_ILPmodel( VarDiv, _RLTVars, _RLTEdges );
    stats.tILPSET += cpuclock();
  }

  // Identify valid reduction constraints with all candidate multiplier / divider variables
  stats.tCTRGEN -= cpuclock();
  _solve_ILPmodel();
  _extract_RLT();
  stats.tCTRGEN += cpuclock();

#else
    throw std::runtime_error( "Simultaneous RLT reduction using integer programming is unavailable.\n" );
#endif
}

inline void
RLTRED::_bigraph_RLT
( const t_RLTVar& VarRed, t_Vars& RLTVars, t_Edges& RLTEdges )
const
{
  // Erase those variables participating in nonlinear terms with variable <a>idVar</a>
  switch( VarRed.second ){
  case FFOp::TIMES:
    _bigraph_RLT_mult( VarRed.first, RLTVars, RLTEdges ); break;
  case FFOp::DIV:
    _bigraph_RLT_div( VarRed.first, RLTVars, RLTEdges );  break;
  default:
    throw std::runtime_error( "Internal error" );
  }
}

inline void
RLTRED::_bigraph_RLT_mult
( const FFVar*VarMult, t_Vars &linVars, t_Edges &linEdges )
const
{
  // Erase constants since their multiplication with <a>idMult</a> does not introduce a new nonlinear term
  for( auto it = linVars.begin(); it != linVars.end(); ){
    if( (*it)->id().first != FFVar::CINT && (*it)->id().first != FFVar::CREAL
     && !(*it)->cst() ) ++it;
    else{ linEdges.erase( *it ); it = linVars.erase( it ); }
  }

  // Erase those variables participating with <a>idMult</a> in bilinear terms
  for( auto&& term : _bilTerms ){
    if( *term->pops[0] == *VarMult ){
      linVars.erase( term->pops[1] ); linEdges.erase( term->pops[1] );
    }
    else if( *term->pops[1] == *VarMult ){
      linVars.erase( term->pops[0] ); linEdges.erase( term->pops[0] );
    }
  }

  // Erase those variables participating with <a>idMult</a> in square terms
  for( auto&& term : _sqrTerms ){
    if( *term->pops[0] == *VarMult ){
      linVars.erase( term->pops[0] ); linEdges.erase( term->pops[0] );
    }
  }

  // Erase those variables participating with <a>idMult</a> in square root terms
  for( auto it = _sqrtTerms.begin(); it != _sqrtTerms.end(); ++it ){
    if( *(*it)->pres == *VarMult ){
      linVars.erase( (*it)->pres ); linEdges.erase( (*it)->pres );
    }
  }

//  // Erase those variables participating with <a>idMult</a> in power terms
//  for( auto&& term : _ipowTerms ){
//    if( *term->pops[1] > 2 ){ // positive exponent
//      if( *term->pops[0] == *VarMult )
//       


//    if( *term->pops[0] == *VarMult ){
//      // search for 
//      linVars.erase( term->pops[1] ); linEdges.erase( term->pops[1] );
//    }
//    else if( *term->pops[1] == *VarMult ){
//      linVars.erase( term->pops[0] ); linEdges.erase( term->pops[0] );
//    }
//  }

  // Erase those variables participating with <a>idMult</a> in fractional terms
  for( auto it = _divTerms.begin(); it != _divTerms.end(); ++it ){
    if( *(*it)->pres == *VarMult ){
      linVars.erase( (*it)->pops[1] ); linEdges.erase( (*it)->pops[1] );
    }
    else if( *(*it)->pops[1] == *VarMult ){
      linVars.erase( (*it)->pres ); linEdges.erase( (*it)->pres );
    }
  }

  // Erase those variables participating with <a>idMult</a> in inverse terms
  for( auto it = _invTerms.begin(); it != _invTerms.end(); ++it ){
    if( *(*it)->pres == *VarMult ){
      linVars.erase( (*it)->pops[1] ); linEdges.erase( (*it)->pops[1] );
    }
    else if( *(*it)->pops[1] == *VarMult ){
      linVars.erase( (*it)->pres ); linEdges.erase( (*it)->pres );
    }
  }

  // INCORPORATE POW TERMS!

  if( options.DISPLAY >= 2 ){
    std::cout << "\nVARIABLES IN RLT TERMS NOT YET PARTICIPATING IN ANY PRODUCT TERM WITH "
              << VarMult->name() << ":\n";
    for( auto itv = linVars.begin(); itv != linVars.end(); ++itv )
      std::cout << *(*itv) << std::endl;
    std::cout << "\nBI-PARTITE GRAPH FOR " << VarMult->name()
              << ":\nOPERATION <-> VARIABLE\n";
    for( auto itvo = linEdges.begin(); itvo != linEdges.end(); ++itvo )
      std::cout << *(*itvo).second << "  <->  " << *(*itvo).first << std::endl;
  }
}

inline void
RLTRED::_bigraph_RLT_div
( const FFVar* VarDiv, t_Vars &linVars, t_Edges &linEdges )
const
{
  // Erase <a>idDiv</a> since division by itself yield 1
  linVars.erase( VarDiv ); linEdges.erase( VarDiv );
  //for( auto it = linVars.begin(); it != linVars.end(); ++it )
  //  if( *(*it) == *VarDiv ){ linVars.erase( *it ); break; }

  // Erase those variables participating with <a>idDiv</a> in bilinear terms
  for( auto it = _bilTerms.begin(); it != _bilTerms.end(); ++it ){
    if( *(*it)->pops[0] == *VarDiv ){
      linVars.erase( (*it)->pres ); linEdges.erase( (*it)->pres );
    }
    else if( *(*it)->pops[1] == *VarDiv ){
      linVars.erase( (*it)->pres ); linEdges.erase( (*it)->pres );
    }
  }

  // Erase those variables participating with <a>idDiv</a> in fractional terms
  for( auto it = _divTerms.begin(); it != _divTerms.end(); ++it ){
    if( *(*it)->pres == *VarDiv ){
      linVars.erase( (*it)->pops[0] ); linEdges.erase( (*it)->pops[0] );
    }
    else if( *(*it)->pops[1] == *VarDiv ){
      linVars.erase( (*it)->pops[0] ); linEdges.erase( (*it)->pops[0] );
    }
  }

  // Erase those variables participating with <a>idDiv</a> in inverse terms
  for( auto it = _invTerms.begin(); it != _invTerms.end(); ++it ){
    if( *(*it)->pops[1] == *VarDiv ){
      linVars.erase( (*it)->pres ); linEdges.erase( (*it)->pres );
    }
  }

  // Erase those variables participating with <a>idDiv</a> in square terms
  for( auto it = _sqrTerms.begin(); it != _sqrTerms.end(); ++it ){
    if( *(*it)->pops[0] == *VarDiv ){
      linVars.erase( (*it)->pres ); linEdges.erase( (*it)->pres );
    }
  }

  // Erase those variables participating with <a>idDiv</a> in square root terms
  for( auto it = _sqrtTerms.begin(); it != _sqrtTerms.end(); ++it ){
    if( *(*it)->pres == *VarDiv ){
      linVars.erase( (*it)->pops[0] ); linEdges.erase( (*it)->pops[0] );
    }
  }

  if( options.DISPLAY >= 2 ){
    std::cout << "\nVARIABLES IN RLT TERMS NOT YET PARTICIPATING IN ANY FRACTIONAL TERM WITH "
              << VarDiv->name() << ":\n";
    for( auto itv = linVars.begin(); itv != linVars.end(); ++itv )
      std::cout << *(*itv) << std::endl;
    std::cout << "\nBI-PARTITE GRAPH FOR " << VarDiv->name()
              << ":\nOPERATION <-> VARIABLE\n";
    for( auto itvo = linEdges.begin(); itvo != linEdges.end(); ++itvo )
      std::cout << *(*itvo).second << "  <->  " << *(*itvo).first << std::endl;
  }
}

#if defined(MC__USE_CPLEX)
inline typename RLTRED::t_ILPVar::iterator
RLTRED::_set_ILPvar
( const t_RLTVar& VarRed, const FFVar* pVarLin )
{
  //stats.tILPAUX -= cpuclock();
  auto itvar = _ILPvar.find( std::make_pair( VarRed, pVarLin ) );
  //stats.tILPAUX += cpuclock();
  if( itvar == _ILPvar.end() ){
    IloNumVar var( *_ILOenv, 0., 1., ILOBOOL );
    //stats.tILPAUX -= cpuclock();
    *_ILPobj += (1.+options.VARTOL)*var;
    //stats.tILPAUX += cpuclock();
    //stats.tILPAUX -= cpuclock();
    itvar = _ILPvar.insert( std::make_pair( std::make_pair( VarRed, pVarLin ), var ) ).first;
    //stats.tILPAUX += cpuclock();
  }
  return itvar;
}

inline typename RLTRED::t_ILPCtr::iterator
RLTRED::_set_ILPctr
( const t_RLTVar& VarRed, const FFOp* pOpLin )
{
  //stats.tILPAUX -= cpuclock();
  auto itctr = _ILPctr.find( std::make_pair( VarRed, pOpLin ) );
  //stats.tILPAUX += cpuclock();
  if( itctr == _ILPctr.end() ){
    IloNumVar var( *_ILOenv, 0., 1., ILOBOOL );
    //stats.tILPAUX -= cpuclock();
    *_ILPobj -= var;
    //stats.tILPAUX += cpuclock();
    //stats.tILPAUX -= cpuclock();
    itctr = _ILPctr.insert( std::make_pair( std::make_pair( VarRed, pOpLin ), var ) ).first;
    //stats.tILPAUX += cpuclock();
  }
  return itctr;
}

inline void
RLTRED::_add_ILPedge
( const IloNumVar& varbin, const IloNumVar& ctrbin )
{
  //stats.tILPAUX -= cpuclock();
  _ILOmodel->add( ctrbin >= varbin );
  //stats.tILPAUX += cpuclock();
}

#elif defined(MC__USE_GUROBI)
inline typename RLTRED::t_ILPVar::iterator
RLTRED::_set_ILPvar
( const t_RLTVar& VarRed, const FFVar* pVarLin )
{
  auto itvar = _ILPvar.find( std::make_pair( VarRed, pVarLin ) );
  if( itvar == _ILPvar.end() ){
    GRBVar var = _GRBmodel->addVar( 0., 1., 1.+options.VARTOL, GRB_BINARY );
    itvar = _ILPvar.insert( std::make_pair( std::make_pair( VarRed, pVarLin ), var ) ).first;
  }
  return itvar;
}

inline typename RLTRED::t_ILPCtr::iterator
RLTRED::_set_ILPctr
( const t_RLTVar& VarRed, const FFOp* pOpLin )
{
  auto itctr = _ILPctr.find( std::make_pair( VarRed, pOpLin ) );
  if( itctr == _ILPctr.end() ){
    GRBVar var = _GRBmodel->addVar( 0., 1., -1., GRB_BINARY );
    itctr = _ILPctr.insert( std::make_pair( std::make_pair( VarRed, pOpLin ), var ) ).first;
  }
  return itctr;
}

inline void
RLTRED::_add_ILPedge
( const GRBVar& varbin, const GRBVar& ctrbin )
{
  _GRBmodel->addConstr( GRBLinExpr( ctrbin,  1. ), GRB_GREATER_EQUAL, GRBLinExpr( varbin,  1. ) );
}
#endif

#if defined(MC__USE_CPLEX) || defined(MC__USE_GUROBI)
inline void
RLTRED::_append_ILPmodel
( const t_RLTVar& VarRed, t_Vars& RLTVars, t_Edges& RLTEdges )
{
  for( auto&& edge : RLTEdges ){
    // Append reduction operation and term in maps _ILPctr and _ILPvar
    auto itctr = _set_ILPctr( VarRed, edge.second );
    auto itvar = _set_ILPvar( VarRed, edge.first );
    // Append edge constraint in ILP model
    _add_ILPedge( itctr->second, itvar->second );
  }
}

inline void
RLTRED::_reset_ILPmodel
()
{
#if defined(MC__USE_CPLEX)
  delete _ILPobj;
  delete _ILOmodel; delete _ILOcplex;
  _ILOmodel = new IloModel( *_ILOenv );
  _ILOcplex = new IloCplex( *_ILOenv );
  _ILPobj = new IloExpr( *_ILOenv );
#elif defined(MC__USE_GUROBI)
  delete _GRBmodel;
  _GRBmodel = new GRBModel( *_GRBenv );
#endif
  _ILPvar.clear(); _ILPctr.clear();
}

inline void
RLTRED::_solve_ILPmodel
()
{
  _set_ILPoptions();
  _ILPexcpt = false;
#if defined(MC__USE_CPLEX)
  try{
    _ILOmodel->add( IloMinimize( *_ILOenv, *_ILPobj ) );
    if( options.MIPFILE != "" )
      _ILOcplex->exportModel( options.MIPFILE.c_str() );
    _ILOcplex->solve();
    if( options.MIPDISPLAY )
      std::cout << "  #reduction constraints: " << std::ceil(-_ILOcplex->getObjValue()) << std::endl;
  }
  catch(IloException& e){
    if( options.MIPDISPLAY )
      std::cout << "Error code = " << e.getMessage() << std::endl;
    _ILPexcpt = true;
  }
#elif defined(MC__USE_GUROBI)
  try{
    _GRBmodel->update();
    if( options.MIPFILE != "" )
      _GRBmodel->write( options.MIPFILE );
    fedisableexcept(FE_ALL_EXCEPT);
    _GRBmodel->optimize();
    if( options.MIPDISPLAY )
      std::cout << "  #reduction constraints: " << std::ceil(-_GRBmodel->get( GRB_DoubleAttr_ObjVal )) << std::endl;
  }
  catch(GRBException& e){
    if( options.MIPDISPLAY )
      std::cout << "Error code = " << e.getErrorCode() << std::endl
                << e.getMessage() << std::endl;
    _ILPexcpt = true;
  }
#endif
}

inline void
RLTRED::_set_ILPoptions
()
{
#if defined(MC__USE_CPLEX)
  // CPLEX options
  _ILOcplex->extract(*_ILOmodel);
  _ILOcplex->setWarning( options.MIPDISPLAY? std::cout: _ILOenv->getNullStream() );
  _ILOcplex->setOut( options.MIPDISPLAY? std::cout: _ILOenv->getNullStream() );
  _ILOcplex->setParam( IloCplex::RootAlg, options.LPALGO );
  _ILOcplex->setParam( IloCplex::EpOpt,   options.LPOPTIMTOL );
  _ILOcplex->setParam( IloCplex::EpRHS,   options.LPFEASTOL );
  _ILOcplex->setParam( IloCplex::EpGap,   options.MIPRELGAP );
  _ILOcplex->setParam( IloCplex::EpAGap,  options.MIPABSGAP );
  _ILOcplex->setParam( IloCplex::TiLim,   options.MAXCPU );
  _ILOcplex->setParam( IloCplex::PreInd,  options.LPPRESOLVE?true:false );
#elif defined(MC__USE_GUROBI)
  // Gurobi options
  _GRBmodel->getEnv().set( GRB_IntParam_OutputFlag,        options.MIPDISPLAY );
  _GRBmodel->getEnv().set( GRB_IntParam_Method,            options.LPALGO );
  _GRBmodel->getEnv().set( GRB_DoubleParam_FeasibilityTol, options.LPFEASTOL );
  _GRBmodel->getEnv().set( GRB_DoubleParam_OptimalityTol,  options.LPOPTIMTOL );
  _GRBmodel->getEnv().set( GRB_DoubleParam_MIPGap,         options.MIPRELGAP );
  _GRBmodel->getEnv().set( GRB_DoubleParam_MIPGapAbs,      options.MIPABSGAP );
  _GRBmodel->getEnv().set( GRB_DoubleParam_TimeLimit,      options.MAXCPU );
  _GRBmodel->getEnv().set( GRB_IntParam_Presolve,          options.LPPRESOLVE  );
#endif
}

inline void
RLTRED::_extract_RLT
()
{
  // Need checking successfull opimization first
  for( auto&& binctr : _ILPctr ){
#if defined(MC__USE_CPLEX)
    if( !isequal( _ILOcplex->getValue( binctr.second ), 0 ) )
#elif defined(MC__USE_GUROBI)
    if( !isequal( binctr.second.get(GRB_DoubleAttr_X), 0 ) )
#endif
      _RLTRed.insert( binctr.first );
  }
}

#endif

inline void
RLTRED::_reduction_RLT
( const t_RLTVar& VarRed )
{
  // Matching initially empty
  FFOp* NA = 0;
  _VarAssigned.assign( _RLTVars.size(), std::make_pair( NA, 0 ) );

  // Try to construct an augmenting emanating from each constraint
  auto ito = _linTerms.begin();
  for( unsigned io=0; ito != _linTerms.end(); ++ito, io++ ){
    _TermVisited.assign( _linTerms.size(), 0 );
    _VarVisited.assign( _RLTVars.size(), 0 );
    if( _augpath_RLT( std::make_pair(*ito,io) ) ) continue;

    // Append new reduction constraints
    auto jto = _linTerms.begin();
    for( unsigned jo=0; jto != _linTerms.end(); ++jto, jo++ )
      if( _TermVisited[jo] ) _RLTRed.insert( std::make_pair( VarRed, *jto ) );
  }
}

inline bool
RLTRED::_augpath_RLT
( const pt_Op& pOp )
{
  _TermVisited[pOp.second] = true;

  // Try and find immediate assignment
  auto itv = _RLTVars.begin();
  for( unsigned iv=0; itv != _RLTVars.end(); ++itv, iv++ ){
    if( _VarAssigned[iv].first ) continue;
    auto rangeOp = _RLTEdges.equal_range( *itv );
    for( auto itvo = rangeOp.first; itvo != rangeOp.second; ++itvo ){
      if( (*itvo).second != pOp.first ) continue;
      _VarAssigned[iv] = pOp;
      return true;
    }
  }

  // Try and find an augmenting path starting from another variable
  itv = _RLTVars.begin();
  for( unsigned iv=0; itv != _RLTVars.end(); ++itv, iv++ ){
    if( _VarVisited[iv] ) continue;
    auto rangeOp = _RLTEdges.equal_range( *itv );
    for( auto itvo = rangeOp.first; itvo != rangeOp.second; ++itvo ){
      if( (*itvo).second != pOp.first ) continue;
      _VarVisited[iv] = true;
      if( !_augpath_RLT( _VarAssigned[iv] ) ) continue;
      _VarAssigned[iv] = pOp;
      return true;
    }
  }

  // Failed to find an augmenting path
  return false;
}

inline typename RLTRED::t_Ops
RLTRED::_subset_op
( const typename RLTRED::l_Ops& Ops, const unsigned int nOp, const typename FFOp::TYPE*typeOp )
const
{
  // Create subset of operations of given type
  t_Ops subOps;
  for( unsigned i=0; i<nOp; i++ ){
    for( auto&& op : Ops ){
      if( op->type != typeOp[i] ) continue;
      subOps.insert( op );
    }
  }

  return subOps;
}

inline typename RLTRED::t_Ops
RLTRED::_subset_op
( const typename RLTRED::l_Ops& Ops, const typename FFOp::TYPE&typeOp )
const
{
  // Create subset of operations of given type
  t_Ops subOps;
  for( auto&& op : Ops ){
    if( op->type != typeOp ) continue;
    subOps.insert( op );
  }

  return subOps;
}

inline void
RLTRED::_subset_var_edge
( const std::vector<const FFVar*>& vlhs, const typename RLTRED::t_Ops& Ops,
  typename RLTRED::t_Vars& Vars, typename RLTRED::t_Edges& Edges )
const
{
  Vars.clear();
  Edges.clear();

  // Create subset of variables participating in a set of operations Ops
  for( auto&& operation : Ops ){
    for( auto&& operand : operation->pops ){
      Vars.insert( operand ); // Keep constants in variable list for division
      Edges.insert( std::make_pair( operand, operation ) );
    }
    auto result = operation->pres;
    if( !result || _is_rhs( vlhs, *result ) ) continue;
    Vars.insert( result ); // Keep constants in variable list for division
    Edges.insert( std::make_pair( result, operation ) );
  }
}

inline bool
RLTRED::_is_rhs
( const std::vector<const FFVar*>& vlhs, const FFVar& varRes )
const
{
  for( auto&& lhs : vlhs )
    if( *lhs == varRes ) return true;
  return false;
}

} // namespace mc

#endif
