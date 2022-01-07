// Copyright (C) 2020 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_SQUAD Decomposition of Sparse Polynomials into Quadratic Forms
\author Benoit Chachuat, Tanuj Karia & OMEGA Research Group (http://www3.imperial.ac.uk/environmentenergyoptimisation)
\date 2020
\bug No known bugs.

The class mc::SQuad defined in <tt>squad.hpp</tt> enables the reformulation of sparse multivariate polynomials into a set of quadratic forms via the introduction of auxiliary variables [<A href="https://doi.org/10.1007/BF01070233">Shor, 1987</A>; <A href="https://doi.org/10.1080/00986449208936033">Manousiouthakis & Sourlas, 1992</A>; <A href="http://www.biomedcentral.com/1752-0509/4/69">Rumschinski <i>et al.</i>, 2010</A>]

\section sec_SQUAD_algo What is the Algorithm Running Behind the Reformulation?

Consider the following sparse multivariate polynomials:
\f{align*}
p_i(\mathbf{x}) :=\ & \sum_{j=1}^{J_i} a_{i,j} m_{i,j}(\mathbf{x}),\ \ i=1,\ldots, I\\
\text{with}\ m_{i,j}(\mathbf{x}) :=\ & x_1^{\alpha_{i,j,1}}\cdots x_{n_x}^{\alpha_{i,j,n_x}},\ \ j = 1,\ldots,J_i 
\f}
with coefficients \f$a_{i,j} \in \mathbb{R}\f$ and exponents \f$\alpha_{i,j,k}\in\mathbb{Z}^+\f$. The reformulation entails the construction of a list \f$\Xi\f$ of monomials in the variables \f$\mathbf{x}\f$, with the following properties [<A href="http://www.biomedcentral.com/1752-0509/4/69">Rumschinski <i>et al.</i>, 2010</A>]:
-# For every monomial \f$m\f$ in any polynomial function \f$p_i\f$, there exists \f$\xi_1, \xi_2 \in \Xi\f$ such that \f$m = \xi_1\cdot \xi_2\f$
-# For every monomial \f$m \in \Xi\f$ of degree higher than 1, there exists \f$\xi_1, \xi_2 \in \Xi\f$ such that \f$m = \xi_1\cdot \xi_2\f$
-# The constant monomial \f$1 \in \Xi\f$

From Property 1, any polynomial function \f$p_i\f$ can be rewritten into quadratic form as:
\f{align*}
p_i(\mathbf{x},\mathbf{y}) =\ & \boldsymbol{\xi}^\intercal \mathbf{Q}_i \boldsymbol{\xi},\ i = 0,\ldots,n_{\rm t} \label{eq:mainqf}
\f}
where the variable vectors \f$\boldsymbol{\xi}\f$ correspond to the elements of the list \f$\Xi\f$; and \f$\mathbf{Q}_i\in\mathbb{S}_+^{n_\xi}\f$ are suitable real symmetric matrices with \f$n_\xi=|\Xi|\f$.

From Property 2, decomposing the elements in \f$\Xi\f$ with degree higher than 1 into products of lower-order monomials also in \f$\Xi\f$ gives rise to auxiliary quadratic constraints:
\f{align*}
\boldsymbol{\xi}^\intercal \mathbf{A}_i \boldsymbol{\xi} = 0,\ i = 0,\ldots,n_{\rm a} \label{eq:auxqf}
\f}
for suitable real symmetric matrices \f$\mathbf{A}_i\in\mathbb{S}_+^{n_\xi}\f$.

<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html Quadratization_algo.png width=35%</TD>
</TR>
</TABLE></CENTER>

The decomposition procedure implemented in mc::SQuad consists of the recursive application of Algorithm 1 to each monomial term \f$m\f$ in any polynomial functions \f$p_i\f$, starting with \f$1\f$ and any participating variables as the list \f$\Xi\f$. The monomials are processed either in increasing or decreasing order, using the <A href="https://en.wikipedia.org/wiki/Monomial_order#Graded_lexicographic_order">graded lexicographic order</a> (grlex) of monomials by default. The function <TT>Decompose</TT> proceeds by decomposing a monomial \f$m\f$ into the monomial product \f$\xi_1\cdot\xi_2\f$, and keeps processing the monomial factors \f$\xi_1\f$ and \f$\xi_2\f$ until they or their own factors are all found in the monomial set \f$\Xi\f$. The latter set is grown dynamically in the function <TT>SubExpression</TT> by appending any factor \f$\xi_1\f$ or \f$\xi_2\f$ not initially present (line 22), and so is the set \f$\Gamma\f$ of auxiliary quadratic constraints (line 21) so the quadratization is exact.

Optionally, extra quadratic equality constraints may be generated between each newly inserted monomials and other monomials in \f$\Xi\f$ (lines 23-24). These constraints arise due to multiple possible decompositions of higher-order monomials in terms of existing factors in \f$\Xi\f$ and are structurally redundant in the sense that the quadratization remains exact without them. Such redundant constraints may strengthen a QCP relaxation of a polynomial optimization program since they do not add extra auxiliary variables [<A href="https://doi.org/10.1007/s10898-006-9005-4">Liberti & Pantelides, 2006</A>; <A href="https://doi.org/10.1016/j.compchemeng.2011.01.035">Ruiz & Grossmann, 2011</A>]. However, they introduce new bilinear or square terms in the QCP relaxation which reduce its sparsity as a result.

Algorithm 1 starts by testing if \f$m\f$ can be decomposed in terms of monomials that are already present in \f$\Xi\f$ (line 2), in an attempt to reducing the number of auxiliary variables and constraints in the reformulation. Several such decompositions may be possible, in which case preference is given to the trivial decomposition \f$1\cdot m\f$ if \f$m\in\Xi\f$; otherwise, a decomposition \f$\xi_1\cdot\xi_2=m\f$ is chosen whereby \f$\xi_1\f$ and \f$\xi_2\f$ are as close as possible in the grlex order. Failing this, the algorithm tests if part of \f$m\f$ comprise an existing monomial \f$\xi_1\in\Xi\f$ with total degree greater than one, with ties broken in favor of larger monomial in grlex order (line 4). The other factor \f$\xi_2\f$ is further decomposed, as necessary (line 5). If no such decompositions are possible for the current list \f$\Xi\f$ yet \f$m\f$ is multilinear (line 6), a decomposition is chosen whereby \f$\xi_1\f$ and \f$\xi_2\f$ are as close as possible in the grlex sense (line 7), then both \f$\xi_1\f$ and \f$\xi_2\f$ are further processed. In case \f$m\f$ is not multilinear but contains one or more variables with odd degree (line 10), it is split into a multilinear term &mdash; which is further processed as earlier &mdash; and a monomial where all the variable degrees are even (line 11). This latter type of monomials are simply processed as square terms (line 15). 


\section sec_SQUAD_process How Do I Reformulate a (Set of) Sparse Polynomial Model(s)?

For illustration, consider the problem to quadratize the following polynomials:
\f{align*}
  p_0(x_0,x_1,x_2) :=\ & \left(x_0+x_1^2-2 x_2\right)^3\\
  p_1(x_1) :=\ & 2 x_1^2-1
\f}

It is convenient to leverage MC++'s sparse polynomial inclusion classes, e.g. in <A>scmodel.hpp</a>
\code
      #include "interval.hpp"
      typedef mc::Interval I;

      #include "scmodel.hpp"
      typedef mc::SCModel<I> SCM;
      typedef mc::SCVar<I> SCV;
\endcode
to help construct the polynomials in the desired format:
\code
      defined( TEST_DOXYGEN )
      SCM modSCM( 6 );
      modSCM.options.BASIS = SCM::Options::MONOM;
      unsigned const NX = 3, NP = 2;
      SCV X[NX], P[NP];
      for( unsigned i=0; i<NX; i++ )
        X[i].set( &modSCM, i, I(-1,1) );
      P[0] = pow( X[0] + sqr( X[1] ) - 2 * X[2], 3 );
      P[1] = 2 * sqr( X[1] ) - 1;

      cout << "\nSparse multivariate polynomials:";
      for( unsigned i=0; i<NP; i++ )
        cout << P[i];
\endcode
The last three line display the resulting sparse polynomial expressions, here adopting a monomial basis representation:
\verbatim
    Sparse multivariate polynomials:
    
     1.0000000e+00   3  [0]^3
    -6.0000000e+00   3  [0]^2·[2]
     1.2000000e+01   3  [0]·[2]^2
    -8.0000000e+00   3  [2]^3
     3.0000000e+00   4  [0]^2·[1]^2
    -1.2000000e+01   4  [0]·[1]^2·[2]
     1.2000000e+01   4  [1]^2·[2]^2
     3.0000000e+00   5  [0]·[1]^4
    -6.0000000e+00   5  [1]^4·[2]
     1.0000000e+00   6  [1]^6
       R     =  [ 0.0000000e+00, 0.0000000e+00]
       B     =  [-2.7000000e+01, 6.4000000e+01]
       I     =  { 0 1 2 }

    -1.0000000e+00   0  1
     2.0000000e+00   2  [1]^2
       R     =  [ 0.0000000e+00, 0.0000000e+00]
       B     =  [-1.0000000e+00, 1.0000000e+00]
       I     =  { 1 }
\endverbatim

The quadratization requires the header file <tt>squad.hpp</tt> to be included:
\code
      #include "squad.hpp"
\endcode
An environment <a>mc::SQuad</a> is defined and the sparse multivariate polynomials are quadratized by calling the method <a>mc::SQuad::process</a> for each sparse polynomial. This method takes a coefficient map (or an array of maps) as argument, which is of type <a>std::map<mc::SPolyMon,double,mc::lt_SPolyMon></a> where the class mc::SPolyMon is sued to store and manipulate sparse monomials and mc::lt_SPolyMon implements the grlex order:
\code
      SQuad SQF;
      for( unsigned i=0; i<NP; i++ )
        SQF.process( P[i].coefmon() );
      std::cout << "\nSparse quadratic forms:\n" << SQF;
\endcode
By default, the quadratization expects polynomials in monomial basis, the monomials are processed in decreased grlex order, and no redundant constraints are appended (see \ref sec_SQUAD_opt). The final line displays the equivalent quadratic forms:
\verbatim
Sparse quadratic forms:

  12 Monomials: [ 1 [0] [1] [2] [0]^2 [0]·[1] [0]·[2] [1]^2 [1]·[2] [2]^2 [1]^3 [1]^4 ]

  Quadratic form for P[0]:
    1.00000e+00    ( [0] ; [0]^2 )
   -6.00000e+00    ( [0] ; [0]·[2] )
    1.20000e+01    ( [0] ; [2]^2 )
    3.00000e+00    ( [0] ; [1]^4 )
   -8.00000e+00    ( [2] ; [2]^2 )
   -6.00000e+00    ( [2] ; [1]^4 )
    3.00000e+00    ( [0]·[1] ; [0]·[1] )
   -1.20000e+01    ( [0]·[2] ; [1]^2 )
    1.20000e+01    ( [1]·[2] ; [1]·[2] )
    1.00000e+00    ( [1]^3 ; [1]^3 )

  Quadratic form for P[1]:
   -1.00000e+00    ( 1 ; 1 )
    2.00000e+00    ( 1 ; [1]^2 )

  Auxiliary quadratic form #1:
   -1.00000e+00    ( 1 ; [1]^2 )
    1.00000e+00    ( [1] ; [1] )

  Auxiliary quadratic form #2:
   -1.00000e+00    ( 1 ; [1]^3 )
    1.00000e+00    ( [1] ; [1]^2 )

  Auxiliary quadratic form #3:
   -1.00000e+00    ( 1 ; [1]^4 )
    1.00000e+00    ( [1] ; [1]^3 )

  Auxiliary quadratic form #4:
   -1.00000e+00    ( 1 ; [1]·[2] )
    1.00000e+00    ( [1] ; [2] )

  Auxiliary quadratic form #5:
   -1.00000e+00    ( 1 ; [0]·[2] )
    1.00000e+00    ( [0] ; [2] )

  Auxiliary quadratic form #6:
   -1.00000e+00    ( 1 ; [0]·[1] )
    1.00000e+00    ( [0] ; [1] )

  Auxiliary quadratic form #7:
   -1.00000e+00    ( 1 ; [2]^2 )
    1.00000e+00    ( [2] ; [2] )

  Auxiliary quadratic form #8:
   -1.00000e+00    ( 1 ; [0]^2 )
    1.00000e+00    ( [0] ; [0] )
\endverbatim

These results show that a total of 12 monomials participate in the quadratic forms: The constant monomial \f$1\f$; the participating variables \f$x_0,x_1,x_2\f$; and 6 lifted monomials \f$x_3:=x_0^2, x_4:=x_0\cdot x_1, x_5:=x_0\cdot x_2, x_6:=x_1^2, x_7:=x_1\cdot x_2, x_8:=x_2^2, x_9:=x_1^3, x_10:=x_1^4\f$. A reformulation of the factorable function \f${\bf f}\f$ in terms of these monomials is given by:
\f{align*}
  q_0(x_0,\ldots,x_10) = x_0 x_3 - 6 x_0 x_5 + 12 x_0 x_8 + 3 x_0 x_{10} - 8 x_2 x_8 - 6 x_2 x_{10} + 3 x_4^2 - 12 x_5 x_6 + 12 x_7^2 + x_9^2\\
  q_1(x_0,\ldots,x_10) = - 1 + 2 x_1^2
\f}
with the following 8 auxiliary quadratic constraints are also generated:
\f{align*}
  x_0^2 - x_3 =\ & 0\\
  x_0x_1 - x_4 =\ & 0\\
  x_0x_2 - x_5 =\ & 0\\
  x_1^2 - x_6 =\ & 0\\
  x_1x_2 - x_7 =\ & 0\\
  x_2^2 - x_8 =\ & 0\\
  x_1x_6 - x_9 =\ & 0\\
  x_1x_9 - x_{10} =\ & 0\\

\f}

With the redundant constraint option activated before the quadratization:
\code
      SQuad SQF;
      SQuad::options.REDUC = true;
      for( unsigned i=0; i<NP; i++ )
        SQF.process( P[i].coefmon() );
      std::cout << "\nSparse quadratic forms:\n" << SQF;
\endcode
the final display changes to:
\verbatim
Sparse quadratic forms:

  12 Monomials: [ 1 [0] [1] [2] [0]^2 [0]·[1] [0]·[2] [1]^2 [1]·[2] [2]^2 [1]^3 [1]^4 ]

  Quadratic form for P[0]:
    1.00000e+00    ( [0] ; [0]^2 )
   -6.00000e+00    ( [0] ; [0]·[2] )
    1.20000e+01    ( [0] ; [2]^2 )
    3.00000e+00    ( [0] ; [1]^4 )
   -8.00000e+00    ( [2] ; [2]^2 )
   -6.00000e+00    ( [2] ; [1]^4 )
    3.00000e+00    ( [0]·[1] ; [0]·[1] )
   -1.20000e+01    ( [0]·[2] ; [1]^2 )
    1.20000e+01    ( [1]·[2] ; [1]·[2] )
    1.00000e+00    ( [1]^3 ; [1]^3 )

  Quadratic form for P[1]:
   -1.00000e+00    ( 1 ; 1 )
    2.00000e+00    ( 1 ; [1]^2 )

  Auxiliary quadratic form #1:
   -1.00000e+00    ( 1 ; [1]^2 )
    1.00000e+00    ( [1] ; [1] )

  Auxiliary quadratic form #2:
   -1.00000e+00    ( 1 ; [1]^3 )
    1.00000e+00    ( [1] ; [1]^2 )

  Auxiliary quadratic form #3:
   -1.00000e+00    ( [1] ; [1]^3 )
    1.00000e+00    ( [1]^2 ; [1]^2 )

  Auxiliary quadratic form #4:
   -1.00000e+00    ( 1 ; [1]^4 )
    1.00000e+00    ( [1] ; [1]^3 )

  Auxiliary quadratic form #5:
   -1.00000e+00    ( [1] ; [1]^4 )
    1.00000e+00    ( [1]^2 ; [1]^3 )

  Auxiliary quadratic form #6:
   -1.00000e+00    ( [1]^2 ; [1]^4 )
    1.00000e+00    ( [1]^3 ; [1]^3 )

  Auxiliary quadratic form #7:
   -1.00000e+00    ( 1 ; [1]·[2] )
    1.00000e+00    ( [1] ; [2] )

  Auxiliary quadratic form #8:
   -1.00000e+00    ( [1] ; [1]·[2] )
    1.00000e+00    ( [2] ; [1]^2 )

  Auxiliary quadratic form #9:
    1.00000e+00    ( [2] ; [1]^3 )
   -1.00000e+00    ( [1]^2 ; [1]·[2] )

  Auxiliary quadratic form #10:
    1.00000e+00    ( [2] ; [1]^4 )
   -1.00000e+00    ( [1]·[2] ; [1]^3 )

  Auxiliary quadratic form #11:
   -1.00000e+00    ( 1 ; [0]·[2] )
    1.00000e+00    ( [0] ; [2] )

  Auxiliary quadratic form #12:
    1.00000e+00    ( [0] ; [1]·[2] )
   -1.00000e+00    ( [1] ; [0]·[2] )

  Auxiliary quadratic form #13:
   -1.00000e+00    ( 1 ; [0]·[1] )
    1.00000e+00    ( [0] ; [1] )

  Auxiliary quadratic form #14:
    1.00000e+00    ( [0] ; [1]^2 )
   -1.00000e+00    ( [1] ; [0]·[1] )

  Auxiliary quadratic form #15:
    1.00000e+00    ( [1] ; [0]·[2] )
   -1.00000e+00    ( [2] ; [0]·[1] )

  Auxiliary quadratic form #16:
    1.00000e+00    ( [0] ; [1]·[2] )
   -1.00000e+00    ( [2] ; [0]·[1] )

  Auxiliary quadratic form #17:
    1.00000e+00    ( [0] ; [1]^3 )
   -1.00000e+00    ( [0]·[1] ; [1]^2 )

  Auxiliary quadratic form #18:
    1.00000e+00    ( [0] ; [1]^4 )
   -1.00000e+00    ( [0]·[1] ; [1]^3 )

  Auxiliary quadratic form #19:
   -1.00000e+00    ( 1 ; [2]^2 )
    1.00000e+00    ( [2] ; [2] )

  Auxiliary quadratic form #20:
   -1.00000e+00    ( [0] ; [2]^2 )
    1.00000e+00    ( [2] ; [0]·[2] )

  Auxiliary quadratic form #21:
   -1.00000e+00    ( [1] ; [2]^2 )
    1.00000e+00    ( [2] ; [1]·[2] )

  Auxiliary quadratic form #22:
   -1.00000e+00    ( [0]·[1] ; [2]^2 )
    1.00000e+00    ( [0]·[2] ; [1]·[2] )

  Auxiliary quadratic form #23:
   -1.00000e+00    ( [1]^2 ; [2]^2 )
    1.00000e+00    ( [1]·[2] ; [1]·[2] )

  Auxiliary quadratic form #24:
   -1.00000e+00    ( 1 ; [0]^2 )
    1.00000e+00    ( [0] ; [0] )

  Auxiliary quadratic form #25:
    1.00000e+00    ( [0] ; [0]·[1] )
   -1.00000e+00    ( [1] ; [0]^2 )

  Auxiliary quadratic form #26:
    1.00000e+00    ( [0] ; [0]·[2] )
   -1.00000e+00    ( [2] ; [0]^2 )
\endverbatim
The quadratization still introduces 12 monomials but 26 auxiliary quadratic constraints are now identified, including the following 18 redundant quadratic constraints:
\f{align*}
  x_6^2 - x_1x_9 =\ & 0\\
  x_6x_9 - x_1x_{10} =\ & 0\\
  x_9^2 - x_6x_{10} =\ & 0\\
  x_1x_7 - x_2x_6 =\ & 0\\
  x_2x_9 - x_6x_7 =\ & 0\\
  x_2x_{10} - x_7x_9 =\ & 0\\
  x_0x_7 - x_1x_5 =\ & 0\\
  x_0x_6 - x_1x_4 =\ & 0\\
  x_1x_5 - x_2x_4 =\ & 0\\
  x_0x_7 - x_2x_4 =\ & 0\\
  x_4x_6 - x_0x_9 =\ & 0\\
  x_4x_9 - x_0x_{10} =\ & 0\\
  x_2x_5 - x_0x_8 =\ & 0\\
  x_2x_7 - x_1x_8 =\ & 0\\
  x_5x_7 - x_4x_8 =\ & 0\\
  x_7^2 - x_6x_8 =\ & 0\\
  x_0x_4 - x_1x_3 =\ & 0\\
  x_0x_5 - x_2x_3 =\ & 0\\
\f}

The quadratization may also be performed in Chebyshev basis instead of the default monomial basis. And the monomials may be processed in increasing grlex order rather than the default decreasing grlex order.

Pointers of type mc::SPolyMon to the monomials participating in the quadratic forms can be retrieved with the method mc::SQuad::SetMon, which is of type std::set<SPolyMon const*,lt_pSPolyMon>.

The sparse coefficient matrices defining the main quadratic forms for each processed multivariate polynomial and the corresponding auxiliary quadratic constraints can be retrieved with the methods mc::SQuad::MatFct and mc::SQuad::MatRed, respectively. These coefficient matrices are of the type std::map<std::pair<SPolyMon const*,SPolyMon const*>,double,lt_SQuad> where the comparison operator mc::lt_SQuad orders monomial pairs in grlex order for both elements. 


\section sec_SQUAD_opt What Are the Options in mc::SQuad and How Are They Set?

The public static class member mc::SQuad::options that can be used to set/modify the options; e.g.,

\code
      mc::SQuad::options.BASIS = mc::SQuad::Options::CHEB;
      mc::SQuad::options.ORDER = mc::SQuad::Options::DEC;
\endcode

The available options are the following:

<TABLE border="1">
 <TR><TH><b>Name</b>  <TD><b>Type</b> <TD><b>Default</b> <TD><b>Description</b>
 <TR><TH><tt>mc::SQuad::Options::BASIS</tt> <TD><tt>int</tt> <TD>mc::SQuad::Options::MONOM <TD>Basis representation of the multivariate polynomial - mc::SQuad::Options::MONOM: monomial basis, mc::SQuad::Options::CHEB: Chebyshev basis
 <TR><TH><tt>mc::SQuad::Options::ORDER</tt> <TD><tt>int</tt> <TD>mc::SQuad::Options::DEC <TD>Processing order for the monomial terms - mc::SQuad::Options::INC: increasing grlex monomial order, mc::SQuad::Options::DEC: decreasing grlex monomial order
 <TR><TH><tt>mc::SQuad::Options::REDUC</tt> <TD><tt>bool</tt> <TD>false <TD>Whether to search for and append extra reduction constraints
 <TR><TH><tt>mc::SQuad::Options::DISPLEN</tt> <TD><tt>unsigned int</tt> <TD>5 <TD>Number of digits in output stream
</TABLE>


\section sec_SQUAD_err What Errors Can Be Encountered during the Quadratization of a Multivariate Polynomial?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type mc::SQuad::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during a quadratization, and then make the appropriate changes. Should an exception be thrown and not caught, the program will stop.

Possible errors encountered during quadratization of a multivariate polynomial are:

<TABLE border="1">
 <TR><TH><b>Number</b> <TD><b>Description</b>
 <TR><TH><tt>-33</tt> <TD>Internal error
</TABLE>

\section sec_SQUAD_refs References

- L Liberti, CC Pantelides, <A href="https://doi.org/10.1007/s10898-006-9005-4">An exact reformulation algorithm for large nonconvex NLPs involving bilinear terms</A>, <I>Journal of Global Optimization</I> <B>36</B>:161–189, 2006
- V Manousiouthakis, D Sourlas, <A href="https://doi.org/10.1080/00986449208936033">A global optimization approach to rationally constrained rational programming</A>, <I>Chemical Engineering Communications</I> <B>115</B>:127–147, 1992
- JP Ruiz, IE Grossmann, <A href="https://doi.org/10.1016/j.compchemeng.2011.01.035">Using redundancy to strengthen the relaxation for the global optimization of MINLP problems</A>, <I>Computers & Chemical Engineering</I> <B>35</B>:2729–2740, 2011
- P Rumschinski, S Borchers, S Bosio, R Weismantel, R Findeisen, <A href="http://www.biomedcentral.com/1752-0509/4/69">Set-base dynamical parameter estimation and model invalidation for biochemical reaction networks</A>, <I>BMC Systems Biology</I> <b>4</b>:69, 2010
- NZ Shor, <A href="https://doi.org/10.1007/BF01070233">Class of global minimum bounds of polynomial functions</A>, <I>Cybernetics</I> <B>23</B>:731–734, 1987 
- <A href="https://en.wikipedia.org/w/index.php?title=Sum-of-squares_optimization&oldid=929488354">Sum-of-squares optimization</A>, <i>Wikipedia, accessed: 27-Jan-2020
.
*/

#ifndef MC__SQUAD_H
#define MC__SQUAD_H

#include <list>
#include "spolymon.hpp"
#include "spolyexpr.hpp"

#define MC__SQUAD_CHECK
#undef  MC__SQUAD_PROCESS_DEBUG

namespace mc
{
//! @brief C++ structure for ordering of monomial pairs
struct lt_SQuad
{
  // Comparison operator
  bool operator
    ()
    ( std::pair< SPolyMon const*, SPolyMon const* > const& pMon1,
      std::pair< SPolyMon const*, SPolyMon const* > const& pMon2 )
    const
    {
      // Order based on first monomial first
      if( lt_SPolyMon()( *pMon1.first, *pMon2.first ) ) return true;
      if( lt_SPolyMon()( *pMon2.first, *pMon1.first ) ) return false;
      // Order based on second monomial next
      if( lt_SPolyMon()( *pMon1.second, *pMon2.second ) ) return true;
      if( lt_SPolyMon()( *pMon2.second, *pMon1.second ) ) return false;
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
class SQuad
////////////////////////////////////////////////////////////////////////
{
  friend  std::ostream& operator<< ( std::ostream&, SQuad const& );

public:

  typedef std::set< unsigned > t_SPolyVar;
  typedef std::set< SPolyMon, lt_SPolyMon > t_SPolyMon;
  typedef std::set< SPolyMon const*, lt_pSPolyMon > t_pSPolyMon;
  typedef std::map< SPolyMon, double, lt_SPolyMon > t_SPolyMonCoef;
  typedef std::pair< SPolyMon const*, SPolyMon const* > key_SQuad;
  typedef std::map< key_SQuad, double, lt_SQuad > t_SQuad;
  
  //! @brief Options of mc::SQuad
  static struct Options
  {
    //! @brief Constructor
    Options():
      BASIS(MONOM), ORDER(DEC), REDUC(false), DISPLEN(5)
      {}
    //! @brief Assignment of mc::SQuad::Options
    Options& operator=
      ( Options& opt ){
        BASIS   = opt.BASIS;
        ORDER   = opt.ORDER;
        REDUC   = opt.REDUC;
        DISPLEN = opt.DISPLEN;
        return *this;
      }
    //! @brief Available basis representations
    enum BASIS_TYPE{
      MONOM=0,	//!< Monomial basis
      CHEB	//!< Chebyshev basis
    };
    //! @brief Available processing order
    enum ORDER_TYPE{
      INC=0,	//!< By increasing order
      DEC	//!< By decreasing order
    };
    //! @brief Basis representation of the quadratic form
    int BASIS;
    //! @brief Processing order for the monomial terms
    int ORDER;
    //! @brief Whether to search for and append extra reduction constraints
    bool REDUC;
    //! @brief Number of digits in output stream for sparse polynomial coefficients
    unsigned DISPLEN;
  } options;

  //! @brief Default Constructor
  SQuad
    ()
    {}

  //! @brief Destructor
  virtual ~SQuad
    ()
    { _reset(); }
  
  //! @brief Process the sparse polynomials in array <a>pPol</a> indexed by <a>ndxSPol</a>
  double process
    ( std::set<unsigned> const& ndxSPol, SQuad::t_SPolyMonCoef const* pSPol,
      int const BASIS=options.BASIS, bool const CHECK=false );

  //! @brief Process the <a>nPol</a> sparse polynomials in array <a>pPol</a>
  double process
    ( unsigned const nSPol, SQuad::t_SPolyMonCoef const* pSPol,
      int const BASIS=options.BASIS, bool const CHECK=false );

  //! @brief Process the sparse polynomial <a>Pol</a>
  double process
    ( SQuad::t_SPolyMonCoef const& SPol, int const BASIS=options.BASIS,
      bool const CHECK=false );

  //! @brief Decompose the quadratic expression <a>mat</a> into separable expressions
  std::list< t_SQuad > separate
    ( t_SQuad const& mat )
    const;

  //! @brief Factorize the quadratic expression <a>mat</a> using eigenvalue decomposition
  std::multimap< double, t_SPolyMonCoef > factorize
    ( SQuad::t_SQuad const& mat )
    const;

  //! @brief Generate positive semi-definite cuts to tighten the quadratic reformulation
  void tighten
    ( bool const threevar=false );
//    ();

  //! @brief Retreive reference to vector of sparse coefficient matrices defining the main quadratic forms
  std::vector<t_SQuad> const& MatFct
    ()
    const
    { return _MatFct; }

  //! @brief Retreive reference to vector of sparse coefficient matrices defining the auxiliary quadratic forms
  std::vector<t_SQuad> const& MatRed
    ()
    const
    { return _MatRed; }

  //! @brief Retreive reference to vector of sparse coefficient matrices defining the positive semi-definite cuts
  std::vector<t_SQuad> const& MatPSD
    ()
    const
    { return _MatPSD; }

  //! @brief Retreive reference to set of monomials in quadratic forms
  t_SPolyMon const& SetMon
    ()
    const
    { return _SetMon; }

  //! @brief Reset quadratic form expressions
  void reset
    ()
    { _reset(); }
    
  //! @brief Exceptions of mc::SQuad
  class Exceptions
  {
   public:
    //! @brief Enumeration type for SQuad exception handling
    enum TYPE{
      INTERNAL = -33  //!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
       case INTERNAL:
       default:
        return "mc::SQuad\t Internal error";
      }
    }
   private:
    TYPE _ierr;
  };

protected:

  //! @brief Set of monomials in quadratic forms
  t_SPolyMon _SetMon;

  //! @brief Vector of sparse coefficient matrices defining the main quadratic forms
  std::vector<t_SQuad> _MatFct;

  //! @brief Vector of sparse coefficient matrices defining the auxiliary quadratic forms
  std::vector<t_SQuad> _MatRed;

  //! @brief Vector of sparse coefficient matrices defining the postiive semi-defninite cuts
  std::vector<t_SQuad> _MatPSD;

  //! @brief Reorder entries in a monomial pair
  std::pair< SPolyMon const*, SPolyMon const* > _reorder
    ( std::pair< SPolyMon const*, SPolyMon const* > pMon )
    const;

  //! @brief Create set of monomials for Chebyshev product between <a>mon1</a> and <a>mon2</a>
  std::set< SPolyMon, lt_SPolyMon > _prodmon
    ( SPolyMon const& mon1, SPolyMon const& mon2 )
    const;

  //! @brief Insert new entry in matrix <a>mat</a> corresponding to the mononial <a>pMon</a> in a product with the constant monomial and with corresponding coefficient <a>coef</a>
  bool _insert
    ( t_SQuad& mat, SPolyMon const* pMon, double const coef, bool const add=false );//, bool const rem=false );

  //! @brief Insert new entry in matrix <a>mat</a> corresponding to the product between mononials <a>pLMon</a> and <a>pRMon</a> with corresponding coefficient <a>coef</a>
  bool _insert
    ( t_SQuad& mat, SPolyMon const* pLMon, SPolyMon const* pRMon,
      double const coef, bool const add=false );

  //! @brief Insert new entry in matrix <a>mat</a> corresponding to the product between mononials <a>pLMon</a> and <a>pRMon</a> with corresponding coefficient <a>coef</a>, and insert the associated low-order terms in monomial map <a>coefmon</a>
  bool _insert
    ( t_SQuad& mat, std::map< SPolyMon, double, lt_SPolyMon >& coefmon, 
      SPolyMon const* pLMon, SPolyMon const* pRMon, double const coef,
      bool const add=false );
      
  //! @brief Decompose <a>mon</a> into product of two monomials in <a>_SetMon</a> by expanding <a>_SetMon</a> as necessary
  std::pair< SPolyMon const*, SPolyMon const* > _decompose
    ( SPolyMon const& mon );
    
  //! @brief Search for <a>mon</a> in <a>_SetMon</a> and append it to <a>_SetMon</a> if missing along reduction constraints in <a>_MatRed</a>
  typename t_SPolyMon::iterator _subexpression
    ( SPolyMon const& mon );

  //! @brief Search for extra reduction constraints for <a>mon</a> and append them to <a>_MatRed</a>
  void _reduction
    ( SPolyMon const& mon );
    
  //! @brief Check quadratic form expressions by comparing with Chebyshev model
  //template <typename T>
  double check
    ( unsigned const nSPol, SQuad::t_SPolyMonCoef const* pSPol, int const BASIS )
    const;

  //! @brief Check if a term exist in current quadratic forms 
  bool _find
    ( SPolyMon const* pMon1, SPolyMon const* pMon2 ) 
    const;

  //! @brief Check possible decompositions using existing monomial in _SetMon
  void _candidates
    ( std::set< std::pair< SPolyMon const*, SPolyMon const* >, lt_SQuad >& CandidateDec,
      SPolyMon const& mon )
    const;

  //! @brief Check connections between a quadratic entry and a quadratic form
  bool _isconnected
    ( std::pair< SPolyMon const*, SPolyMon const* > const& entry, 
      t_SQuad const& mat )
    const;
    
  //! @brief Append participating variables in polynomial array pSPol into SVar - return value indicates whether or not there are any participating variables
  static bool _deps
    ( unsigned const nSPol, SQuad::t_SPolyMonCoef const* pSPol, SQuad::t_SPolyVar& SVar );

  //! @brief Append participating variables in polynomial SPol into SVar - return value indicates whether or not there are any participating variables
  static bool _deps
    ( SQuad::t_SPolyMonCoef const& SPol, SQuad::t_SPolyVar& SVar );

  //! @brief Express polynomial SPol as a univariate polynomial in variable ivar
  static void _unipol
    ( unsigned const ivar, SQuad::t_SPolyMonCoef const& SPol,
      std::vector< SQuad::t_SPolyMonCoef >& veccoef );

  //! @brief Determine the maximal order of variable ivar in polynomial SPol - return value is the maximum order
  static unsigned _maxord
    ( unsigned const ivar, SQuad::t_SPolyMonCoef const& SPol );

  //! @brief Convert univariate polynomial from monomial to Chebyshev basis
  static void _pow2cheb
    ( std::vector< SQuad::t_SPolyMonCoef >& veccoef );
    
  //! @brief Convert univariate polynomial from Chebyshev to monomial basis
  static void _cheb2pow
    ( std::vector< SQuad::t_SPolyMonCoef >& veccoef );

  //! @brief Convert univariate polynomial from Chebyshev to monomial basis
  static void _convert
    ( SQuad::t_SPolyMonCoef& SPol, SQuad::t_SPolyVar const& SVar,
      int const BASIS );

private:

  //! @brief Reset the quadratic form
  void _reset
    ();
};

inline SQuad::Options SQuad::options;

////////////////////////////////////////////////////////////////////////

inline std::ostream&
operator<<
( std::ostream& out, SQuad const& quad )
{
  const int BASIS = quad.options.BASIS;
  const unsigned DISPLEN = quad.options.DISPLEN;
  out << std::scientific << std::setprecision(DISPLEN)
      << std::right;

  // Output set of monomials
  out << std::endl << "  " << quad._SetMon.size() << " Monomials: [ ";
  for( auto const& mon : quad._SetMon )
    out << mon.display( BASIS ) << " ";
  out << "]" << std::endl;

  // Output quadratic forms
  unsigned count = 0;
  for( auto const& mat : quad._MatFct ){
    out << std::endl << "  Quadratic form for P[" << count++ << "]:" << std::endl;
    for( auto const& term : mat )
      out << "   " << std::right << std::setw(DISPLEN+7) << term.second << "   "
          << " ( " << term.first.first->display( BASIS )
          << " ; " << term.first.second->display( BASIS ) << " )"
          << std::endl;
  }
  count = 0;
  for( auto const& mat : quad._MatRed ){
    out << std::endl << "  Auxiliary quadratic form #" << ++count << ":" << std::endl;
    for( auto const& term : mat )
      out << "   " << std::right << std::setw(DISPLEN+7) << term.second << "   "
          << " ( " << term.first.first->display( BASIS )
          << " ; " << term.first.second->display( BASIS ) << " )"
          << std::endl;
  }
  count = 0;
  for( auto const& mat : quad._MatPSD ){
    out << std::endl << "  Positive semi-definite cut #" << ++count << ":" << std::endl;
    for( auto const& term : mat )
      out << "   " << std::right << std::setw(DISPLEN+7) << term.second << "   "
          << " ( " << term.first.first->display( BASIS )
          << " ; " << term.first.second->display( BASIS ) << " )"
          << std::endl;
  }
  return out;
}

inline bool
SQuad::_find
( SPolyMon const* pMon1, SPolyMon const* pMon2 ) 
const
{
  for( auto const& mat : _MatFct )
    if( mat.count( std::make_pair( pMon1, pMon2 ) ) ) return true;
  for( auto const& mat : _MatRed )
    if( mat.count( std::make_pair( pMon1, pMon2 ) ) ) return true;
  return false;
}


//inline void
//SQuad::tighten
//( unsigned const maxsize )
//{
//  _MatPSD.clear();
//  std::vector< SPolyMon const* > psdmon;
//  auto itmon = _SetMon.cbegin();
//  for( ; itmon != _SetMon.cend(); ++itmon ){

//    // check square monomial *itmon1 participating
//    if( itmon->tord && !_find( &*itmon, &*itmon ) ) continue;  
//    psdmon.assign( 1, &*itmon );
//    
//    // find next compatible monomial
//    _tighten( maxsize, psdmon, ++itmon );
//  }
//}


//inline void
//SQuad::_tighten
//( unsigned const maxsize, std::vector< SPolyMon const* >& psdmon,
//  t_SPolyMon::iterator itmon, std::vector<t_SQuad>::iterator itpsd )
//{
//  for( ; itmon != _SetMon.cend(); ++itmon ){

//    // check square monomial *itmon participating
//    if( !_find( &*itmon, &*itmon ) continue;

//    // check cross-term with *itmon participating
//    bool candidate = true;
//    for( auto const& pmon : psdmon ){
//      if( _find( pmon, &*itmon ) continue;
//      candidate = false;
//      break;
//    }
//    if( !candidate ) continue;

//    // add positive-definite cuts for current level
//    _MatPSD.push_back( t_SQuad() );
//    auto& mat1 = _MatPSD.back();
//    if( psdmon.size() == 1 )
//      assert( _insert( mat1, psdmon.back(), psdmon.back(), 1. ) ); 
//    assert( _insert( mat1, &*itmon, &*itmon,  1. ) );
//    for( auto const& pmon : psdmon )
//      assert( _insert( mat1, pmon, &*itmon, -2. ) );

//      if( options.BASIS == Options::MONOM && !itmon1->gcexp()%2 && !itmon2->gcexp()%2 )
//        continue;
//      _MatPSD.push_back( t_SQuad() );
//      auto& mat2 = _MatPSD.back();
//      assert( _insert( mat2, &*itmon1, &*itmon1,  1. )
//           && _insert( mat2, &*itmon2, &*itmon2,  1. )
//           && _insert( mat2, &*itmon1, &*itmon2,  2. ) );
//    }
//  }
//}


inline void
SQuad::tighten
( bool const threevar )
{
  //std::cout << "ENTER TIGHTEN\n";
  _MatPSD.clear();
  auto itmon1 = _SetMon.cbegin();
  for( ; itmon1 != _SetMon.cend(); ++itmon1 ){

    // check square monomial *itmon1 participating
    if( itmon1->tord && !_find( &*itmon1, &*itmon1 ) ) continue;
    auto itmon2 = itmon1;
    for( ++itmon2; itmon2 != _SetMon.cend(); ++itmon2 ){

      // check square monomial *itmon2 and cross-term *itmon1.*itmon2 participating
      if( (itmon1->tord && !_find( &*itmon1, &*itmon2 )) || !_find( &*itmon2, &*itmon2 ) ) continue;

      // add positive-definite cuts
      unsigned pos = _MatPSD.size();
      _MatPSD.push_back( t_SQuad() );
      auto& mat1 = _MatPSD.back();
      assert( _insert( mat1, &*itmon1, &*itmon1,  1. )
           && _insert( mat1, &*itmon2, &*itmon2,  1. )
           && _insert( mat1, &*itmon1, &*itmon2, -2. ) );
      if( !threevar && options.BASIS == Options::MONOM
       && !itmon1->gcexp()%2 && !itmon2->gcexp()%2 ) continue;
      _MatPSD.push_back( t_SQuad() );
      auto& mat2 = _MatPSD.back();
      assert( _insert( mat2, &*itmon1, &*itmon1,  1. )
           && _insert( mat2, &*itmon2, &*itmon2,  1. )
           && _insert( mat2, &*itmon1, &*itmon2,  2. ) );

      if( !threevar ) continue;
      auto itmon3 = itmon2;
      for( ++itmon3; itmon3 != _SetMon.cend(); ++itmon3 ){

        // check square monomial *itmon2 and cross-term *itmon1.*itmon2 participating
        if( (itmon1->tord && !_find( &*itmon1, &*itmon3 )) || !_find( &*itmon2, &*itmon3 )
         || !_find( &*itmon3, &*itmon3 ) ) continue;

        // add positive-definite cuts
        _MatPSD.push_back( _MatPSD[pos] );
        auto& mat3 = _MatPSD.back();
        assert( _insert( mat3, &*itmon3, &*itmon3,  1. )
             && _insert( mat3, &*itmon1, &*itmon3, -2. )
             && _insert( mat3, &*itmon2, &*itmon3,  2. ) );
        _MatPSD.push_back( _MatPSD[pos] );
        auto& mat4 = _MatPSD.back();
        assert( _insert( mat4, &*itmon3, &*itmon3,  1. )
             && _insert( mat4, &*itmon1, &*itmon3,  2. )
             && _insert( mat4, &*itmon2, &*itmon3, -2. ) );
        _MatPSD.push_back( _MatPSD[pos+1] );
        auto& mat5 = _MatPSD.back();
        assert( _insert( mat5, &*itmon3, &*itmon3,  1. )
             && _insert( mat5, &*itmon1, &*itmon3, -2. )
             && _insert( mat5, &*itmon2, &*itmon3, -2. ) );
        if( options.BASIS == Options::MONOM && !itmon1->gcexp()%2
         && !itmon2->gcexp()%2 && !itmon3->gcexp()%2 )
          continue;
        _MatPSD.push_back( _MatPSD[pos+1] );
        auto& mat6 = _MatPSD.back();
        assert( _insert( mat6, &*itmon3, &*itmon3,  1. )
             && _insert( mat6, &*itmon1, &*itmon3,  2. )
             && _insert( mat6, &*itmon2, &*itmon3,  2. ) );
      }
    }
  }
}

//inline void
//SQuad::tighten
//()
//{
//  _MatPSD.clear();
//  std::list< std::pair< SPolyMon const*, SPolyMon const* > > psdlist;
//  auto itmon1 = _SetMon.cbegin();
//  for( ; itmon1 != _SetMon.cend(); ++itmon1 ){

//    // check square monomial *itmon1 participating
//    if( itmon1->tord && !_find( &*itmon1, &*itmon1 ) ) continue;
//    auto itmon2 = itmon1;
//    for( ++itmon2; itmon2 != _SetMon.cend(); ++itmon2 ){

//      // check square monomial *itmon2 and cross-term *itmon1.*itmon2 participating
//      if( !_find( &*itmon1, &*itmon2 ) || !_find( &*itmon2, &*itmon2 ) ) continue;

//      // add positive-definite cuts
//      _MatPSD.push_back( t_SQuad() );
//      auto& mat1 = _MatPSD.back();
//      assert( _insert( mat1, &*itmon1, &*itmon1,  1. )
//           && _insert( mat1, &*itmon2, &*itmon2,  1. )
//           && _insert( mat1, &*itmon1, &*itmon2, -2. ) );
//      if( options.BASIS == Options::MONOM && !itmon1->gcexp()%2 && !itmon2->gcexp()%2 )
//        continue;
//      _MatPSD.push_back( t_SQuad() );
//      auto& mat2 = _MatPSD.back();
//      assert( _insert( mat2, &*itmon1, &*itmon1,  1. )
//           && _insert( mat2, &*itmon2, &*itmon2,  1. )
//           && _insert( mat2, &*itmon1, &*itmon2,  2. ) );
//    }
//  }
//}


inline std::list< SQuad::t_SQuad >
SQuad::separate
( SQuad::t_SQuad const& mat )
const
{
  std::list< t_SQuad > submat;
  for( auto itmat : mat ){
    std::pair< SPolyMon const*, SPolyMon const* > ijmon = itmat.first;
    double coef = itmat.second;
    auto itsubmat = submat.end();
//    // Search for alternative decompositions for terms multiplying monomial '1'
//    // NEED CORRECTION IN ORDER TO ACCOUNT THAT DECOMPOSED TERM MAY ALREADY BE PRESENT IN MAT!
//    if( options.BASIS == Options::MONOM && !ijmon.first->tord && ijmon.second->tord ){
//      std::set< std::pair< SPolyMon const*, SPolyMon const* >, lt_SQuad > CandidateDec;
//      _candidates( CandidateDec, *ijmon.second );
//      ijmon.first  = CandidateDec.crbegin()->first;
//      ijmon.second = CandidateDec.crbegin()->second;
//    }
    // Search connections between entry <a>ijmon</a> and existing terms
    for( auto it=submat.begin(); it!=submat.end(); ){
      // not connected to term
      if( !_isconnected( ijmon, *it ) ){
        ++it;
      }
      // merge terms if more than one connection 
      else if( itsubmat != submat.end() ){
        itsubmat->insert( it->begin(), it->end() );
        it = submat.erase( it );
      }
      // first connection to term
      else{
        (*it)[ijmon] = coef;
        itsubmat = it;
        ++it;
      }
    }
    if( itsubmat != submat.end() ) continue;
    // Create new independent term  if no connections
    submat.push_back( t_SQuad() );
    submat.back()[ijmon] = coef;
  }
  
#ifdef MC__SQUAD_DEBUG_SEPARATE
  unsigned i = 0;
  for( auto const& mat : submat ){
    std::cout << std::endl << "  Separable terms #" << ++i << ":" << std::endl;
    for( auto const& term : mat )
      std::cout << "   " << std::right << std::setw(options.DISPLEN+7) << term.second << "   "
                << " ( " << term.first.first->display( options.BASIS )
                << " ; " << term.first.second->display( options.BASIS ) << " )"
                << std::endl;
  }
#endif
  return submat;
}

inline bool
SQuad::_isconnected
( std::pair< SPolyMon const*, SPolyMon const* > const& entry, 
  SQuad::t_SQuad const& mat )
const
{
  for( auto const& [ijmon,coef] : mat ){
    if( entry.first->tord && ( entry.first == ijmon.first || entry.first == ijmon.second ) )
    //if( entry.first->tord && ( entry.first->inter( *ijmon.first ) || entry.first->inter( *ijmon.second ) ) )
      return true;
    if( entry.second->tord && ( entry.second == ijmon.first || entry.second == ijmon.second ) )
    //if( entry.second->tord && ( entry.second->inter( *ijmon.first ) || entry.second->inter( *ijmon.second ) ) )
      return true;
  }
  return false;  
}

inline std::multimap< double, std::map< SPolyMon, double, lt_SPolyMon > >
SQuad::factorize
( SQuad::t_SQuad const& mat )
const
{
  // Populate sparse symmetric coefficient matrix
  CPPL::dssmatrix coefmat;
  std::map< SPolyMon, unsigned, lt_SPolyMon > indexmap;
  unsigned index = 0;
  for( auto const& [ijmon,coef] : mat ){
    auto itimon = indexmap.find( *ijmon.first );
    if( itimon == indexmap.end() ){
      itimon = indexmap.insert( std::make_pair( *ijmon.first, index++ ) ).first;
      coefmat.stretch( 1 );
    }
    // Diagonal entry
    if( ijmon.first == ijmon.second ){
      coefmat.put( itimon->second, itimon->second, coef );
    }
    // Off-diagonal entry
    else{
      auto itjmon = indexmap.find( *ijmon.second );
      if( itjmon == indexmap.end() ){
        itjmon = indexmap.insert( std::make_pair( *ijmon.second, index++ ) ).first;
        coefmat.stretch( 1 );
      }
      coefmat.put( itimon->second, itjmon->second, coef/2e0 );
    }
  }
  
  // Perform eigenvalue decomposition
  std::multimap< double, std::map< SPolyMon, double, lt_SPolyMon > > eigdec;
  std::vector<double> eigval;
  std::vector<CPPL::dcovector> eigvec;
  CPPL::dsymatrix dcoefmat = coefmat.to_dsymatrix();
  if( dcoefmat.dsyev( eigval, eigvec ) ) return eigdec;

  // Populate decomposition map
  auto itval = eigval.begin();
  auto itvec = eigvec.begin();
  for( ; itval != eigval.end(); ++itval, ++itvec ){
    std::map< SPolyMon, double, lt_SPolyMon > eigterm;
    for( auto const& [mon,index] : indexmap ){
      //if( isequal( itvec->array[index], 0. ) ) continue;
      eigterm[mon] = itvec->array[index];
    }
    eigdec.insert( std::make_pair( *itval, eigterm ) );
  }

#ifdef MC__SQUAD_DEBUG_SEPARATE
  unsigned i = 0;
  for( auto const& [eigval,eigterm] : eigdec ){
    std::cout << std::endl << "  Eigen-direction #" << ++i << ": " << std::scientific 
              << std::right << std::setw(options.DISPLEN+7) << eigval << std::endl;
    for( auto const& [mon,coord] : eigterm )
      std::cout << "  " << mon.display(options.BASIS) << ": "
                << std::right << std::setw(options.DISPLEN+7) << coord << std::endl;
  }
#endif
  return eigdec;
}

inline double
SQuad::check
( unsigned const nSPol, SQuad::t_SPolyMonCoef const* pSPol, int const BASIS )
const
{ 
  assert( nSPol <= _MatFct.size() );
  SPolyExpr::options.BASIS = options.BASIS;
  double sumdiff = 0e0;
#ifdef MC__SQUAD_DEBUG_CHECK
  std::cout << *this << std::endl;
#endif
  
  // Process entries in _SetMon
  FFGraph dag;
  std::map< unsigned, FFVar > mapvar;
  std::map< SPolyMon, SPolyExpr, lt_SPolyMon > mapmon; 
  for( auto const& mon : _SetMon ){
    if( mon.tord == 0 ){
      mapmon[mon] = SPolyExpr( 1e0 );
    }
    else if( mon.tord == 1 ){
      unsigned const ivar = mon.expr.begin()->first;
      mapvar[ivar] = FFVar( &dag );
      mapmon[mon] = std::make_pair( FFMon( mapvar[ivar] ), 1e0 );
    }
    else{
      FFMon newmon;
      for( auto const& [ivar,iord] : mon.expr )
        assert( newmon.expr.insert( std::make_pair( &mapvar[ivar], iord ) ).second );
      newmon.tord = mon.tord;
      mapmon[mon] = std::make_pair( newmon, 1e0 );
    }
#ifdef MC__SQUAD_DEBUG_CHECK
    std::cout << mapmon[mon];
#endif
  }

  // Check entries in _MatFct
  auto itmat = _MatFct.begin();
  std::advance( itmat, _MatFct.size()-nSPol );
  for( unsigned i=0; nSPol && itmat != _MatFct.end(); ++itmat, ++i ){
    SPolyExpr spe( 0e0 );
    auto coefmon = pSPol[i];
    if( BASIS != options.BASIS ){
      t_SPolyVar ndxvar;
      _deps( coefmon, ndxvar );
      _convert( coefmon, ndxvar, BASIS );
    }
    for( auto const& [mon,coef] : coefmon ){
      if( mon.tord == 0 ){
        spe += coef;
        continue;
      }
      auto itmon = mapmon.find( mon );
      if( itmon != mapmon.end() ){
        spe += coef * itmon->second;
      }
      else{ 
        FFMon newmon;
        for( auto const& [ivar,iord] : mon.expr )
          assert( newmon.expr.insert( std::make_pair( &mapvar[ivar], iord ) ).second );
        newmon.tord = mon.tord;
        spe += SPolyExpr( std::make_pair( newmon, coef ) );
      }
    }

    for( auto const& [ijmon,coef] : *itmat )
      spe -= coef * mapmon[*ijmon.first] * mapmon[*ijmon.second];
#ifdef MC__SQUAD_DEBUG_CHECK
    std::cout << "\n  Quadratic form of P[" << _MatFct.size()-nSPol+i << "]:" << spe;
#endif
    double locdiff = 0;
    for( auto const& [mon,coef] : spe.mapmon() )
      locdiff += std::fabs( coef ); 
    if( std::fabs(locdiff) > 1e-10 )
      std::cerr << "\n  Error in quadratic form of P[" << _MatFct.size()-nSPol+i << "]:" << spe;
    sumdiff += locdiff;  
  }

  // Check entries in _MatRed
  itmat = _MatRed.begin();
  for( unsigned i=0; itmat != _MatRed.end(); ++itmat, ++i ){
    SPolyExpr spe( 0e0 );
    for( auto const& [ijmon,coef] : *itmat ){
      //std::cout << mapmon[*ijmon.first] << " * " << mapmon[*ijmon.second];
      spe += coef * mapmon[*ijmon.first] * mapmon[*ijmon.second];
    }
#ifdef MC__SQUAD_DEBUG_CHECK
    std::cout << " \n Auxiliary quadratic form #" << i << spe;
#endif
    double locdiff = 0;
    for( auto const& [mon,coef] : spe.mapmon() )
      locdiff += std::fabs( coef ); 
    if( std::fabs(locdiff) > 1e-10 )
      std::cerr << " \n Error in auxiliary quadratic form #" << i << spe;
    sumdiff += locdiff;  
  }

  // Check entries in _MatPSD
  itmat = _MatPSD.begin();
  for( unsigned i=0; itmat != _MatPSD.end(); ++itmat, ++i ){
    SPolyExpr spe( 0e0 );
    for( auto const& [ijmon,coef] : *itmat ){
      //std::cout << mapmon[*ijmon.first] << " * " << mapmon[*ijmon.second];
      spe += coef * mapmon[*ijmon.first] * mapmon[*ijmon.second];
    }
#ifdef MC__SQUAD_DEBUG_CHECK
    std::cout << " \n Positive semi-definite quadratic form #" << i << spe;
#endif
    double locdiff = 0;
    for( auto const& [mon,coef] : spe.mapmon() )
      locdiff += std::fabs( coef ); 
    if( std::fabs(locdiff) > 1e-10 )
      std::cerr << " \n Error in auxiliary quadratic form #" << i << spe;
    sumdiff += locdiff;  
  }

  return sumdiff;
}

inline void
SQuad::_pow2cheb
( std::vector< SQuad::t_SPolyMonCoef >& veccoef )
{
  // Converts COEFF in monomial basis to the Chebyshev basis.
  // Based on SCONMC function (http://www.netlib.org/math/MATH77/dconmc.f, http://www.netlib.org/math/docpdf/ch11-03.pdf)
  int const N = veccoef.size()-1;
  if( N <= 1 ) return;
  // TP = .5D0**(N-1)
  // COEFF(N) = TP * COEFF(N)
  // COEFF(N-1) = TP * COEFF(N-1)
  double TP = std::pow(0.5,N-1);
  for( auto& [mon,coef] : veccoef[N] )   coef *= TP;
  for( auto& [mon,coef] : veccoef[N-1] ) coef *= TP;
  //    do 20 J = N-2, 0, -1
  //      TP = 2.D0 * TP
  //      COEFF(J) = TP * COEFF(J)
  //      COEFF(J+1) = 2.D0 * COEFF(J+1)
  for( int J=N-2; J>=0; --J ){
    TP *= 2e0;
    for( auto& [mon,coef] : veccoef[J] )   coef *= TP;
    for( auto& [mon,coef] : veccoef[J+1] ) coef *= 2e0;
    //    do 10 I = J, N-2
    //      COEFF(I) = COEFF(I) + COEFF(I+2)
    // 10    continue
    for( int I=J; I<=N-2; ++I )
      for( auto& [mon,coef] : veccoef[I+2] ){
        auto [itmon,ins] = veccoef[I].insert( std::make_pair( mon, coef ) );
        if( !ins ) itmon->second += coef;
        if( itmon->second == 0e0 ) veccoef[I].erase( itmon );
      }
  // 20 continue
  }
  //    return
  //    end  
}

inline void
SQuad::_cheb2pow
( std::vector< SQuad::t_SPolyMonCoef >& veccoef )
{
  // Converts COEFF in Chebyshev basis to the monomial basis.
  // Based on SCONCM function (http://www.netlib.org/math/MATH77/dconcm.f, http://www.netlib.org/math/docpdf/ch11-03.pdf)
  int const N = veccoef.size()-1;
  if( N <= 1 ) return;
  //    TP = 1.D0
  double TP = 1e0;
  //    do 20 J = 0, N-2
  for( int J=0; J<=N-2; ++J ){
  //       do 10 I = N-2, J, -1
  //          COEFF(I) = COEFF(I) - COEFF(I+2)
  // 10    continue
    for( int I=N-2; I>=J; --I )
      for( auto& [mon,coef] : veccoef[I+2] ){
        auto [itmon,ins] = veccoef[I].insert( std::make_pair( mon, -coef ) );
        if( !ins ) itmon->second -= coef;
        if( itmon->second == 0e0 ) veccoef[I].erase( itmon );
      }
  //       COEFF(J+1) = .5D0 * COEFF(J+1)
  //       COEFF(J) = TP * COEFF(J)
  //       TP = 2.D0 * TP
  // 20 continue
    for( auto& [mon,coef] : veccoef[J+1] ) coef /= 2e0;
    for( auto& [mon,coef] : veccoef[J] )   coef *= TP;
    TP *= 2e0;
  }
  //    COEFF(N) = TP * COEFF(N)
  //    COEFF(N-1) = TP * COEFF(N-1)
  for( auto& [mon,coef] : veccoef[N] )   coef *= TP;
  for( auto& [mon,coef] : veccoef[N-1] ) coef *= TP;
  //    return
  //    end  
}

inline void
SQuad::_unipol
( unsigned const ivar, SQuad::t_SPolyMonCoef const& SPol,
  std::vector< SQuad::t_SPolyMonCoef >& veccoef )
{
  // Get univariate polynomial in the variable ivar
  veccoef.clear();
  unsigned const N = _maxord( ivar, SPol );
  veccoef.assign( N+1, t_SPolyMonCoef() );
  
  // No dependence
  if( !N ){
    veccoef[0] = SPol;
    return;
  }
  
  // Construct first- and higher-order terms
  for( auto const& [mon,coef] : SPol ){
    auto ie = mon.expr.find( ivar );
    if( ie == mon.expr.end() ) // no dependence on variable ivar 
      veccoef[ 0 ].insert( std::make_pair( mon, coef ) );
    else{
      auto const& [ivar,iord] = *ie;
      SPolyMon monmod( mon.tord - iord, mon.expr );
      monmod.expr.erase( ivar ); // remove ivar entry
      veccoef[ iord ].insert( std::make_pair( monmod, coef ) );
    }
  }
}

inline unsigned
SQuad::_maxord
( unsigned const ivar, SQuad::t_SPolyMonCoef const& SPol )
{
  unsigned maxord = 0;
  for( auto const& [mon,coef] : SPol ){
    auto ie = mon.expr.find( ivar );
    if( ie != mon.expr.end() && ie->second > maxord )
      maxord = ie->second;
  }
  return maxord;
}

inline bool
SQuad::_deps
( unsigned const nSPol, SQuad::t_SPolyMonCoef const* pSPol, SQuad::t_SPolyVar& SVar )
{
  for( unsigned i=0; i<nSPol; i++ )
    _deps( pSPol[i], SVar );
  return !SVar.empty();
}

inline bool
SQuad::_deps
( SQuad::t_SPolyMonCoef const& SPol, SQuad::t_SPolyVar& SVar )
{
  for( auto const& [mon,coef] : SPol )
    for( auto const& [ivar,iord] : mon.expr )
      SVar.insert( ivar );
  return !SVar.empty();
}

inline void
SQuad::_convert
( SQuad::t_SPolyMonCoef& SPol, SQuad::t_SPolyVar const& SVar,
  int const BASIS )
{
  // Constant or linear polynomial - nothing to do
  if( SPol.empty() || SVar.empty() || SPol.crbegin()->first.tord <= 1 ) return;

  // Apply conversion variable by variable
  std::vector< SQuad::t_SPolyMonCoef > veccoef;
  for( unsigned ivar : SVar ){
    // Separate contribution of variable ivar, and convert to Chebyshev form
    _unipol( ivar, SPol, veccoef );
    switch( BASIS ){
      case Options::MONOM: _pow2cheb( veccoef ); break;
      case Options::CHEB:  _cheb2pow( veccoef ); break;
    }
    // Merge back into multivariate polynomial
    SPol = veccoef[0];
    for( unsigned iord=1; iord<veccoef.size(); iord++ ){
      for( auto& [mon,coef] : veccoef[iord] ){
        SPolyMon monmod( mon.tord + iord, mon.expr );
        monmod.expr[ ivar ] = iord;
        SPol[ monmod ] = coef;
      }
    }
    //_simplify( coefmon );
  }
}

inline double
SQuad::process
( unsigned const nSPol, SQuad::t_SPolyMonCoef const* pSPol,
  int const BASIS, bool const CHECK )
{
  for( unsigned i=0; i<nSPol; i++ )
    process( pSPol[i], BASIS, false );
  return( CHECK? check( nSPol, pSPol, BASIS ): 0. );
}

inline double
SQuad::process
( std::set<unsigned> const& ndxSPol, SQuad::t_SPolyMonCoef const* pSPol,
  int const BASIS, bool const CHECK )
{
  std::vector<t_SPolyMonCoef> vpSPol;
  vpSPol.reserve( ndxSPol.size() );
  for( unsigned const& i : ndxSPol ) vpSPol.push_back( pSPol[i] );
  return process( ndxSPol.size(), vpSPol.data(), BASIS, CHECK );
}

inline double
SQuad::process
( SQuad::t_SPolyMonCoef const& SPol, int const BASIS,
  bool const CHECK )
{
  // Append entry in <a>MatFct</a>
  unsigned ndxmat = _MatFct.size();
  _MatFct.push_back( t_SQuad() );
  auto& mat = _MatFct[ndxmat];

  // Always insert constant monomial 1
  _SetMon.insert( SPolyMon() );

  // Identify and insert variables participating in monomials
  t_SPolyVar ndxvar;
  _deps( SPol, ndxvar );
  for( auto ivar : ndxvar ) _SetMon.insert( SPolyMon( ivar ) );

  // Local copy and conversion to desired basis
  auto coefmon = SPol;
  if( BASIS != options.BASIS ) _convert( coefmon, ndxvar, BASIS );

  // Iterate through monomial terms
  //for( auto it=coefmon.crbegin(); it!=coefmon.crend(); ++it ){
    //const auto& [mon,coef] = *it;
  for( ; !coefmon.empty(); ){
    // Local copy of next monomial, then erase
    auto const [mon,coef] = options.ORDER==Options::INC? *coefmon.cbegin(): *coefmon.crbegin();
    coefmon.erase( mon );

    // Insert variables participating in monomial
    //for( auto const& [ivar,iord] : mon.expr )
    //  _SetMon.insert( SPolyMon( ivar ) );

    // Monomial already present in _SetMon
    auto itmon = _SetMon.find( mon );
    if( itmon != _SetMon.end() ){
#ifdef MC__SQUAD_DEBUG_DECOMP
      std::cout << "Inserted: " << itmon->display(options.BASIS)
                << " = " << _SetMon.cbegin()->display(options.BASIS)
                << " · " << itmon->display(options.BASIS)
                << std::endl;
#endif
      bool ins = _insert( mat, &(*itmon), coef, true );  
#ifdef MC__SQUAD_CHECK
      assert( ins );
#endif
      continue;
    }
      
    // Mononial needs further decomposition
    auto const& [plmon,prmon] = _decompose( mon );
#ifdef MC__SQUAD_CHECK
    assert( plmon && prmon );
#endif
#ifdef MC__SQUAD_DEBUG_DECOMP
  std::cout << "Inserted: " << mon.display(options.BASIS)
            << " = " << plmon->display(options.BASIS)
            << " · " << prmon->display(options.BASIS)
            << std::endl;
#endif
    bool ins = _insert( mat, coefmon, plmon, prmon, coef );
    //bool ins = _insert( _MatFct[ndxmat], coefmon, plmon, prmon, coef );
#ifdef MC__SQUAD_CHECK
    assert( ins );
#endif
  }

  return( CHECK? check( 1, &SPol, BASIS ): 0. );
}

inline std::pair< SPolyMon const*, SPolyMon const* >
SQuad::_reorder
( std::pair< SPolyMon const*, SPolyMon const* > pMon )
const
{
  if( lt_SPolyMon()( *pMon.second, *pMon.first ) )
    std::swap( pMon.first, pMon.second );
  return pMon;
}

inline std::set< SPolyMon, lt_SPolyMon >
SQuad::_prodmon
( SPolyMon const& mon1, SPolyMon const& mon2 )
const
{
#ifdef MC__SQUAD_DEBUG_PRODMON
  std::cout << "mon1:\n" << mon1.display(options.BASIS) << std::endl;
  std::cout << "mon2:\n" << mon2.display(options.BASIS) << std::endl;
#endif
  std::set<SPolyMon,lt_SPolyMon> prodmon;
  switch( options.BASIS ){
   // Monomial basis representation
   case Options::MONOM:
     prodmon.insert( mon1+mon2 );
     break;

   // Chebyshev basis representation
   case Options::CHEB:
    prodmon.insert( SPolyMon() );
    for( auto const& [ivar,iord] : (mon1+mon2).expr ){
      std::set<SPolyMon,lt_SPolyMon> prodmon2;
      auto&& it1 = mon1.expr.find( ivar );
      auto&& it2 = mon2.expr.find( ivar );
      if( it1 != mon1.expr.end() && it2 != mon2.expr.end() ){
        unsigned const& iord1 = it1->second;
        unsigned const& iord2 = it2->second;
        assert( iord1 && iord2 );
        for( auto const& mon3 : prodmon ){
          prodmon2.insert( mon3 + SPolyMon( iord1+iord2, {std::make_pair( ivar, iord1+iord2 )} ) );
          if( iord1 > iord2 )
            prodmon2.insert( mon3 + SPolyMon( iord1-iord2, {std::make_pair( ivar, iord1-iord2 )} ) );
          else if( iord1 < iord2 )
            prodmon2.insert( mon3 + SPolyMon( iord2-iord1, {std::make_pair( ivar, iord2-iord1 )} ) );
          else
            prodmon2.insert( mon3 );
        }
      }
      else if( it1 != mon1.expr.end() ){
        unsigned const& iord1 = it1->second;
        assert( iord1 );
        for( auto const& mon3 : prodmon )
          prodmon2.insert( mon3 + SPolyMon( iord1, {std::make_pair( ivar, iord1 )} ) );
      }    
      else{
        unsigned const& iord2 = it2->second;
        assert( iord2 );
        for( auto const& mon3 : prodmon )
          prodmon2.insert( mon3 + SPolyMon( iord2, {std::make_pair( ivar, iord2 )} ) );
      }    
      std::swap( prodmon, prodmon2 );
    }
    break;
  }
#ifdef MC__SQUAD_DEBUG_PRODMON
  std::cout << "monprod:\n";
  for( auto const& mon : prodmon )
    std::cout << mon.display(options.BASIS) << std::endl;
#endif
  return prodmon;
}

inline bool
SQuad::_insert
( t_SQuad& mat, SPolyMon const* pMon, double const coef, bool const add )//, bool const rem )
{
  // New entry in quadratic form as product with constant monomial
#ifdef MC__SQUAD_DEBUG_DECOMP
  std::cerr << "SQuad::_insert, &mat = " << &mat << std::endl;
#endif
  auto [itmat,ins] = mat.insert( std::make_pair( std::make_pair( &(*_SetMon.cbegin()), pMon ), coef ) );
  if( !ins && add ){
    itmat->second += coef;
    if( itmat->second == 0. ) mat.erase( itmat );
  }
  return ins || add;
  //if( !ins && rem ) mat.erase( itmat );
  //return ins || rem;
}

inline bool
SQuad::_insert
( t_SQuad& mat, SPolyMon const* pLMon, SPolyMon const* pRMon,
  double const coef, bool const add )
{
  // New entry in quadratic form
  auto [itmat,ins] = mat.insert( std::make_pair( _reorder( std::make_pair( pLMon, pRMon ) ), coef ) );
  if( !ins && add ){
    itmat->second += coef;
    if( itmat->second == 0. ) mat.erase( itmat );
  }
  return ins || add;
}

inline bool
SQuad::_insert
( t_SQuad& mat, std::map< SPolyMon, double, lt_SPolyMon >& coefmon,
  SPolyMon const* pLMon, SPolyMon const* pRMon, double const coef,
  bool const add )
{
  // Extra lower-order terms generate by Chebyshev product
  auto&& prodmon = _prodmon( *pLMon, *pRMon );
  unsigned const nprod = prodmon.size();
  auto&& itmon = prodmon.crbegin(); 
  for( ++itmon; itmon != prodmon.crend(); ++itmon ){
    auto [itcmon,ins] = coefmon.insert( std::make_pair( *itmon, -coef ) );
    if( !ins ) itcmon->second -= coef; 
    if( itcmon->second == 0. ) coefmon.erase( itcmon );
  }
  
  // New entry in quadratic form
  auto [itmat,ins] = mat.insert( std::make_pair( _reorder( std::make_pair( pLMon, pRMon ) ), coef * nprod ) );
  if( !ins && add ){
    itmat->second += coef;
    if( itmat->second == 0. ) mat.erase( itmat );
  }
  return ins || add;
}

inline void
SQuad::_candidates
( std::set< std::pair< SPolyMon const*, SPolyMon const* >, lt_SQuad >& CandidateDec,
  SPolyMon const& mon )
const
{
  for( auto&& mon2 : _SetMon ){
    if( !mon2.subseteq( mon ) ) continue;
    auto&& itmon3 = _SetMon.find( mon - mon2 );
    if( itmon3 == _SetMon.end() || ( &mon2 != &*itmon3 && lt_SPolyMon()( mon2, *itmon3 ) ) ) continue;
    //if( itmon3 == _SetMon.end() || lt_SPolyMon()( *itmon3, mon2 ) ) continue;
#ifdef MC__SQUAD_DEBUG_DECOMP
    std::cout << "Candidate: " << mon.display(options.BASIS)
              << " = " << mon2.display(options.BASIS)
              << " · " << itmon3->display(options.BASIS)
              << std::endl;
#endif
    CandidateDec.insert( std::make_pair( &mon2, &(*itmon3) ) );
  }
}

inline std::pair< SPolyMon const*, SPolyMon const* >
SQuad::_decompose
( SPolyMon const& mon )
{
  // Possible decompositions using existing monomial in _SetMon
  std::set< std::pair< SPolyMon const*, SPolyMon const* >, lt_SQuad > CandidateDec;
  _candidates( CandidateDec, mon );
  
  // Case 1: Monomial can be decomposed in terms of existing monomial in _SetMon
  if( !CandidateDec.empty() ){
#ifdef MC__SQUAD_DEBUG_DECOMP
    std::cout << "Decomposed: " << mon.display(options.BASIS)
              << " = " << CandidateDec.rbegin()->first->display(options.BASIS)
              << " · " << CandidateDec.rbegin()->second->display(options.BASIS)
              << std::endl;
#endif
    // Prefered candidate is based on order defined by lt_SQuad
    // **unless** monomial already present 'as is'
    if( !CandidateDec.begin()->first->tord ) return *CandidateDec.begin();
    return *CandidateDec.rbegin();
  }

  // Case 2: Monomial is linear in all of the variables (multilinear)
  if( mon.gexp() == 1 ){
    t_pSPolyMon CandidateMon;
    for( auto&& mon2 : _SetMon ){
      if( !mon2.subseteq( mon ) ) continue;
#ifdef MC__SQUAD_CHECK
      assert( _SetMon.find( mon - mon2 ) == _SetMon.end() ); // Covered by Case 1
#endif
      CandidateMon.insert( &mon2 );
#ifdef MC__SQUAD_DEBUG_DECOMP
      std::cout << "Candidate: " << (mon-mon2).display(options.BASIS)
                << " = " << mon.display(options.BASIS)
                << " / " << mon2.display(options.BASIS)
                << std::endl;
#endif
    }
#ifdef MC__SQUAD_CHECK
    assert( !CandidateMon.empty() ); // _SetMon comprises the participating variables
#endif

    // Case 2a: Use existing non-trivial monomial component in _SetMon
    if( (*CandidateMon.rbegin())->tord > 1 ){
      SPolyMon const* pmon2 = *CandidateMon.rbegin();
      SPolyMon mon3( mon - *pmon2 );
#ifdef MC__SQUAD_DEBUG_DECOMP
      std::cout << "Decomposed: " << mon.display(options.BASIS)
                << " = " << pmon2->display(options.BASIS)
                << " · " << mon3.display(options.BASIS)
                << std::endl;
#endif
      // Decompose mon3
      auto&& itmon3 = _subexpression( mon3 );
      //return _reorder( std::make_pair( pmon2, &(*itmon3) ) );
      return std::make_pair( pmon2, &(*itmon3) );
    }

    // Case 2b: Split monomial into two monomials of similar total order
    unsigned count = 0;
    SPolyMon mon2;
    for( auto const& [ivar,iord] : mon.expr ){
#ifdef MC__SQUAD_CHECK
      assert( iord == 1 );
#endif
      mon2 += SPolyMon( ivar );
      if( ++count >= mon.tord / 2 + mon.tord % 2 ) break;
    }
    SPolyMon mon3( mon - mon2 );
#ifdef MC__SQUAD_DEBUG_DECOMP
    std::cout << "Decomposed: " << mon.display(options.BASIS)
              << " = " << mon2.display(options.BASIS)
              << " · " << mon3.display(options.BASIS)
              << std::endl;
#endif
    // Decompose mon2 and mon3
    auto&& itmon2 = _subexpression( mon2 );
    auto&& itmon3 = _subexpression( mon3 );
    //return _reorder( std::make_pair( &(*itmon2), &(*itmon3) ) );
    return std::make_pair( &(*itmon2), &(*itmon3) );
  }

  // Case 3: Monomial has partial order >1 in all of the variables w/ some odd partial order
  if( mon.gcexp() % 2 ){
  //if( mon.lexp() == 1 && mon.gexp() > 1  ){
    SPolyMon mon2;
    for( auto&& [ivar,iord] : mon.expr )
      if( iord % 2 )
        mon2 += SPolyMon( ivar );
#ifdef MC__SQUAD_DEBUG_DECOMP
    std::cout << "Decomposed: " << mon.display(options.BASIS)
              << " = " << mon2.display(options.BASIS)
              << " · " << (mon-mon2).display(options.BASIS)
              << std::endl;
#endif
    SPolyMon mon3( mon - mon2 );
    // Decompose mon2 and mon3
#ifdef MC__SQUAD_DEBUG_DECOMP
      std::cout << "Decomposing: " << mon2.display(options.BASIS) << std::endl;
#endif
    auto&& itmon2 = _SetMon.find( mon2 );
    if( itmon2 == _SetMon.end() ) itmon2 = _subexpression( mon2 );
#ifdef MC__SQUAD_DEBUG_DECOMP
      std::cout << "Decomposing: " << mon3.display(options.BASIS) << std::endl;
#endif
    auto&& itmon3 = _SetMon.find( mon3 );
    if( itmon3 == _SetMon.end() ) itmon3 = _subexpression( mon3 );
    //return _reorder( std::make_pair( &(*itmon2), &(*itmon3) ) );
    return std::make_pair( &(*itmon2), &(*itmon3) );
  }
  
  // Case 4: Monomial has even partial order in all of the variables
  SPolyMon mon2 = mon / 2;
  // Decompose mon2
  auto&& itmon2 = _SetMon.find( mon2 );
  if( itmon2 == _SetMon.end() ) itmon2 = _subexpression( mon2 );
  //return _reorder( std::make_pair( &(*itmon2), &(*itmon2) ) );   
  return std::make_pair( &(*itmon2), &(*itmon2) );
}

inline typename SQuad::t_SPolyMon::iterator
SQuad::_subexpression
( SPolyMon const& mon )
{
  // Monomial mon already in _SetMon - may happen for constant monomial and variables
  auto itmon0 = _SetMon.find( mon );
  if( itmon0 != _SetMon.end() ){
    //std::cerr << "SQuad::_subexpression: Existing monomial " << mon.display(options.BASIS) << std::endl;
    return itmon0;
  }

  // Perform further decomposition of monomial <a>mon</a>
  auto const& [plmon,prmon] = _decompose( mon );
#ifdef MC__SQUAD_CHECK
  assert( plmon && prmon );
#endif

  // Append new reduction constraint for <a>mon</a>
  unsigned ndxmat = _MatRed.size();
  _MatRed.push_back( t_SQuad() );
  auto& mat = _MatRed.back();
#ifdef MC__SQUAD_DEBUG_DECOMP
  std::cerr << "SQuad::_subexpression, &mat = " << &mat << std::endl;
#endif
  auto itmon = _SetMon.insert( mon ).first;
  std::map< SPolyMon, double, lt_SPolyMon > coefmon;
  bool ins = _insert( mat, &(*itmon), -1. )
          && _insert( mat, coefmon, plmon, prmon, 1. );
#ifdef MC__SQUAD_CHECK
    assert( ins );
#endif
  //for( auto it=coefmon.crbegin(); it!=coefmon.crend(); ++it ){
    //auto const& [monlow,coeflow] = *it;
  for( ; !coefmon.empty(); ){
    // Local copy of next monomial, then erase
    auto const [monlow,coeflow] = options.ORDER==Options::INC? *coefmon.cbegin(): *coefmon.crbegin();
    coefmon.erase( monlow );

    auto itmonlow = _subexpression( monlow );
    ins = _insert( _MatRed[ndxmat], &(*itmonlow), coeflow );
#ifdef MC__SQUAD_CHECK
    assert( ins );
#endif
  }
  
  // Search for extra reduction constraints for <a>mon</a>
  if( options.REDUC ) _reduction( *itmon );
   
  return itmon;
}

inline void
SQuad::_reduction
( SPolyMon const& mon )
{
  // Search for extra reduction constraints for <a>mon</a>
  // that don't add additional low-order monomials
  auto itmon  = _SetMon.find( mon );
  auto itmon2 = _SetMon.begin();
  for( ++itmon2; itmon2 != _SetMon.end(); ++itmon2 ){
  
    // Low-order terms generated by Chebyshev product
    auto&& prodmon12 = _prodmon( mon, *itmon2 );
    bool monmis = false;
    auto itprodmon12 = prodmon12.crbegin(); 
    // Check that all low-order terms are present -> OR INTRODUCE THEM?
    for( ++itprodmon12; !monmis && itprodmon12 != prodmon12.crend(); ++itprodmon12 )
       if( _SetMon.find( *itprodmon12 ) == _SetMon.end() ) 
         //_subexpression( *itprodmon12 );
         monmis = true;
    if( monmis ) continue;

    // Find alternative decompositions
    SPolyMon montot( mon + *itmon2 );
    auto itmon3 = itmon2;
    for( ++itmon3; itmon3 != _SetMon.end(); ++itmon3 ){
      if( !itmon3->subseteq( montot ) ) continue;
      auto itmon4 = _SetMon.find( montot - *itmon3 );

      // Prevent duplication
      switch( options.ORDER ){
       case Options::INC:      
        if( itmon4 == _SetMon.end() || lt_SPolyMon()( *itmon4, *itmon3 ) || itmon == itmon4 ) continue;
        break;
       case Options::DEC:
        if( itmon4 == _SetMon.end() || lt_SPolyMon()( *itmon3, *itmon4 ) || itmon == itmon3 ) continue;
        break;
      }

      // Low-order terms generated by Chebyshev product
      auto&& prodmon34 = _prodmon( *itmon3, *itmon4 );
      bool monmis = false;
      auto itprodmon34 = prodmon34.crbegin(); 

      // Check that all low-order terms are present -> OR INTRODUCE THEM?
      for( ++itprodmon34; !monmis && itprodmon34 != prodmon34.crend(); ++itprodmon34 )
         if( _SetMon.find( *itprodmon34 ) == _SetMon.end() ) 
           //_subexpression( *itprodmon34 );
           monmis = true;
      if( monmis ) continue;
#ifdef MC__SQUAD_DEBUG_REDUC
      std::cout << "Reduction: ";
      if( prodmon12.size() > 1 ) std::cout << prodmon12.size() << " · ";
      std::cout << mon.display(options.BASIS) << " · " << itmon2->display(options.BASIS);
      itprodmon12 = prodmon12.crbegin(); 
      for( ++itprodmon12; itprodmon12 != prodmon12.crend(); ++itprodmon12 )
        std::cout << " - " << _SetMon.find( *itprodmon12 )->display(options.BASIS);
      std::cout << " == ";
      if( prodmon34.size() > 1 ) std::cout << prodmon34.size() << " · ";
      std::cout << itmon3->display(options.BASIS) << " · " << itmon4->display(options.BASIS);
      itprodmon34 = prodmon34.crbegin(); 
      for( ++itprodmon34; itprodmon34 != prodmon34.crend(); ++itprodmon34 )
        std::cout << " - " << _SetMon.find( *itprodmon34 )->display(options.BASIS); 
      std::cout << std::endl;
#endif

      // Populate reduction constraint
      _MatRed.push_back( t_SQuad() );
      auto& mat = _MatRed.back();
      // Insert product terms
      bool ins = _insert( mat, &mon,       &(*itmon2), -(double)prodmon12.size() )
              && _insert( mat, &(*itmon3), &(*itmon4),  (double)prodmon34.size() );
#ifdef MC__SQUAD_CHECK
      assert( ins );
#endif
      // Insert low-order terms 
      itprodmon12 = prodmon12.crbegin(); 
      for( ++itprodmon12; itprodmon12 != prodmon12.crend(); ++itprodmon12 ){
        ins = _insert( mat, &(*_SetMon.find( *itprodmon12 )),  1. );
#ifdef MC__SQUAD_CHECK
        assert( ins );
#endif
      }
      itprodmon34 = prodmon34.crbegin(); 
      for( ++itprodmon34; itprodmon34 != prodmon34.crend(); ++itprodmon34 ){
        ins = _insert( mat, &(*_SetMon.find( *itprodmon34 )), -1., true ); // Element removed if already present
#ifdef MC__SQUAD_CHECK
        assert( ins );
#endif
      }
    }
  }
}

inline void
SQuad::_reset
()
{
  _SetMon.clear();
  _MatFct.clear();
  _MatRed.clear();
  _MatPSD.clear();
}

} // namespace mc

#endif
