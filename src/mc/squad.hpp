// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_SQUAD Decomposition of Sparse Polynomials into Quadratic Forms
\author Benoit Chachuat, Tanuj Karia & OMEGA Research Group (http://www3.imperial.ac.uk/environmentenergyoptimisation)
\date 2020
\bug mc::SQuad::Optimize only working with default argument OnlyminOrd=2 presently.

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

It is convenient to leverage MC++'s sparse polynomial class in <A>spoly.hpp</a>
\code
      #include "spoly.hpp"
\endcode
to help construct the polynomials in the desired format:
\code
      unsigned const NX = 3, NP = 2;
      mc::SPoly<> X[NX], P[NP];
      for( unsigned i=0; i<NX; i++ ) X[i].var( i );
      P[0] = pow( X[0] + sqr( X[1] ) - 2 * X[2], 3 );
      P[1] = 2 * sqr( X[1] ) - 1;

      std::cout << "\nSparse multivariate polynomials:";
      for( unsigned i=0; i<NP; i++ )
        std::cout << P[i];
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
An environment <a>mc::SQuad</a> is defined and the sparse multivariate polynomials are quadratized by calling the method <a>mc::SQuad::process</a> for each sparse polynomial. This method takes a coefficient map (or an array of maps) as argument, which is of type <a>std::map<mc::SMon,double,mc::lt_SMon></a> where the class mc::SMon is used to store and manipulate sparse monomials and mc::lt_SMon implements the grlex order:
\code
      mc::SQuad<> SQF;
      SQF.process( NP, P, &mc::SPoly<>::mapmon );
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

These results show that a total of 12 monomials participate in the quadratic forms: The constant monomial \f$\xi_0=1\f$; the participating variables \f$\xi_1:=x_0\f$, \f$\xi_2:=x_1\f$, \f$\xi_3=x_2\f$; and 8 lifted monomials \f$x_4:=x_0^2\f$, \f$\xi_5:=x_0\cdot x_1\f$, \f$\xi_6:=x_0\cdot x_2\f$, \f$\xi_7:=x_1^2\f$, \f$\xi_8:=x_1\cdot x_2\f$, \f$\xi_9:=x_2^2\f$, \f$\xi_{10}:=x_1^3\f$, \f$\xi_{11}:=x_1^4\f$. A reformulation of the polynomials \f$p_0,p_1\f$ in terms of these monomials is given by:
\f{align*}
  q_0(\xi_0,\ldots,\xi_{11}) =\ & \xi_1\cdot \xi_4 - 6 \xi_1\cdot \xi_6 + 12 \xi_1\cdot \xi_9 + 3 \xi_1\cdot \xi_{11} - 8 \xi_3\cdot \xi_9 - 6 \xi_3\cdot \xi_{11} + 3 \xi_5^2 - 12 \xi_6\cdot \xi_7 + 12 \xi_8^2 + \xi_{10}^2\\
  q_1(\xi_0,\ldots,\xi_{11}) =\ & - \xi_0^2 + 2 \xi_0\cdot\xi_7
\f}
with the following 8 auxiliary quadratic constraints are also generated:
\f{align*}
  \xi_1^2 - \xi_0\cdot\xi_4 =\ & 0\\
  \xi_1\cdot\xi_2 - \xi_0\cdot\xi_5 =\ & 0\\
  \xi_1\cdot\xi_3 - \xi_0\cdot\xi_6 =\ & 0\\
  \xi_2^2 - \xi_0\cdot\xi_7 =\ & 0\\
  \xi_2\cdot\xi_3 - \xi_0\cdot\xi_8 =\ & 0\\
  \xi_3^2 - \xi_0\cdot\xi_9 =\ & 0\\
  \xi_2\cdot\xi_7 - \xi_0\cdot\xi_{10} =\ & 0\\
  \xi_2\cdot\xi_{10} - \xi_0\cdot\xi_{11} =\ & 0\\
\f}

With the redundant constraint option activated before the quadratization:
\code
      SQuad<> SQF;
      SQuad<>::options.REDUC = true;
      SQF.process( NP, P, &mc::SPoly<>::mapmon );
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
  \xi_7^2 - \xi_2\cdot\xi_{10} =\ & 0\\
  \xi_7\cdot\xi_{10} - \xi_2\cdot\xi_{11} =\ & 0\\
  \xi_{10}^2 - \xi_7\cdot\xi_{11} =\ & 0\\
  \xi_2\cdot\xi_8 - \xi_3\cdot\xi_7 =\ & 0\\
  \xi_3\cdot\xi_{10} - \xi_7\cdot\xi_8 =\ & 0\\
  \xi_3\cdot\xi_{11} - \xi_8\cdot\xi_{10} =\ & 0\\
  \xi_1\cdot\xi_8 - \xi_2\cdot\xi_6 =\ & 0\\
  \xi_1\cdot\xi_7 - \xi_2\cdot\xi_5 =\ & 0\\
  \xi_2\cdot\xi_6 - \xi_3\cdot\xi_5 =\ & 0\\
  \xi_1\cdot\xi_8 - \xi_3\cdot\xi_5 =\ & 0\\
  \xi_5\cdot\xi_7 - \xi_1\cdot\xi_{10} =\ & 0\\
  \xi_5\cdot\xi_{10} - \xi_1\cdot\xi_{11} =\ & 0\\
  \xi_3\cdot\xi_6 - \xi_1\cdot\xi_9 =\ & 0\\
  \xi_3\cdot\xi_8 - \xi_2\cdot\xi_9 =\ & 0\\
  \xi_6\cdot\xi_8 - \xi_5\cdot\xi_9 =\ & 0\\
  \xi_8^2 - \xi_7\cdot\xi_9 =\ & 0\\
  \xi_1\cdot\xi_5 - \xi_2\cdot\xi_4 =\ & 0\\
  \xi_1\cdot\xi_6 - \xi_3\cdot\xi_4 =\ & 0\\
\f}

The quadratization may also be performed in Chebyshev basis instead of the default monomial basis. And the monomials may be processed in increasing grlex order rather than the default decreasing grlex order.

Pointers of type mc::SMon to the monomials participating in the quadratic forms can be retrieved with the method mc::SQuad::SetMon, which is of type std::set<SMon const*,lt_pSMon>.

The sparse coefficient matrices defining the main quadratic forms for each processed multivariate polynomial and the corresponding auxiliary quadratic constraints can be retrieved with the methods mc::SQuad::MatFct and mc::SQuad::MatRed, respectively. These coefficient matrices are of the type std::map<std::pair<SMon const*,SMon const*>,double,lt_SQuad> where the comparison operator mc::lt_SQuad orders monomial pairs in grlex order for both elements. 


\section sec_SQUAD_optim How do I minimize the number of auxiliary variables in my reformulation?

Having decomposed a set of multivariate polynomials into quadratic forms using Algorithm 1, one may seek to further improve the decomposition by minimizing the number of auxiliary variables needed. Altough computing such a minimal decomposition is NP-hard in general, it can be automated using mixed-integer programming (MIP) and may be tractable for simple multivariate polynomials. The MIP model implemented in class mc::SQuad is the following:

\f{align*}
\displaystyle\min_{\boldsymbol{z},\boldsymbol{\nu}^{\rm L},\boldsymbol{\nu}^{\rm R},\boldsymbol{\beta},\boldsymbol{\omega}^{\rm L},\boldsymbol{\omega}^{\rm R}}\ & \sum_{k=1}^{n_a} z_k\\
\displaystyle\text{s.t.}\ \ \ & z_k \geq z_{k+1},\ \ k=1\ldots n_a-1\\
& \sum_{i=1}^{n_x+k-1} \nu^{\rm L}_{k,i} = z_k,\ \ k=1\ldots n_a\\
& \sum_{i=1}^{n_x+k-1} \nu^{\rm R}_{k,i} = z_k,\ \ k=1\ldots n_a\\
& \sum_{i=1}^{n_x+n_a} \omega^{\rm L}_{j,i} = 1,\ \ j=1\ldots n_m\\
& \sum_{i=1}^{n_x+n_a} \omega^{\rm R}_{j,i} \leq 1,\ \ j=1\ldots n_m\\
& \beta_{k,i} = \nu^{\rm L}_{k,i} + \nu^{\rm R}_{k,i} + \sum_{l=1}^{k-1} \left(\nu^{\rm L}_{k,l} + \nu^{\rm R}_{k,l}\right) \beta_{l,i},\ \ i=1\ldots n_x,\ \ k=1\ldots n_a\\
& \alpha_{j,i} = \omega^{\rm L}_{j,i} + \omega^{\rm R}_{j,i} + \sum_{l=1}^{n_a} \left(\omega^{\rm L}_{j,l} + \omega^{\rm R}_{j,l}\right) \beta_{l,i},\ \ i=1\ldots n_x,\ \ j=1\ldots n_m\\
& \beta_{k,i} \leq z_k\ \max\{\alpha_{j,i}: j=1\ldots n_m\},\ \ i=1\ldots n_x,\ \ k=1\ldots n_a\\
& \sum_{i=1}^{n_x} \beta_{k,i} \leq \max\left\{\sum_{i=1}^{n_x}\alpha_{j,i}: j=1\ldots n_m\right\}-1,\ \ k=1\ldots n_a\\
& z_k, \nu^{\rm L}_{k,i}, \nu^{\rm R}_{k,i}, \omega^{\rm L}_{j,i}, \omega^{\rm R}_{j,i}\in\{0,1\},\ \ \beta_{k,i}\geq 0,\ \ i=1\ldots n_x+n_a,\ \ j=1\ldots n_m,\ \ k=1\ldots n_a
\f}
where the following sets and variables are used:
- \f$n_m\f$, number of monomial terms to decompose, given by \f$m_{j} := \xi_1^{\alpha_{j,1}}\cdots \xi_{n_x}^{\alpha_{j,n_x}},\ j = 1\ldots n_m\f$
- \f$n_x\f$, original number of variables, \f$x_i=\xi_i,\ i=1\ldots n_x\f$ 
- \f$n_a\f$, maximal number of auxiliary variables; e.g., determined using the heuristic approach in Algorithm 1
- \f$z_k\f$, whether the auxiliary \f$\xi_{n_x+k},\ k=1\ldots n_a\f$ is used in decomposition
- \f$\nu^{\rm L}_{k,i},\nu^{\rm R}_{k,i}\f$, whether the variable or auxiliary \f$\xi_{i},\ i=1\ldots n_x+k\f$ decomposes the auxiliary \f$\xi_{n_x+k},\ k=1\ldots n_a\f$
- \f$\beta_{k,i}\f$, integer power of variable \f$\xi_{i},\ i=1\ldots n_x\f$ in the auxiliary \f$\xi_{n_x+k},\ k=1\ldots n_a\f$
- \f$\omega^{\rm L}_{j,i}, \omega^{\rm R}_{j,i}\f$, whether the variable or auxiliary \f$\xi_{1+i}, i=1\ldots n_x+k\f$ decomposes the monomial \f$m_{j},\ j = 1\ldots n_m\f$
.

Continuing the illustrative example in the previous section, a minimal quadratic form is computed by calling the method <a>mc::SQuad::optimize</a>:
\code
      SQF.optimize();
      std::cout << "\nSparse quadratic forms:\n" << SQF;
\endcode
By default, this optimization may be warm-started with the heuristic decomposition obtained from the method <a>mc::SQuad::process</a> implementing Algorithm 1. In particular, the redundant constraints are ignored. The final line displays the minimal quadratic forms:
\verbatim
Sparse quadratic forms:

  9 Monomials: [ 1 [0] [1] [2] [0]^2 [1]^2 [2]^2 [0]·[1]^2 [1]^4 ]

  Quadratic form for P[0]:
    1.00000e+00    ( [0] ; [0]^2 )
    1.20000e+01    ( [0] ; [2]^2 )
    3.00000e+00    ( [0] ; [0]·[1]^2 )
    3.00000e+00    ( [0] ; [1]^4 )
   -6.00000e+00    ( [2] ; [0]^2 )
   -8.00000e+00    ( [2] ; [2]^2 )
   -1.20000e+01    ( [2] ; [0]·[1]^2 )
   -6.00000e+00    ( [2] ; [1]^4 )
    1.20000e+01    ( [1]^2 ; [2]^2 )
    1.00000e+00    ( [1]^2 ; [1]^4 )

  Quadratic form for P[1]:
   -1.00000e+00    ( 1 ; 1 )
    2.00000e+00    ( [1] ; [1] )

  Auxiliary quadratic form #1:
    1.00000e+00    ( 1 ; [1]^2 )
   -1.00000e+00    ( [1] ; [1] )

  Auxiliary quadratic form #2:
    1.00000e+00    ( 1 ; [0]·[1]^2 )
   -1.00000e+00    ( [0] ; [1]^2 )

  Auxiliary quadratic form #3:
    1.00000e+00    ( 1 ; [2]^2 )
   -1.00000e+00    ( [2] ; [2] )

  Auxiliary quadratic form #4:
    1.00000e+00    ( 1 ; [1]^4 )
   -1.00000e+00    ( [1]^2 ; [1]^2 )

  Auxiliary quadratic form #5:
    1.00000e+00    ( 1 ; [0]^2 )
   -1.00000e+00    ( [0] ; [0] )
\endverbatim

These results, therfore, show that a minimum of 5 auxiliary variables only are needed to decompose the multivariate polynomials into quadratic forms: \f$\xi_4:=x_0^2, \xi_5:=x_1^2, \xi_6:=x_2^2, \xi_7:=x_0\cdot x_1^2, \xi_8:=x_1^4\f$. The resulting decomposition is given by:
\f{align*}
  q_0(\xi_0,\ldots,\xi_8) =\ & \xi_1\cdot\xi_4 + 12 \xi_1\cdot\xi_6 + 3 \xi_1\cdot\xi_7 + 3 \xi_1\cdot\xi_8 - 6 \xi_3\cdot\xi_4 - 8 \xi_3\cdot\xi_6 - 12 \xi_3\cdot\xi_7 - 6 \xi_3\cdot\xi_8 + 12 \xi_3 \xi_7 + 12 \xi_5\cdot\xi_6 + \xi_5\cdot\xi_8\\
  q_1(\xi_0,\ldots,\xi_8) =\ & - \xi_0^2 + 2 \xi_2^2
\f}
with the following 8 auxiliary quadratic constraints are also generated:
\f{align*}
  \xi_0\cdot\xi_4 - \xi_1^2 =\ & 0\\
  \xi_0\cdot\xi_5 - \xi_2^2 =\ & 0\\
  \xi_0\cdot\xi_6 - \xi_3^2 =\ & 0\\
  \xi_0\cdot\xi_7 - \xi_1\cdot\xi_5 =\ & 0\\
  \xi_0\cdot\xi_8 - \xi_5^2 =\ & 0\\
\f}

Note that this minimal decomposition is currently implemented for polynomials expressed in monomial basis only.


\section sec_SQUAD_opt What are the options in mc::SQuad and how do I set them?

The public static class member mc::SQuad::options that can be used to set/modify the options; e.g.,

\code
      mc::SQuad::options.BASIS = mc::SQuad::Options::CHEB;
      mc::SQuad::options.ORDER = mc::SQuad::Options::DEC;
\endcode

The available options are the following:

<TABLE border="1">
 <TR><TH><b>Name</b>  <TD><b>Type</b> <TD><b>Default</b> <TD><b>Description</b>
 <TR><TD><tt>mc::SQuad::Options::BASIS</tt> <TD><tt>int</tt> <TD>mc::SQuad::Options::MONOM <TD>Basis representation of the multivariate polynomial - mc::SQuad::Options::MONOM: monomial basis, mc::SQuad::Options::CHEB: Chebyshev basis
 <TR><TD><tt>mc::SQuad::Options::ORDER</tt> <TD><tt>int</tt> <TD>mc::SQuad::Options::DEC <TD>Processing order for the monomial terms - mc::SQuad::Options::INC: increasing grlex monomial order, mc::SQuad::Options::DEC: decreasing grlex monomial order
 <TR><TD><tt>mc::SQuad::Options::REDUC</tt> <TD><tt>bool</tt> <TD>false <TD>Whether to search for and append extra reduction constraints
 <TR><TD><tt>mc::SQuad::Options::CHKTOL</tt> <TD><tt>double</tt> <TD>1e-10 <TD>Tolerance for checking exactness of quadratic forms
 <TR><TD><tt>mc::SQuad::Options::LPALGO</tt> <TD><tt>int</tt> <TD>-1 <TD>LP algorithm used by MIP solver
 <TR><TD><tt>mc::SQuad::Options::LPPRESOLVE</tt> <TD><tt>int</tt> <TD>-1 <TD>LP presolve strategy in MIP solver
 <TR><TD><tt>mc::SQuad::Options::LPFEASTOL</tt> <TD><tt>double</tt> <TD>1e-9 <TD>Tolerance on LP feasibility in MIP solver
 <TR><TD><tt>mc::SQuad::Options::LPOPTIMTOL</tt> <TD><tt>double</tt> <TD>1e-9 <TD>Tolerance on LP optimality in MIP solver
 <TR><TD><tt>mc::SQuad::Options::MIPRELGAP</tt> <TD><tt>double</tt> <TD>1e-7 <TD>Tolerance on relative gap in MIP solver
 <TR><TD><tt>mc::SQuad::Options::MIPABSGAP</tt> <TD><tt>double</tt> <TD>1e-7 <TD>Tolerance on absolute gap in MIP solver
 <TR><TD><tt>mc::SQuad::Options::MIPTHREADS</tt> <TD><tt>int</tt> <TD>0 <TD>Number of threads used by MIP solver
 <TR><TD><tt>mc::SQuad::Options::MIPCONCURRENT</tt> <TD><tt>int</tt> <TD>1 <TD>Number of independent MIP solves in parallel
 <TR><TD><tt>mc::SQuad::Options::MIPFOCUS</tt> <TD><tt>int</tt> <TD>0 <TD>MIP high-level solution strategy
 <TR><TD><tt>mc::SQuad::Options::MIPHEURISTICS</tt> <TD><tt>double</tt> <TD>0.2 <TD>Fraction of time spent in MIP heuristics
 <TR><TD><tt>mc::SQuad::Options::MIPDISPLEVEL</tt> <TD><tt>int</tt> <TD>1 <TD>Display level for MIP solver
 <TR><TD><tt>mc::SQuad::Options::MIPOUTPUTFILE</tt> <TD><tt>std::string</tt> <TD>"" <TD>Name of output file for MIP model
 <TR><TD><tt>mc::SQuad::Options::MIPTIMELIMIT</tt> <TD><tt>double</tt> <TD>600 <TD>Maximum MIP runtime (seconds)
 <TR><TD><tt>mc::SQuad::Options::DISPLEN</tt> <TD><tt>unsigned int</tt> <TD>5 <TD>Number of digits in output stream
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
- <A href="https://en.wikipedia.org/w/index.php?title=Sum-of-squares_optimization&oldid=929488354">Sum-of-squares optimization</A>, <i>Wikipedia</i>, accessed: 27-Jan-2020
.
*/

#ifndef MC__SQUAD_H
#define MC__SQUAD_H

#include <list>
#include <tuple>
#include "spoly.hpp"
#include "mclapack.hpp"

//#define MC__USE_GUROBI
#if defined(MC__USE_GUROBI)
 #include "gurobi_c++.h"
 extern "C"{
  #include <fenv.h>
  int fedisableexcept( int );
 }
#endif

#define MC__SQUAD_CHECK
#undef  MC__SQUAD_PROCESS_DEBUG

namespace mc
{
//! @brief C++ structure for ordering of monomial pairs
template <typename COMP=std::less<unsigned>>
struct lt_SQuad
{
  // Comparison operator
  template <typename KEY>
  bool operator
    ()
    ( std::pair< SMon<KEY,COMP> const*, SMon<KEY,COMP> const* > const& pMon1,
      std::pair< SMon<KEY,COMP> const*, SMon<KEY,COMP> const* > const& pMon2 )
    const
    {
      // Order based on first monomial first
      if( lt_SMon<COMP>()( *pMon1.first, *pMon2.first ) ) return true;
      if( lt_SMon<COMP>()( *pMon2.first, *pMon1.first ) ) return false;
      // Order based on second monomial next
      if( lt_SMon<COMP>()( *pMon1.second, *pMon2.second ) ) return true;
      if( lt_SMon<COMP>()( *pMon2.second, *pMon1.second ) ) return false;
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
class SQuad
////////////////////////////////////////////////////////////////////////
{
  template <typename K, typename C> friend std::ostream& operator<< ( std::ostream&, SQuad<K,C> const& );
 
public:

  typedef SMon<KEY,COMP> t_SMon;
  typedef SPoly<KEY,COMP> t_SPoly;
  typedef std::set< KEY, COMP > set_SVar;
  typedef std::set< t_SMon, lt_SMon<COMP> > set_SMon;
  typedef std::set< t_SMon const*, lt_pSMon<COMP> > set_pSMon;
  typedef std::map< t_SMon, double, lt_SMon<COMP> > map_SPoly;
  typedef std::pair< t_SMon const*, t_SMon const* > key_SQuad;
  typedef std::set< key_SQuad, lt_SQuad<COMP> > set_SQuad;
  typedef std::map< key_SQuad, double, lt_SQuad<COMP> > map_SQuad;
  typedef std::map< unsigned, t_SPoly > vecmap_SPoly;

  //! @brief Options of mc::SQuad
  static struct Options
  {
    //! @brief Constructor
    Options():
      BASIS(MONOM), ORDER(DEC), REDUC(false), CHKTOL(1e-10),
#if defined(MC__USE_GUROBI)
      LPALGO( LPALGO_DEFAULT ), LPPRESOLVE(-1),
      LPFEASTOL(1e-9), LPOPTIMTOL(1e-9), MIPRELGAP(1e-7), MIPABSGAP(1e-7),
      MIPTHREADS(0), MIPCONCURRENT(1), MIPFOCUS(0), MIPHEURISTICS(0.2),
      MIPDISPLEVEL(1), MIPOUTPUTFILE(""), MIPTIMELIMIT(600),
#endif
      DISPLEN(5)
      {}
    //! @brief Assignment of mc::SQuad::Options
    Options& operator=
      ( Options& opt ){
        BASIS          = opt.BASIS;
        ORDER          = opt.ORDER;
        REDUC          = opt.REDUC;
        CHKTOL         = opt.CHKTOL;
#if defined(MC__USE_GUROBI)
        LPALGO         = opt.LPALGO;
        LPPRESOLVE     = opt.LPPRESOLVE;
        LPFEASTOL      = opt.LPFEASTOL;
        LPOPTIMTOL     = opt.LPOPTIMTOL;
        MIPRELGAP      = opt.MIPRELGAP;
        MIPABSGAP      = opt.MIPABSGAP;
        MIPTHREADS     = opt.MIPTHREADS;
        MIPCONCURRENT  = opt.MIPCONCURRENT;
        MIPFOCUS       = opt.MIPFOCUS;
        MIPHEURISTICS  = opt.MIPHEURISTICS;
        MIPDISPLEVEL   = opt.MIPDISPLEVEL;
        MIPOUTPUTFILE  = opt.MIPOUTPUTFILE;
        MIPTIMELIMIT   = opt.MIPTIMELIMIT;
#endif
        DISPLEN        = opt.DISPLEN;
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
    //! @brief Tolerance for checking exactness of quadratic forms
    double CHKTOL;
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
    //! @brief Display level for MIP solver
    int MIPDISPLEVEL;
    //! @brief Name of output file for MIP model
    std::string MIPOUTPUTFILE;
    //! @brief Maximum MIP run time (seconds)
    double MIPTIMELIMIT;
#endif
    //! @brief Number of digits in output stream for sparse polynomial coefficients
    unsigned DISPLEN;

  //! @brief Default option for LP solver
#if defined(MC__USE_GUROBI)
    static const int LPALGO_DEFAULT = -1;
#endif
  } options;
    
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
  set_SMon _SetMon;

  //! @brief Vector of sparse coefficient matrices defining the main quadratic forms
  std::vector<map_SQuad> _MatFct;

  //! @brief Vector of sparse coefficient matrices defining the auxiliary quadratic forms
  std::vector<map_SQuad> _MatRed;

  //! @brief Vector of sparse coefficient matrices defining the postiive semi-defninite cuts
  std::vector<map_SQuad> _MatPSD;
  
#if defined(MC__USE_GUROBI)
  //! @brief whether the MIP solver has sent an exception
  bool _MIPexcpt;
  //! @brief Set of monomials in MIP optimization model
  set_SMon _MIP_SetMon;
  //! @brief Dummy monomial storing partipating variables and their higest degrees in MIP optimization model
  t_SMon _MIP_VarDeg;
  //! @brief MIP environment
  GRBEnv* _GRBenv;
  //! @brief MIP model
  GRBModel* _GRBmodel;
  //! @brief vector of binary variables indicating active auxiliary variables
  std::vector<GRBVar> _MIP_auxbin;
  //! @brief vector of continuous variables describing the monomial structure of auxiliary variables
  std::vector<std::vector<GRBVar>> _MIP_auxexp;
  //! @brief vectors of binary variables describing the quadratic decomposition of auxiliary variables
  std::vector<std::vector<GRBVar>> _MIP_auxdec1;
  std::vector<std::vector<GRBVar>> _MIP_auxdec2;
  //! @brief vectors of binary variables describing the quadratic decomposition of monomials
  std::vector<std::vector<GRBVar>> _MIP_mondec1;
  std::vector<std::vector<GRBVar>> _MIP_mondec2;
#endif

public:

  //! @brief Default Constructor
  SQuad
    ()
    {
#if defined(MC__USE_GUROBI)
      _GRBenv   = new GRBEnv();
      _GRBmodel = nullptr;
#endif
    }

  //! @brief Destructor
  virtual ~SQuad
    ()
    {
      _reset();
#if defined(MC__USE_GUROBI)
      delete _GRBmodel;
      delete _GRBenv;
#endif
    }
  
  //! @brief Process the sparse polynomials in array <a>pPol</a> indexed by <a>ndxSPol</a>
  template <typename POL>
  double process
    ( std::set<unsigned> const& ndxSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
      int const BASIS=options.BASIS, bool const CHECK=false );

  //! @brief Process the <a>nPol</a> sparse polynomials in array <a>pPol</a>
  template <typename POL>
  double process
    ( unsigned const nSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
      int const BASIS=options.BASIS, bool const CHECK=false );

  //! @brief Process the sparse polynomial <a>Pol</a>
  template <typename POL>
  double process
    ( POL const& SPol, map_SPoly const& (POL::*mapmon)() const,
      int const BASIS=options.BASIS, bool const CHECK=false );

  //! @brief Optimize the quadratic expressions for minimum number of auxiliaries 
  void optimize
    ( bool const warmStart=true );

  //! @brief Decompose the quadratic expression <a>mat</a> into separable expressions
  std::list< map_SQuad > separate
    ( map_SQuad const& mat )
    const;

  //! @brief Factorize the quadratic expression <a>mat</a> using eigenvalue decomposition
  std::multimap< double, map_SPoly > factorize
    ( map_SQuad const& mat )
    const;

  //! @brief Generate positive semi-definite cuts to tighten the quadratic reformulation
  void tighten
    ( bool const threevar=false );

  //! @brief Check quadratic form expressions by comparing with Chebyshev model
  template <typename POL>
  double check
    ( unsigned const nSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
      int const BASIS )
    const
    { return _check( nSPol, pSPol, mapmon, BASIS ); }

  //! @brief Retreive reference to vector of sparse coefficient matrices defining the main quadratic forms
  std::vector<map_SQuad> const& MatFct
    ()
    const
    { return _MatFct; }

  //! @brief Retreive reference to vector of sparse coefficient matrices defining the auxiliary quadratic forms
  std::vector<map_SQuad> const& MatRed
    ()
    const
    { return _MatRed; }

  //! @brief Retreive reference to vector of sparse coefficient matrices defining the positive semi-definite cuts
  std::vector<map_SQuad> const& MatPSD
    ()
    const
    { return _MatPSD; }

  //! @brief Retreive reference to set of monomials in quadratic forms
  set_SMon const& SetMon
    ()
    const
    { return _SetMon; }

  //! @brief Reset quadratic form expressions
  void reset
    ()
    { _reset(); }

protected:

  //! @brief Reorder entries in a monomial pair
  key_SQuad& _reorder
    ( key_SQuad&& pMon )
    const;

  //! @brief Create set of monomials for Chebyshev product between <a>mon1</a> and <a>mon2</a>
  set_SMon _prodmon
    ( t_SMon const& mon1, t_SMon const& mon2 )
    const;

  //! @brief Insert new entry in matrix <a>mat</a> corresponding to the mononial <a>pMon</a> in a product with the constant monomial and with corresponding coefficient <a>coef</a>
  bool _insert
    ( map_SQuad& mat, t_SMon const* pMon, double const coef, bool const add=false );//, bool const rem=false );

  //! @brief Insert new entry in matrix <a>mat</a> corresponding to the product between mononials <a>pLMon</a> and <a>pRMon</a> with corresponding coefficient <a>coef</a>
  bool _insert
    ( map_SQuad& mat, t_SMon const* pLMon, t_SMon const* pRMon,
      double const coef, bool const add=false );

  //! @brief Insert new entry in matrix <a>mat</a> corresponding to the product between mononials <a>pLMon</a> and <a>pRMon</a> with corresponding coefficient <a>coef</a>, and insert the associated low-order terms in monomial map <a>mapmon</a>
  bool _insert
    ( map_SQuad& mat, map_SPoly& mapmon, t_SMon const* pLMon,
      t_SMon const* pRMon, double const coef,
      bool const add=false );
      
  //! @brief Decompose <a>mon</a> into product of two monomials in <a>_SetMon</a> by expanding <a>_SetMon</a> as necessary
  key_SQuad _decompose
    ( t_SMon const& mon );

  //! @brief Search for <a>mon</a> in <a>_SetMon</a> and append it to <a>_SetMon</a> if missing along reduction constraints in <a>_MatRed</a>
  typename set_SMon::iterator _subexpression
    ( t_SMon const& mon );

  //! @brief Search for extra reduction constraints for <a>mon</a> and append them to <a>_MatRed</a>
  void _reduction
    ( t_SMon const& mon );

  //! @brief Check quadratic form expressions by comparing with Chebyshev model
  template <typename POL>
  double _check
    ( unsigned const nSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
      int const BASIS )
    const;

  //! @brief Check if a term exist in current quadratic forms 
  bool _find
    ( t_SMon const* pMon1, t_SMon const* pMon2 ) 
    const;

  //! @brief Check possible decompositions using existing monomial in _SetMon
  void _candidates
    ( set_SQuad& CandidateDec, t_SMon const& mon )
    const;

  //! @brief Check connections between a quadratic entry and a quadratic form
  bool _isconnected
    ( key_SQuad const& entry, map_SQuad const& mat )
    const;

  //! @brief Check if monomial maps are all empty
  bool _empty
    ( vecmap_SPoly const& vecSPol )
    const;

  //! @brief Select next monomial to be processed in vector of monomial maps
  std::tuple< unsigned, t_SMon, double > _next
    ( vecmap_SPoly const& vecSPol )
    const;

  //! @brief Encode MIP optimization model for minimal decomposition
  void _MIP_encode
    ( unsigned const maxAux, unsigned const minOrd );

  //! @brief Encode MIP optimization model for minimal decomposition
  void _MIP_decode
    ( unsigned const minOrd );

  //! @brief Reset variable vectors in MIP optimization model
  void _MIP_reset
    ();

  //! @brief Set options in MIP optimization model
  void _MIP_options
    ();

  //! @brief Initialize MIP optimization model with current decomposition
  void _MIP_initialize
    ( unsigned const minOrd );

  //! @brief Solve MIP optimization model
  void _MIP_solve
    ();

  //! @brief Display MIP current point
  void _MIP_display
    ( std::ostream& os=std::cout )
    const;

private:

  //! @brief Reset the quadratic form
  void _reset
    ();
};

template <typename KEY, typename COMP>
inline typename SQuad<KEY,COMP>::Options SQuad<KEY,COMP>::options;

#if defined(MC__USE_GUROBI)
template <typename KEY, typename COMP>
inline int const SQuad<KEY,COMP>::Options::LPALGO_DEFAULT;
#endif

////////////////////////////////////////////////////////////////////////

template <typename KEY, typename COMP>
inline std::ostream&
operator<<
( std::ostream& out, SQuad<KEY,COMP> const& quad )
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

template <typename KEY, typename COMP>
inline bool
SQuad<KEY,COMP>::_find
( t_SMon const* pMon1, t_SMon const* pMon2 ) 
const
{
  for( auto const& mat : _MatFct )
    if( mat.count( std::make_pair( pMon1, pMon2 ) ) ) return true;
  for( auto const& mat : _MatRed )
    if( mat.count( std::make_pair( pMon1, pMon2 ) ) ) return true;
  return false;
}

template <typename KEY, typename COMP>
inline void
SQuad<KEY,COMP>::tighten
( bool const threevar )
{
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
      _MatPSD.push_back( map_SQuad() );
      auto& mat1 = _MatPSD.back();
      assert( _insert( mat1, &*itmon1, &*itmon1,  1. )
           && _insert( mat1, &*itmon2, &*itmon2,  1. )
           && _insert( mat1, &*itmon1, &*itmon2, -2. ) );
      if( !threevar && options.BASIS == Options::MONOM
       && !itmon1->gcexp()%2 && !itmon2->gcexp()%2 ) continue;
      _MatPSD.push_back( map_SQuad() );
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

template <typename KEY, typename COMP>
inline std::list<typename SQuad<KEY,COMP>::map_SQuad>
SQuad<KEY,COMP>::separate
( map_SQuad const& mat )
const
{
  std::list<map_SQuad> submat;
  for( auto itmat : mat ){
    key_SQuad ijmon = itmat.first;
    double coef = itmat.second;
    auto itsubmat = submat.end();
//    // Search for alternative decompositions for terms multiplying monomial '1'
//    // NEED CORRECTION IN ORDER TO ACCOUNT THAT DECOMPOSED TERM MAY ALREADY BE PRESENT IN MAT!
//    if( options.BASIS == Options::MONOM && !ijmon.first->tord && ijmon.second->tord ){
//      set_SQuad CandidateDec;
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
    submat.push_back( map_SQuad() );
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

template <typename KEY, typename COMP>
inline bool
SQuad<KEY,COMP>::_isconnected
( key_SQuad const& entry, map_SQuad const& mat )
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

template <typename KEY, typename COMP>
inline std::multimap<double,typename SQuad<KEY,COMP>::map_SPoly>
SQuad<KEY,COMP>::factorize
( map_SQuad const& mat )
const
{
  // Populate sparse symmetric coefficient matrix
  CPPL::dssmatrix coefmat;
  map_SPoly indexmap;
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
  std::multimap<double,map_SPoly> eigdec;
  std::vector<double> eigval;
  std::vector<CPPL::dcovector> eigvec;
  CPPL::dsymatrix dcoefmat = coefmat.to_dsymatrix();
  if( dcoefmat.dsyev( eigval, eigvec ) ) return eigdec;

  // Populate decomposition map
  auto itval = eigval.begin();
  auto itvec = eigvec.begin();
  for( ; itval != eigval.end(); ++itval, ++itvec ){
    map_SPoly eigterm;
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

template <typename KEY, typename COMP>
template< typename POL >
inline double
SQuad<KEY,COMP>::_check
( unsigned const nSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
  int const BASIS )
const
{
  assert( nSPol <= _MatFct.size() );
  double sumdiff = 0e0;
#ifdef MC__SQUAD_DEBUG_CHECK
  std::cout << *this << std::endl;
#endif
  
  // Process entries in _SetMon
  std::map< t_SMon, t_SPoly, lt_SMon<COMP> > mmon; 
  for( auto const& mon : _SetMon ){
    mmon[mon] = std::make_pair( mon, 1e0 );
#ifdef MC__SQUAD_DEBUG_CHECK
    std::cout << mmon[mon];
#endif
  }

  // Check entries in _MatFct
  auto itmat = _MatFct.begin();
  std::advance( itmat, _MatFct.size()-nSPol );
  for( unsigned i=0; nSPol && itmat != _MatFct.end(); ++itmat, ++i ){
    SPoly<KEY,COMP>::options.BASIS = BASIS;
    t_SPoly spe( (pSPol[i].*mapmon)() ); spe.convert( options.BASIS );
    SPoly<KEY,COMP>::options.BASIS = options.BASIS;
    for( auto const& [ijmon,coef] : *itmat )
      spe -= coef * mmon[*ijmon.first] * mmon[*ijmon.second];
#ifdef MC__SQUAD_DEBUG_CHECK
    std::cout << "\n  Quadratic form of P[" << _MatFct.size()-nSPol+i << "]:" << spe;
#endif
    double locdiff = 0;
    for( auto const& [mon,coef] : spe.mapmon() )
      locdiff += std::fabs( coef ); 
    if( std::fabs(locdiff) > options.CHKTOL )
      std::cerr << "\n  Error in quadratic form of P[" << _MatFct.size()-nSPol+i << "]:" << spe;
    sumdiff += locdiff;  
  }

  // Check entries in _MatRed
  itmat = _MatRed.begin();
  for( unsigned i=0; itmat != _MatRed.end(); ++itmat, ++i ){
    t_SPoly spe;
    for( auto const& [ijmon,coef] : *itmat ){
      spe += coef * mmon[*ijmon.first] * mmon[*ijmon.second];
    }
#ifdef MC__SQUAD_DEBUG_CHECK
    std::cout << " \n Auxiliary quadratic form #" << i << spe;
#endif
    double locdiff = 0;
    for( auto const& [mon,coef] : spe.mapmon() )
      locdiff += std::fabs( coef ); 
    if( std::fabs(locdiff) > options.CHKTOL )
      std::cerr << " \n Error in auxiliary quadratic form #" << i << spe;
    sumdiff += locdiff;  
  }

  // Check entries in _MatPSD
  itmat = _MatPSD.begin();
  for( unsigned i=0; itmat != _MatPSD.end(); ++itmat, ++i ){
    t_SPoly spe;
    for( auto const& [ijmon,coef] : *itmat ){
      spe += coef * mmon[*ijmon.first] * mmon[*ijmon.second];
    }
#ifdef MC__SQUAD_DEBUG_CHECK
    std::cout << " \n Positive semi-definite quadratic form #" << i << spe;
#endif
    double locdiff = 0;
    for( auto const& [mon,coef] : spe.mapmon() )
      locdiff += std::fabs( coef ); 
    if( std::fabs(locdiff) > options.CHKTOL )
      std::cerr << " \n Error in auxiliary quadratic form #" << i << spe;
    sumdiff += locdiff;  
  }

  return sumdiff;
}

template <typename KEY, typename COMP>
template <typename POL >
inline double
SQuad<KEY,COMP>::process
( std::set<unsigned> const& ndxSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
  int const BASIS, bool const CHK )
{
  std::vector<POL> vpSPol;
  vpSPol.reserve( ndxSPol.size() );
  for( unsigned const& i : ndxSPol ) vpSPol.push_back( pSPol[i] );
  return process( ndxSPol.size(), vpSPol.data(), mapmon, BASIS, CHK );
}

template <typename KEY, typename COMP>
template <typename POL >
inline double
SQuad<KEY,COMP>::process
( POL const& SPol, map_SPoly const& (POL::*mapmon)() const,
  int const BASIS, bool const CHK )
{
  return process( 1, &SPol, mapmon, BASIS, CHK );
}

template <typename KEY, typename COMP>
inline bool
SQuad<KEY,COMP>::_empty
( vecmap_SPoly const& vecSPol )
const
{
  for( auto const& [ndx,SPol] : vecSPol )
    if( !SPol.mapmon().empty() ) return false;
  return true;
}

template <typename KEY, typename COMP>
inline std::tuple< unsigned, typename SQuad<KEY,COMP>::t_SMon, double >
SQuad<KEY,COMP>::_next
( vecmap_SPoly const& vecSPol )
const
{
  // Locate first non-empty monomial map
  auto itnext = vecSPol.cbegin();
  for( ; itnext != vecSPol.cend(); ++itnext ){
    auto const& mapnext = itnext->second.mapmon();
    if( !mapnext.empty() ) break;
  }

  // Locate monomial map with next element
  auto ittry = itnext;
  for( ++ittry; ittry != vecSPol.cend(); ++ittry ){
    auto const& mapnext = itnext->second.mapmon();
    auto const& maptry  = ittry->second.mapmon();
    //t_SMon const& monnext = mapnext.begin()->first;
    if( maptry.empty() ) continue;
    switch( options.ORDER ){
      case Options::INC:
#ifdef MC__SQUAD_DEBUG_NEXT
        std::cout << "SQuad::_next: Comparing " << maptry.cbegin()->first.display(options.BASIS)
                  << " <? " << mapnext.cbegin()->first.display(options.BASIS) << std::endl;
#endif
        if( lt_SMon<COMP>()( maptry.cbegin()->first, mapnext.cbegin()->first ) ) itnext = ittry;
        break;
      case Options::DEC:
#ifdef MC__SQUAD_DEBUG_NEXT
        std::cout << "SQuad::_next: Comparing " << mapnext.crbegin()->first.display(options.BASIS)
                  << " <? " << maptry.crbegin()->first.display(options.BASIS) << std::endl;
#endif
        if( lt_SMon<COMP>()( mapnext.crbegin()->first, maptry.crbegin()->first ) ) itnext = ittry;
        break;
    }
  }
  
  // Return next element
  assert( itnext != vecSPol.end() );
  auto const& mapnext = itnext->second.mapmon();
#ifdef MC__SQUAD_DEBUG_NEXT
  std::cout << "SQuad::_next: Selected " << mapnext.crbegin()->first.display(options.BASIS) << std::endl;
#endif
  switch( options.ORDER ){
    default:
    case Options::INC:
      return std::make_tuple( itnext->first, mapnext.cbegin()->first, mapnext.cbegin()->second );
    case Options::DEC:
      return std::make_tuple( itnext->first, mapnext.crbegin()->first, mapnext.crbegin()->second );
  }
}

template <typename KEY, typename COMP>
template <typename POL>
inline double
SQuad<KEY,COMP>::process
( unsigned const nSPol, POL const* pSPol, map_SPoly const& (POL::*mapmon)() const,
  int const BASIS, bool const CHK )
{
  if( !nSPol ) return( CHK? _check( nSPol, pSPol, mapmon, BASIS ): 0. ); 

  // Initialize monomial vector with constant monomial and participating variables
  _SetMon.insert( t_SMon() );

  // Build monomial multimap
  SPoly<KEY,COMP>::options.BASIS = BASIS;
  vecmap_SPoly vecSPol;
  for( unsigned i=0; i<nSPol; i++ ){
    // Add new quadratic form and define corresponding pointer pmat
    unsigned ndxmat = _MatFct.size();
    _MatFct.push_back( map_SQuad() );

    // Make local copy and convert to desired basis
    vecSPol[ndxmat] = (pSPol[i].*mapmon)(); vecSPol[ndxmat].convert( options.BASIS );

    // Insert all participating variables in monomial set
    _SetMon.insert( t_SMon() );
    for( auto var : vecSPol[ndxmat].setvar() )
      _SetMon.insert( t_SMon( var ) );
  }
  SPoly<KEY,COMP>::options.BASIS = options.BASIS;
  //unsigned const defVar = _SetMon.size();

  // Iterate through monomial terms
  for( ; !_empty( vecSPol ); ){
    // Local copy of next monomial, then erase
    auto const [ ndxmat, mon, coef ] = _next( vecSPol );
    vecSPol[ndxmat].mapmon().erase( mon );
    auto& mat = _MatFct[ndxmat];

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
    }

    // Mononial needs further decomposition
    else{
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
      bool ins = _insert( mat, vecSPol[ndxmat].mapmon(), plmon, prmon, coef );
#ifdef MC__SQUAD_CHECK
      assert( ins );
#endif
    }
  }

  return( CHK? _check( nSPol, pSPol, mapmon, BASIS ): 0. );
}

#if defined(MC__USE_GUROBI)
template <typename KEY, typename COMP>
inline void
SQuad<KEY,COMP>::optimize
( bool const warmStart )
{
  if( options.BASIS != Options::MONOM ) return;

  // Current number of auxiliaries
  unsigned maxAux = 0;
  for( auto const& mon : _SetMon ){
    if( mon.tord < 2 ) continue;
    maxAux++;
  }
  
  // Only working with minOrd = 2 currently
  unsigned const minOrd = 2;

  // Run MIP optimization for a minimal representation
  _MIP_encode( maxAux, minOrd );
  if( warmStart) _MIP_initialize( minOrd );
  _MIP_solve();
  _MIP_decode( minOrd );
}

template <typename KEY, typename COMP>
inline void
SQuad<KEY,COMP>::_MIP_encode
( unsigned const maxAux, unsigned const minOrd )
{
  // Reset variable vectors
  _MIP_reset();

  // Construct a set of monomials of degree minOrd or greater, and corresponding set of participating variables
  _MIP_SetMon.clear();
  _MIP_VarDeg.expr.clear(); _MIP_VarDeg.tord = 0;
  for( auto const& mat : _MatFct ){
    for( auto const& [ijmon,coef] : mat ){
      if( ijmon.first->tord + ijmon.second->tord < minOrd ) continue;
      auto const& [itmon,ins] = _MIP_SetMon.insert( *ijmon.first + *ijmon.second );
      if( ins ) _MIP_VarDeg.hull( *itmon ); 
    }
  }
  unsigned const maxDeg = _MIP_SetMon.rbegin()->tord;

  unsigned const nVar = _MIP_VarDeg.expr.size(); // Number of participating variables
  _MIP_auxbin.reserve( maxAux );   // z[k]: whether auxiliary x[nVar+k], 0<=k<maxAux, used in decomposition
  _MIP_auxexp.reserve( maxAux );   // b[k,i]: power of variable x[i], 0<=i<nVar, in auxiliary x[nVar+k], 0<=k<maxAux
  _MIP_auxdec1.reserve( maxAux );  // v1|2[k,i]: whether variable/auxiliary x[i], 0<=i<nVar+k, in decomposition of 
  _MIP_auxdec2.reserve( maxAux );  //            auxiliary x[nVar+k], 0<=k<maxAux

  for( unsigned k=0; k<maxAux; ++k ){
#if defined(MC__USE_GUROBI)
    // z[k] in {0,1} with corresponding cost coefficient c[k] = 1
    _MIP_auxbin.push_back( _GRBmodel->addVar( 0., 1., 1., GRB_BINARY ) );
    // z[k] <= z[k-1]
    if( k ) _GRBmodel->addConstr( _MIP_auxbin[k-1], GRB_GREATER_EQUAL, _MIP_auxbin[k] );

    // b[k,i] in [0,maxDeg(Var[i])], 0<=i<nVar
    // b[k,i] <= z[k] * maxDeg(Var[i]), 0<=i<nVar
    _MIP_auxexp.push_back( std::vector<GRBVar>() );
    _MIP_auxexp[k].reserve( nVar );
    for( auto const& [Var,maxExp] : _MIP_VarDeg.expr ){
      _MIP_auxexp[k].push_back( _GRBmodel->addVar( 0., maxExp, 0., GRB_CONTINUOUS ) );
      _GRBmodel->addConstr( _MIP_auxexp[k].back(), GRB_LESS_EQUAL, maxExp * _MIP_auxbin[k] );
    }
    
    // v1[k,i], v2[k,i] in {0,1}, 0<=i<nVar+k
    _MIP_auxdec1.push_back( std::vector<GRBVar>() );
    _MIP_auxdec2.push_back( std::vector<GRBVar>() );
    _MIP_auxdec1[k].reserve( nVar+k );
    _MIP_auxdec2[k].reserve( nVar+k );
    // SUM( 0<=i<nVar+k, v1[k,i] ) = z[k]
    // SUM( 0<=i<nVar+k, v2[k,i] ) = z[k]
    // SUM( 0<=i<j, v1[k,i] ) >= SUM( 0<=i<j, v2[k,i] ), 0<=j<nVar+k
    GRBLinExpr sum_auxdec1, sum_auxdec2;
    for( unsigned i=0; i<nVar+k; ++i ){
      _MIP_auxdec1[k].push_back( _GRBmodel->addVar( 0., 1., 0., GRB_BINARY ) );
      _MIP_auxdec2[k].push_back( _GRBmodel->addVar( 0., 1., 0., GRB_BINARY ) );
      sum_auxdec1 += _MIP_auxdec1[k][i];
      sum_auxdec2 += _MIP_auxdec2[k][i];
      //_GRBmodel->addConstr( sum_auxdec1, GRB_GREATER_EQUAL, sum_auxdec2 ); // symmetry breaking
    }
    _GRBmodel->addConstr( sum_auxdec1, GRB_EQUAL, _MIP_auxbin[k] );//1. );
    _GRBmodel->addConstr( sum_auxdec2, GRB_EQUAL, _MIP_auxbin[k] );//1. );

    // b[k,i] = v1[k,i]+v2[k,i] + SUM( 0<=l<k, (v1[k,nVar+l]+v2[k,nVar+l])*b[l,i] ), 0<=i<nVar
    for( unsigned i=0; i<nVar; ++i ){
      GRBQuadExpr sum_auxexp;
      sum_auxexp += _MIP_auxdec1[k][i] + _MIP_auxdec2[k][i];
      for( unsigned l=0; l<k; ++l ){
        sum_auxexp.addTerm( 1., _MIP_auxdec1[k][nVar+l], _MIP_auxexp[l][i] );
        sum_auxexp.addTerm( 1., _MIP_auxdec2[k][nVar+l], _MIP_auxexp[l][i] );
      }
      _GRBmodel->addQConstr( sum_auxexp, GRB_EQUAL, _MIP_auxexp[k][i] );
    }      
#endif
  }

  unsigned const nMon = _MIP_SetMon.size(); // Number of monomials to be decomposed
  std::cout << "No monomials to decompose: " << nMon << std::endl;
  _MIP_mondec1.reserve( nMon );  // w1|2[j,i]: whether variable/auxiliary x[i], 0<=i<nVar+k, in decomposition of 
  _MIP_mondec2.reserve( nMon );  //            monomial m[j], 0<=j<nMon

  auto jmon = _MIP_SetMon.cbegin();
  for( unsigned j=0; j<nMon; ++j, ++jmon ){
#if defined(MC__USE_GUROBI)
    // w1[k,i], w2[k,i] in {0,1}, 0<=i<nMon
    _MIP_mondec1.push_back( std::vector<GRBVar>() );
    _MIP_mondec2.push_back( std::vector<GRBVar>() );
    _MIP_mondec1[j].reserve( nVar+maxAux );
    _MIP_mondec2[j].reserve( nVar+maxAux );
    // SUM( 0<=i<nVar+maxAux, w1[j,i] )  = 1
    // SUM( 0<=i<nVar+maxAux, w2[j,i] ) <= 1  could be 0 if monomial present
    // SUM( 0<=i<k, w1[j,i] ) >= SUM( 0<=i<k, w2[j,i] ), 0<=k<nVar+maxAux
    GRBLinExpr sum_mondec1, sum_mondec2;
    for( unsigned i=0; i<nVar+maxAux; ++i ){
      _MIP_mondec1[j].push_back( _GRBmodel->addVar( 0., 1., 0., GRB_BINARY ) );
      _MIP_mondec2[j].push_back( _GRBmodel->addVar( 0., 1., 0., GRB_BINARY ) );
      sum_mondec1 += _MIP_mondec1[j][i];
      sum_mondec2 += _MIP_mondec2[j][i];
      //_GRBmodel->addConstr( sum_mondec1, GRB_GREATER_EQUAL, sum_mondec2 ); // symmetry breaking
    }
    _GRBmodel->addConstr( sum_mondec1, GRB_EQUAL, 1. );
    _GRBmodel->addConstr( sum_mondec2, GRB_LESS_EQUAL, 1. );

    // b[k,i] = w1[j,i]+w2[j,i] + SUM( 0<=k<maxAux, (w1[j,nVar+k]+w2[j,nVar+k])*b[k,i] ), 0<=i<nVar
    auto imonvar = _MIP_VarDeg.expr.cbegin();
    for( unsigned i=0; i<nVar; ++i, ++imonvar ){
      GRBQuadExpr sum_monexp;
      sum_monexp += _MIP_mondec1[j][i] + _MIP_mondec2[j][i];
      for( unsigned k=0; k<maxAux; ++k ){
        sum_monexp.addTerm( 1., _MIP_mondec1[j][nVar+k], _MIP_auxexp[k][i] );
        sum_monexp.addTerm( 1., _MIP_mondec2[j][nVar+k], _MIP_auxexp[k][i] );
      }
      auto ivar = jmon->expr.find( imonvar->first );
      if( ivar != jmon->expr.cend() )
        _GRBmodel->addQConstr( sum_monexp, GRB_EQUAL, ivar->second );
      else
        _GRBmodel->addQConstr( sum_monexp, GRB_EQUAL, 0. );
    }
#endif
  }
  
  // SUM( 0<=i<nVar, b[l,i] ) <= maxDeg(Mon[j], 0<=j<nMon), 0<=l<nAux
  for( unsigned l=0; l<maxAux; ++l ){
    GRBLinExpr sum_auxexp;
    for( unsigned i=0; i<nVar; ++i )
      sum_auxexp += _MIP_auxexp[l][i];
    _GRBmodel->addConstr( sum_auxexp, GRB_LESS_EQUAL, maxDeg-1 );
  }
}

template <typename KEY, typename COMP>
inline void
SQuad<KEY,COMP>::_MIP_initialize
( unsigned const minOrd )
{
  std::map< t_SMon const*, unsigned, lt_pSMon<COMP> > mapMon;
  std::map< KEY, unsigned, COMP > mapVar;
  unsigned const nVar   = _MIP_VarDeg.expr.size();
  unsigned const maxAux = _MIP_auxbin.size();
  auto const* pone = &*_SetMon.cbegin();
  mapMon[ pone ] = nVar + maxAux;

  // Map participating variables
  unsigned i = 0;
  for( auto const& [var,expmax] : _MIP_VarDeg.expr ){
    mapVar[ var ] = i;
    auto itmon = _SetMon.find( t_SMon( var ) );
    assert( itmon != _SetMon.end() );
    mapMon[ &*itmon ] = i;
    ++i;
  }

  // Activate auxiliary binaries and initialize auxiliary exponents
  unsigned k = 0;
  for( auto const& mon : _SetMon ){
    if( mon.tord < minOrd ) continue;
    _MIP_auxbin[k].set( GRB_DoubleAttr_Start, 1. );
    for( auto& grbvar : _MIP_auxexp[k] )
      grbvar.set( GRB_DoubleAttr_Start, 0. );
    std::cout << "b[" << k << "] =";
    for( auto const& [var,ord] : mon.expr ){
      std::cout << " " << mapVar[var] << ":" << ord;
      _MIP_auxexp[k][ mapVar[var] ].set( GRB_DoubleAttr_Start, ord );
    }
    std::cout << std::endl;
    mapMon[ &mon ] = nVar + k;
    ++k;
  }

  // Initialize monomial decomposition binaries
  for( auto const& mat : _MatFct ){
    for( auto const& [ijmon,coef] : mat ){
      if( ijmon.first->tord + ijmon.second->tord < minOrd ) continue;
      auto itmon = _MIP_SetMon.find( *ijmon.first + *ijmon.second );
      assert( itmon != _MIP_SetMon.cend() );
      unsigned j = 0;
      for( auto jtmon = _MIP_SetMon.cbegin(); jtmon != itmon; ++jtmon, ++j ) continue;
      for( auto& grbvar : _MIP_mondec1[j] )
        grbvar.set( GRB_DoubleAttr_Start, 0. );
      for( auto& grbvar : _MIP_mondec2[j] )
        grbvar.set( GRB_DoubleAttr_Start, 0. );
      std::cout << "Auxiliary monomial " << itmon->display( Options::MONOM ) << ": "
                                         << ijmon.first->display( Options::MONOM ) << ", " 
                                         << ijmon.second->display( Options::MONOM ) << std::endl;
      std::cout << "Auxiliary monomial #" << j << ": " << mapMon[ijmon.first]
                                               << ", " << mapMon[ijmon.second] << std::endl;
      if( ijmon.first != pone ){
        _MIP_mondec1[j][ mapMon[ijmon.first]  ].set( GRB_DoubleAttr_Start, 1. );
        if( ijmon.second != pone )
          _MIP_mondec2[j][ mapMon[ijmon.second] ].set( GRB_DoubleAttr_Start, 1. );
      }
      else if( ijmon.second != pone )
        _MIP_mondec1[j][ mapMon[ijmon.second] ].set( GRB_DoubleAttr_Start, 1. );
    }
  }

  // Initialize auxiliary decomposition binaries
  for( auto const& mat : _MatRed ){
    if( mat.size() != 2 ) continue; // only auxiliary constraints defining new auxiliaries
    //assert( mat.size() == 2 );
    auto const& ijmon = mat.cbegin()->first;
    if( ijmon.first != pone ) continue; // only auxiliary constraints defining new auxiliaries
    //assert( ijmon.first == pone );
    if( mapMon.find( ijmon.second ) == mapMon.end() ) continue;
    unsigned k = mapMon[ijmon.second] - nVar;
    for( auto& grbvar : _MIP_auxdec1[k] )
      grbvar.set( GRB_DoubleAttr_Start, 0. );
    for( auto& grbvar : _MIP_auxdec2[k] )
      grbvar.set( GRB_DoubleAttr_Start, 0. );
    auto const& ijmon2 = mat.crbegin()->first;
    assert( ijmon2.first != pone && ijmon2.second != pone );
    _MIP_auxdec1[k][ mapMon[ijmon2.first]  ].set( GRB_DoubleAttr_Start, 1. );
    _MIP_auxdec2[k][ mapMon[ijmon2.second] ].set( GRB_DoubleAttr_Start, 1. );
  }
}

template <typename KEY, typename COMP>
inline void
SQuad<KEY,COMP>::_MIP_decode
( unsigned const minOrd )
{
  set_SMon SetMonOpt;
  std::vector<map_SQuad> MatFctOpt, MatRedOpt;
  std::vector<t_SMon const*> vecMonOpt;
  unsigned const nVar   = _MIP_VarDeg.expr.size();
  unsigned const optAux = std::round( _GRBmodel->get( GRB_DoubleAttr_ObjVal ) );
  vecMonOpt.resize( nVar+optAux+1 ); // Unit monomial stored in last position
  
  // Add auxiliary monomials to SetMonOpt 
  auto const& [itone,ins] = SetMonOpt.insert( t_SMon() );
  vecMonOpt[nVar+optAux] = &*itone;
  for( unsigned k=0; k<optAux; ++k ){
    t_SMon mon;
    unsigned i = 0;
    for( auto const& [var,expmax] : _MIP_VarDeg.expr ){
      unsigned exp = std::round( _MIP_auxexp[k][i++].get( GRB_DoubleAttr_X ) );
      if( !exp ) continue;
      mon += t_SMon( var, exp );
    }
    std::cout << "Auxiliary monomial #" << k << ": " << mon.display( Options::MONOM ) << std::endl;
    auto const& [itmon,ins] = SetMonOpt.insert( mon );
    //assert( ins ); // Suboptimal solutions may include the same monomial multiple times!
    vecMonOpt[nVar+k] = &*itmon;    
  }
  
  // Duplicate entries in MatFctOpt and SetMonOpt for lower degree monomials
  for( auto const& mat : _MatFct ){
    MatFctOpt.push_back( map_SQuad() );
    auto& matopt = MatFctOpt.back();
    for( auto const& [ijmon,coef] : mat ){
      if( ijmon.first->tord + ijmon.second->tord >= minOrd ) continue;
      auto const& [itmon1,ins1] = SetMonOpt.insert( *ijmon.first );
      auto const& [itmon2,ins2] = SetMonOpt.insert( *ijmon.second );
      _insert( matopt, &*itmon1, &*itmon2, coef );
    }
  }

  // Complete vector entries in vecMonOpt
  unsigned i = 0;
  for( auto const& [var,maxExp] : _MIP_VarDeg.expr ){
    auto const& [itmon,ins] = SetMonOpt.insert( t_SMon( var ) );
    assert( itmon != SetMonOpt.end() );
    vecMonOpt[i++] = &*itmon;
  }

  // Add entries in MatFctOpt for higher degree monomials
  auto itmatopt = MatFctOpt.begin();
  for( auto const& mat : _MatFct ){
    auto& matopt = *itmatopt;
    for( auto const& [ijmon,coef] : mat ){
      if( ijmon.first->tord + ijmon.second->tord < minOrd ) continue;
      auto itmon = _MIP_SetMon.find( *ijmon.first + *ijmon.second );
      assert( itmon != _MIP_SetMon.cend() );
      unsigned j = 0;
      for( auto jtmon = _MIP_SetMon.cbegin(); jtmon != itmon; ++jtmon, ++j )
        continue;
      unsigned i1 = 0, i2 = 0;
      for( auto i1mon = _MIP_mondec1[j].cbegin(); i1 < nVar+optAux; ++i1mon, ++i1 )
        if( i1mon->get(GRB_DoubleAttr_X) > 0.9 ) break;
      for( auto i2mon = _MIP_mondec2[j].cbegin(); i2 < nVar+optAux; ++i2mon, ++i2 )
        if( i2mon->get(GRB_DoubleAttr_X) > 0.9 ) break;
      //std::cout << "j,i1,i2 = " << j << "," << i1 << "," << i2 << " (max: " << nVar+optAux << ")" << std::endl;
      _insert( matopt, vecMonOpt[i1], vecMonOpt[i2], coef );
    }
    ++itmatopt;
  }

  // Add entries in MatRedOpt for auxiliary monomials
  for( unsigned k=0; k<optAux; ++k ){
    MatRedOpt.push_back( map_SQuad() );
    auto& matopt = MatRedOpt.back();
    _insert( matopt, vecMonOpt[nVar+optAux], vecMonOpt[nVar+k], 1. );
    unsigned i1 = 0, i2 = 0;
    for( auto i1mon = _MIP_auxdec1[k].cbegin(); i1mon != _MIP_auxdec1[k].cend(); ++i1mon, ++i1 )
      if( i1mon->get(GRB_DoubleAttr_X) > 0.9 ) break;
    for( auto i2mon = _MIP_auxdec2[k].cbegin(); i2mon != _MIP_auxdec2[k].cend(); ++i2mon, ++i2 )
      if( i2mon->get(GRB_DoubleAttr_X) > 0.9 ) break;
    _insert( matopt, vecMonOpt[i1], vecMonOpt[i2], -1. );
  }

  _SetMon.swap( SetMonOpt );
  _MatFct.swap( MatFctOpt );
  _MatRed.swap( MatRedOpt );
}

template <typename KEY, typename COMP>
inline void
SQuad<KEY,COMP>::_MIP_solve
()
{
  _MIP_options();
  _MIPexcpt = false;
#if defined(MC__USE_GUROBI)
  try{
    _GRBmodel->update();
    //_MIP_display();
    if( options.MIPOUTPUTFILE != "" )
      _GRBmodel->write( options.MIPOUTPUTFILE );
    fedisableexcept(FE_ALL_EXCEPT);
    _GRBmodel->optimize();
    if( options.MIPDISPLEVEL )
      std::cout << "  #auxiliary variables: " << _GRBmodel->get( GRB_DoubleAttr_ObjVal ) << std::endl;
  }
  catch(GRBException& e){
    if( options.MIPDISPLEVEL )
      std::cout << "Error code = " << e.getErrorCode() << std::endl
                << e.getMessage() << std::endl;
    _MIPexcpt = true;
  }
#endif
  _MIP_display();
}

template <typename KEY, typename COMP>
inline void
SQuad<KEY,COMP>::_MIP_display
( std::ostream& os )
const
{
  bool first;
  os << std::fixed << std::setprecision(0);
  os << "z  = ";
  for( auto const& grbvar : _MIP_auxbin )
    os << std::round( grbvar.get(GRB_DoubleAttr_X) ) << " ";
  os << std::endl;
  os << "y1 = ";
  first = true;
  for( auto const& grbvec : _MIP_mondec1 ){
    if( !first ) os << "     ";
    for( auto const& grbvar : grbvec )
      os << std::round( grbvar.get(GRB_DoubleAttr_X) ) << " ";
    os << std::endl;
    first = false;
  }
  os << "y2 = ";
  first = true;
  for( auto const& grbvec : _MIP_mondec2 ){
    if( !first ) os << "     ";
    for( auto const& grbvar : grbvec )
      os << std::round( grbvar.get(GRB_DoubleAttr_X) ) << " ";
    os << std::endl;
    first = false;
  }
  os << "w1 = ";
  first = true;
  for( auto const& grbvec : _MIP_auxdec1 ){
    if( !first ) os << "     ";
    for( auto const& grbvar : grbvec )
      os << std::round( grbvar.get(GRB_DoubleAttr_X) ) << " ";
    os << std::endl;
    first = false;
  }
  os << "w2 = ";
  first = true;
  for( auto const& grbvec : _MIP_auxdec2 ){
    if( !first ) os << "     ";
    for( auto const& grbvar : grbvec )
      os << std::round( grbvar.get(GRB_DoubleAttr_X) ) << " ";
    os << std::endl;
    first = false;
  }
  os << "b  = ";
  first = true;
  for( auto const& grbvec : _MIP_auxexp ){
    if( !first ) os << "     ";
    for( auto const& grbvar : grbvec )
      os << std::round( grbvar.get(GRB_DoubleAttr_X) ) << " ";
    os << std::endl;
    first = false;
  }

}

template <typename KEY, typename COMP>
inline void
SQuad<KEY,COMP>::_MIP_options
()
{
#if defined(MC__USE_GUROBI)
  // Gurobi options
  _GRBmodel->getEnv().set( GRB_IntParam_OutputFlag,        options.MIPDISPLEVEL );
  _GRBmodel->getEnv().set( GRB_IntParam_Method,            options.LPALGO );
  _GRBmodel->getEnv().set( GRB_IntParam_Presolve,          options.LPPRESOLVE );
  _GRBmodel->getEnv().set( GRB_IntParam_Threads,           options.MIPTHREADS );
  _GRBmodel->getEnv().set( GRB_IntParam_ConcurrentMIP,     options.MIPCONCURRENT );
  _GRBmodel->getEnv().set( GRB_IntParam_MIPFocus,          options.MIPFOCUS );
  _GRBmodel->getEnv().set( GRB_DoubleParam_Heuristics,     options.MIPHEURISTICS );
  _GRBmodel->getEnv().set( GRB_DoubleParam_FeasibilityTol, options.LPFEASTOL );
  _GRBmodel->getEnv().set( GRB_DoubleParam_OptimalityTol,  options.LPOPTIMTOL );
  _GRBmodel->getEnv().set( GRB_DoubleParam_MIPGap,         options.MIPRELGAP );
  _GRBmodel->getEnv().set( GRB_DoubleParam_MIPGapAbs,      options.MIPABSGAP );
  _GRBmodel->getEnv().set( GRB_DoubleParam_TimeLimit,      options.MIPTIMELIMIT );
#endif
}

template <typename KEY, typename COMP>
inline void
SQuad<KEY,COMP>::_MIP_reset
()
{
  _MIP_auxbin.clear();
  _MIP_auxexp.clear();
  _MIP_auxdec1.clear();
  _MIP_auxdec2.clear();
  _MIP_mondec1.clear();
  _MIP_mondec2.clear();
#if defined(MC__USE_GUROBI)
  delete _GRBmodel;
  _GRBmodel = new GRBModel( *_GRBenv );
#endif
}
#endif

template <typename KEY, typename COMP>
inline typename SQuad<KEY,COMP>::key_SQuad&
SQuad<KEY,COMP>::_reorder
( key_SQuad&& pMon )
const
{
  if( lt_SMon<COMP>()( *pMon.second, *pMon.first ) )
    std::swap( pMon.first, pMon.second );
  return pMon;
}

template <typename KEY, typename COMP>
inline typename SQuad<KEY,COMP>::set_SMon
SQuad<KEY,COMP>::_prodmon
( t_SMon const& mon1, t_SMon const& mon2 )
const
{
#ifdef MC__SQUAD_DEBUG_PRODMON
  std::cout << "mon1:\n" << mon1.display(options.BASIS) << std::endl;
  std::cout << "mon2:\n" << mon2.display(options.BASIS) << std::endl;
#endif
  set_SMon prodmon;
  switch( options.BASIS ){
   // Monomial basis representation
   case Options::MONOM:
     prodmon.insert( mon1+mon2 );
     break;

   // Chebyshev basis representation
   case Options::CHEB:
    prodmon.insert( t_SMon() );
    for( auto const& [ivar,iord] : (mon1+mon2).expr ){
      set_SMon prodmon2;
      auto&& it1 = mon1.expr.find( ivar );
      auto&& it2 = mon2.expr.find( ivar );
      if( it1 != mon1.expr.end() && it2 != mon2.expr.end() ){
        unsigned const& iord1 = it1->second;
        unsigned const& iord2 = it2->second;
        assert( iord1 && iord2 );
        for( auto const& mon3 : prodmon ){
          prodmon2.insert( mon3 + t_SMon( iord1+iord2, {std::make_pair( ivar, iord1+iord2 )} ) );
          if( iord1 > iord2 )
            prodmon2.insert( mon3 + t_SMon( iord1-iord2, {std::make_pair( ivar, iord1-iord2 )} ) );
          else if( iord1 < iord2 )
            prodmon2.insert( mon3 + t_SMon( iord2-iord1, {std::make_pair( ivar, iord2-iord1 )} ) );
          else
            prodmon2.insert( mon3 );
        }
      }
      else if( it1 != mon1.expr.end() ){
        unsigned const& iord1 = it1->second;
        assert( iord1 );
        for( auto const& mon3 : prodmon )
          prodmon2.insert( mon3 + t_SMon( iord1, {std::make_pair( ivar, iord1 )} ) );
      }    
      else{
        unsigned const& iord2 = it2->second;
        assert( iord2 );
        for( auto const& mon3 : prodmon )
          prodmon2.insert( mon3 + t_SMon( iord2, {std::make_pair( ivar, iord2 )} ) );
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

template <typename KEY, typename COMP>
inline bool
SQuad<KEY,COMP>::_insert
( map_SQuad& mat, t_SMon const* pMon, double const coef, bool const add )
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
}

template <typename KEY, typename COMP>
inline bool
SQuad<KEY,COMP>::_insert
( map_SQuad& mat, t_SMon const* pLMon, t_SMon const* pRMon,
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

template <typename KEY, typename COMP>
inline bool
SQuad<KEY,COMP>::_insert
( map_SQuad& mat, map_SPoly& mapmon, t_SMon const* pLMon,
  t_SMon const* pRMon, double const coef, bool const add )
{
  // Extra lower-order terms generate by Chebyshev product
  auto&& prodmon = _prodmon( *pLMon, *pRMon );
  unsigned const nprod = prodmon.size();
  auto&& itmon = prodmon.crbegin(); 
  for( ++itmon; itmon != prodmon.crend(); ++itmon ){
    auto [itcmon,ins] = mapmon.insert( std::make_pair( *itmon, -coef ) );
    if( !ins ) itcmon->second -= coef; 
    if( itcmon->second == 0. ) mapmon.erase( itcmon );
  }
  
  // New entry in quadratic form
  return _insert( mat, pLMon, pRMon, coef * nprod, add );
}

template <typename KEY, typename COMP>
inline void
SQuad<KEY,COMP>::_candidates
( set_SQuad& CandidateDec, t_SMon const& mon )
const
{
  for( auto&& mon2 : _SetMon ){
    if( !mon2.subseteq( mon ) ) continue;
    auto&& itmon3 = _SetMon.find( mon - mon2 );
    if( itmon3 == _SetMon.end() || ( &mon2 != &*itmon3 && lt_SMon<COMP>()( mon2, *itmon3 ) ) ) continue;
#ifdef MC__SQUAD_DEBUG_DECOMP
    std::cout << "Candidate: " << mon.display(options.BASIS)
              << " = " << mon2.display(options.BASIS)
              << " · " << itmon3->display(options.BASIS)
              << std::endl;
#endif
    CandidateDec.insert( std::make_pair( &mon2, &(*itmon3) ) );
  }
}

template <typename KEY, typename COMP>
inline typename SQuad<KEY,COMP>::key_SQuad
SQuad<KEY,COMP>::_decompose
( t_SMon const& mon )
{
  // Possible decompositions using existing monomial in _SetMon
  set_SQuad CandidateDec;
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
    set_pSMon CandidateMon;
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
      t_SMon const* pmon2 = *CandidateMon.rbegin();
      t_SMon mon3( mon - *pmon2 );
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
    t_SMon mon2;
    for( auto const& [ivar,iord] : mon.expr ){
#ifdef MC__SQUAD_CHECK
      assert( iord == 1 );
#endif
      mon2 += t_SMon( ivar );
      if( ++count >= mon.tord / 2 + mon.tord % 2 ) break;
    }
    t_SMon mon3( mon - mon2 );
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
    t_SMon mon2;
    for( auto&& [ivar,iord] : mon.expr )
      if( iord % 2 )
        mon2 += t_SMon( ivar );
#ifdef MC__SQUAD_DEBUG_DECOMP
    std::cout << "Decomposed: " << mon.display(options.BASIS)
              << " = " << mon2.display(options.BASIS)
              << " · " << (mon-mon2).display(options.BASIS)
              << std::endl;
#endif
    t_SMon mon3( mon - mon2 );
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
  t_SMon mon2 = mon / 2;
  // Decompose mon2
  auto&& itmon2 = _SetMon.find( mon2 );
  if( itmon2 == _SetMon.end() ) itmon2 = _subexpression( mon2 );
  //return _reorder( std::make_pair( &(*itmon2), &(*itmon2) ) );   
  return std::make_pair( &(*itmon2), &(*itmon2) );
}

template <typename KEY, typename COMP>
inline typename SQuad<KEY,COMP>::set_SMon::iterator
SQuad<KEY,COMP>::_subexpression
( t_SMon const& mon )
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
  _MatRed.push_back( map_SQuad() );
  auto& mat = _MatRed.back();
#ifdef MC__SQUAD_DEBUG_DECOMP
  std::cerr << "SQuad::_subexpression, &mat = " << &mat << std::endl;
#endif
  auto itmon = _SetMon.insert( mon ).first;
  map_SPoly mapmon;
  bool ins = _insert( mat, &(*itmon), -1. )
          && _insert( mat, mapmon, plmon, prmon, 1. );
#ifdef MC__SQUAD_CHECK
    assert( ins );
#endif
  //for( auto it=mapmon.crbegin(); it!=mapmon.crend(); ++it ){
    //auto const& [monlow,coeflow] = *it;
  for( ; !mapmon.empty(); ){
    // Local copy of next monomial, then erase
    auto const [monlow,coeflow] = options.ORDER==Options::INC? *mapmon.cbegin(): *mapmon.crbegin();
    mapmon.erase( monlow );

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

template <typename KEY, typename COMP>
inline void
SQuad<KEY,COMP>::_reduction
( t_SMon const& mon )
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
    t_SMon montot( mon + *itmon2 );
    auto itmon3 = itmon2;
    for( ++itmon3; itmon3 != _SetMon.end(); ++itmon3 ){
      if( !itmon3->subseteq( montot ) ) continue;
      auto itmon4 = _SetMon.find( montot - *itmon3 );

      // Prevent duplication
      switch( options.ORDER ){
       case Options::INC:      
        if( itmon4 == _SetMon.end() || lt_SMon<COMP>()( *itmon4, *itmon3 ) || itmon == itmon4 ) continue;
        break;
       case Options::DEC:
        if( itmon4 == _SetMon.end() || lt_SMon<COMP>()( *itmon3, *itmon4 ) || itmon == itmon3 ) continue;
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
      _MatRed.push_back( map_SQuad() );
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

template <typename KEY, typename COMP>
inline void
SQuad<KEY,COMP>::_reset
()
{
  _SetMon.clear();
  _MatFct.clear();
  _MatRed.clear();
  _MatPSD.clear();
}

} // namespace mc

#endif
