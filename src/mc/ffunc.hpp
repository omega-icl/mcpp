// Copyright (C) 2013-2022 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_FFUNC Construction, Manipulation and Evaluation of Factorable Functions
\author Benoit Chachuat & OMEGA Research Group (http://www3.imperial.ac.uk/environmentenergyoptimisation)
\date 2022
\bug No known bugs.

Originally introduced by McCormick [McCormick, 1976] for the development of a convex/concave relaxation arithmetic, <B>factorable functions</B> cover an extremely inclusive class of functions which can be represented finitely on a computer by means of a code list or a computational graph involving atom operations. These are typically unary and binary operations within a library of atom operators, which can be based for example on the C-code library <tt>math.h</tt>. Besides convex/concave relaxations, factorable functions find applications in automatic differentiation (AD) [Naumann, 2009] as well as in interval analysis [Moore <I>et al.</I>, 2009] and other set arithmetics [Chachuat <I>et al.</i>, 2015].

Factorable functions can be represented using <b>directed acyclic graphs (DAGs)</b>, whose nodes are subexpressions and whose directed edges are computational flows [Schichl & Neumaier, 2005]. A key advantage of DAGs is that they can handle common subexpressions shared by several functions effectively.

The classes mc::FFGraph, mc::FFBase, mc::FFVar and mc::FFOp defined in <tt>ffunc.hpp</tt> provide an implementation of DAGs for factorable functions in MC++. They enable their evaluation using any of the arithmetics supported by MC++, including interval ranges (mc::Interval), McCormick relaxations (mc::McCormick), Taylor and Chebyshev models (mc::TVar, mc::CVar, mc::SCVar), spectral bounds (mc::Specbnd), polyhedral relaxations (mc::PolVar), and interval superposition models (mc::ISVar). They also enable their symbolic manipulation, including automatic differentiation, sparse factorization, and quadratization; see:
- \subpage page_SQUAD
- \subpage page_SLIFT
- \subpage page_SELIM
- \subpage page_RLTRED
.

\section sec_FFUNC_dag How do I construct the DAG of a factorable function?

For illustration, suppose we want to construct a DAG for the factorable function \f${\bf f}:\mathbb{R}^4\to\mathbb{R}^2\f$ defined by
\f{align*}
  {\bf f}(x_0,x_1,x_2,x_3) = \left(\begin{array}{c} x_2x_3-x_0\\ x_0(\exp(x_2x_3)+3.0)^4)+x_1\end{array}\right)
\f}

The constructions require the header file <tt>ffunc.hpp</tt> to be included:

\code
      #include "ffunc.hpp"
\endcode

An environment mc::FFGraph is first defined for recording the factorable function DAG. All four variables mc::FFVar participating in that function are then defined in the environment using the method mc::FFVar::set:

\code
      mc::FFGraph DAG;
      const unsigned int NX = 4;
      mc::FFVar X[NX];
      for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
\endcode

The two components of the factorable function can be defined next:

\code
      const unsigned int NF = 2;
      mc::FFVar F[NF]
        = { X[2]*X[3]-X[0],
            X[0]*pow(exp(X[2]*X[3])+3.,4)+X[1] };
      std::cout << DAG;
\endcode

The last line displays the following information about the factorable function DAG:

\verbatim
    DAG VARIABLES:
      X0     => { Z1 Z7 }
      X1     => { Z8 }
      X2     => { Z0 }
      X3     => { Z0 }

    DAG INTERMEDIATES:
      Z0	<=  X2 x X3		 => { Z1 Z2 }
      Z1	<=  Z0 - X0		 => { }
      Z2	<=  EXP( Z0 )		 => { Z4 }
      Z4	<=  Z2 + Z3		 => { Z6 }
      Z6	<=  IPOW( Z4, Z5 )	 => { Z7 }
      Z7	<=  X0 x Z6		 => { Z8 }
      Z8	<=  X1 + Z7		 => { }
      Z5	<=  4(I)		 => { Z6 }
      Z3	<=  3.1(D)		 => { Z4 }
\endverbatim

Observe that 9 auxiliary variables, \f$z_0,\ldots,z_8\f$, have been created in the DAG, which correspond to the various unary and binary operations in the factorable function expression, as well as the (integer or real) participating constants. Observe, in particular, that the common sub-expression \f$x_2x_3\f$ is detected here; that is, the intermediate \f$z_0\f$ is reused to obtain both subsequent auxiliary variables \f$z_1\f$ and \f$z_2\f$.

At this point, the member function <a>mc::FFGraph::subgraph</a> can be used to generate the subgraph of a DAG corresponding to a given subset of dependent variables. This function returns a list of const pointers <a>mc::FFOp*</a> to the operations participating in the subgraph in order of appearance, which can then be displayed using the member function <a>mc::FFGraph::output</a>:

\code
      DAG.output( DAG.subgraph( NF, F ), " F" );
      DAG.output( DAG.subgraph( 1, &F[0] ), " F0" );
\endcode

Here, the first line generates and displays a subgraph of both components of \f${\bf f}\f$, whereas the second line generates and displays a subgraph of the first component \f$f_0\f$ only:

\verbatim
    OPERATIONS IN SUBGRAPH F:
      X2	<=  VARIABLE
      X3	<=  VARIABLE
      Z0	<=  X2 x X3	
      X0	<=  VARIABLE
      Z1	<=  Z0 - X0	
      X1	<=  VARIABLE
      Z2	<=  EXP( Z0 )	
      Z3	<=  3.1(D)	
      Z4	<=  Z2 + Z3	
      Z5	<=  4(I)	
      Z6	<=  IPOW( Z4, Z5 )
      Z7	<=  X0 x Z6	
      Z8	<=  X1 + Z7	

    DEPENDENTS IN SUBGRAPH F:
      0:  Z1
      1:  Z8

    OPERATIONS IN SUBGRAPH F0:
      X2	<=  VARIABLE
      X3	<=  VARIABLE
      Z0	<=  X2 x X3	
      X0	<=  VARIABLE
      Z1	<=  Z0 - X0	

    DEPENDENTS IN SUBGRAPH F0:
      0:  Z1
\endverbatim

The obtained subgraphs can also be depicted using the (open source) graph plotting program <a href="http://www.graphviz.org/">DOT</a>. The dot files <tt>F.dot</tt> and <tt>F1.dot</tt> can be generated for both subgraphs as follows [which requires the header file <tt>fstream.h</tt>]:

\code
      std::ofstream o_F( "F.dot", std::ios_base::out );
      DAG.dot_script( NF, F, o_F );
      o_F.close();

      std::ofstream o_F0( "F0.dot", std::ios_base::out );
      DAG.dot_script( 1, F, o_F0 );
      o_F0.close();
\endcode

The graphs can be visualized, e.g., after generating SVG files using the command line as:

\verbatim
    $ dot -Tpng -O F.dot;  display F.dot.png
    $ dot -Tpng -O F0.dot; display F0.dot.png
\endverbatim

<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html F.png "Figure: Graph for file F.dot"</TD>
<TD>\image html F0.png "Figure: Graph for file F0.dot"</TD>
</TR>
</TABLE></CENTER>


\section sec_FFUNC_FADBAD How do I obtain the DAG of a factorable function's derivatives?

Derivatives of a factorable function in mc::FFGraph can be obtained with the methods mc::FFGraph::FAD and mc::FFGraph::BAD, which implement the forward and reverse mode of automatic differentiation (AD), respectively. It should be noted that mc::FFGraph does <a>not</a> implement these AD methods per se, but uses the classes fadbad::F and fadbad::B as part of <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A>.

In the forward mode of AD, for instance, entries of the Jacobian matrix of the factorable function \f$f\f$ considered in the previous section can be added to the DAG as follows:

\code
      const mc::FFVar* dFdX_FAD = DAG.FAD( NF, F, NX, X );
      std::cout << DAG;
\endcode

The last line displays the following information about the DAG of the factorable function and its Jacobian:

\verbatim
    DAG VARIABLES:
      X0	 => { Z1 Z7 Z17 Z18 }
      X1	 => { Z8 }
      X2	 => { Z0 Z10 }
      X3	 => { Z0 Z9 }

    DAG INTERMEDIATES:
      Z0	<=  X2 * X3		    => { Z1 Z2 }
      Z1	<=  Z0 - X0		    => { }
      Z2	<=  EXP( Z0 )		=> { Z4 Z9 Z10 }
      Z4	<=  Z2 + Z3		    => { Z6 Z12 }
      Z6	<=  POW( Z4, Z5 )	=> { Z7 }
      Z7	<=  X0 * Z6		    => { Z8 }
      Z8	<=  X1 + Z7		    => { }
      Z9	<=  X3 * Z2		    => { Z15 }
      Z10	<=  X2 * Z2		    => { Z16 }
      Z12	<=  POW( Z4, Z11 )	=> { Z14 }
      Z14	<=  Z12 * Z13		=> { Z15 Z16 }
      Z15	<=  Z9 * Z14		=> { Z17 }
      Z16	<=  Z10 * Z14		=> { Z18 }
      Z17	<=  X0 * Z15		=> { }
      Z18	<=  X0 * Z16		=> { }
      Z11	<=  3(I)		    => { Z12 }
      Z5	<=  4(I)		    => { Z6 }
      Z19	<=  -1(D)		    => { }
      Z21	<=  0(D)		    => { }
      Z20	<=  1(D)		    => { }
      Z3	<=  3(D)		    => { Z4 }
      Z13	<=  4(D)		    => { Z14 }
\endverbatim

Observe that 13 extra auxiliary variables, \f$z_9,\ldots,z_{21}\f$, have been created in the DAG after the application of forward AD. The function mc:FFGraph::FAD returns a const array, whose entries correspond to \f$\frac{\partial f_i}{\partial x_j}\f$ of the Jacobian matrix of \f$f\f$ in the DAG, ordered column-wise as \f$\frac{\partial f_1}{\partial x_1},\ldots,\frac{\partial f_1}{\partial x_n},\frac{\partial f_2}{\partial x_1},\ldots,\frac{\partial f_2}{\partial x_n},\ldots\f$. To prevent memory leaks, this const array should be deleted before becoming out of scope.

As previously, subgraphs can be constructed for all or part of the derivatives, and dot file can be generated for these subgraphs too, e.g.:

\code
      std::ofstream o_dFdX_FAD( "dFdX_FAD.dot", std::ios_base::out );
      DAG.dot_script( NX*NF, dFdX_FAD, o_dFdX_FAD );
      o_dFdX_FAD.close();
      
      std::ofstream o_dF1dX3_FAD( "dF1dX3_FAD.dot", std::ios_base::out );
      DAG.dot_script( 1, &dFdX_FAD[NX+3], o_dF1dX3_FAD );
      o_dF1dX3_FAD.close();

      delete[] dFdX_FAD;
\endcode

The first subgraph created above corresponds to both components of the factorable function \f$f\f$ as well as all eight component of its Jacobian matrix \f$\frac{\partial {\bf f}}{\partial {\bf x}}\f$; the second subgraph is for the Jacobian element  \f$\frac{\partial f_1}{\partial x_3}\f$. The corresponding graphs are shown below.

<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html dFdX_FAD.png "Figure: Graph for file dFdX_FAD.dot"</TD>
<TD>\image html dF1dX3_FAD.png "Figure: Graph for file dF1dX3_FAD.dot"</TD>
</TR>
</TABLE></CENTER>

The backward (or adjoint) method of AD can be applied in a likewise manner using the method mc::FFGraph::BAD instead of mc::FFGRAPH::FAD, everything else being the same.

\code
      const mc::FFVar* dFdX_BAD = DAG.BAD( NF, F, NX, X );

      std::ofstream o_dFdX_BAD( "dFdX_BAD.dot", std::ios_base::out );
      DAG.dot_script( NX*NF, dFdX_BAD, o_dFdX_BAD );
      o_dFdX_BAD.close();
      
      std::ofstream o_dF1dX3_BAD( "dF1dX3_BAD.dot", std::ios_base::out );
      DAG.dot_script( 1, &dFdX_BAD[NX+3], o_dF1dX3_BAD );
      o_dF1dX3_BAD.close();

      delete[] dFdX_BAD;
\endcode


The corresponding graphs are shown below. Note that the reverse mode leads to a DAG of the Jacobian matrix with 10 operations (+,*,pow,exp) only, whereas the forward mode needs 12 operations.

<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html dFdX_BAD.png "Figure: Graph for file dFdX_BAD.dot"</TD>
<TD>\image html dF1dX3_BAD.png "Figure: Graph for file dF1dX3_BAD.dot"</TD>
</TR>
</TABLE></CENTER>

The class mc::FFGraph also supports sparse derivatives, both in forward and backward modes, as well as directional derivatives in forward mode.


\section sec_FFUNC_TAD How do I obtain the DAG of the Taylor expansion of ODE solutions?

Consider a dynamic system of the form
\f{align*}
  \dot{\bf x}(t,{\bf p}) \;=\; {\bf f}({\bf x}(t,{\bf p}),{\bf p})\,,
\f}
where \f${\bf x} \in \mathbb{R}^{n_x}\f$ denotes the state variables, and \f${\bf p} \in \mathbb{R}^{n_p}\f$ a set of (time-invariant) parameters. Assuming that the right-hand side function \f${\bf f}\f$ is sufficiently often continuously differentiable on \f$\mathbb{R}^{n_x}\times\mathbb{R}^{n_p}\f$, a \f$q\f$th-order Taylor expansion in time of the ODE solutions \f$x(\cdot,{\bf p})\f$ at a given time \f$t\f$ reads:
\f{align*}
{\bf x}(t+h,{\bf p}) \;=\; \sum_{i=0}^q h^i {\boldsymbol\phi}_i({\bf x}(t,{\bf p})) + \cal{O}(h^{s+1}) \,,
\f}
where \f${\boldsymbol\phi}_0,\ldots,{\boldsymbol\phi}_q\f$ denote the Taylor coefficients of the solution, defined recursively as
\f{align*}
{\boldsymbol\phi}_0({\bf x},{\bf p}) \;:=\; {\bf x} \qquad \text{and} \qquad
{\boldsymbol\phi}_i({\bf x},{\bf p}) \;:=\; \frac{1}{i} \frac{\partial{\boldsymbol\phi}_{i-1}}{\partial {\bf x}}({\bf x},{\bf p})\, {\bf f}({\bf x},{\bf p}) \quad\text{for $i\geq 1$} \,.
\f}
DAGs for these Taylor coefficients can be generated using the method mc::FFGraph::TAD, which relies upon the class fadbad::T of <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A>.

As a simple illustrative example, consider the scalar linear ODE \f$\dot{x}(t) = x(t)\f$, whose solutions are given by \f$x(t+h) = \exp(h)x(t)\f$. Accordingly, the desired Taylor coefficients are:
\f{align*}
\phi_i({\bf x}) \;:=\; \frac{1}{i!}x \quad\text{for all $i\geq 0$} \,.
\f}
A DAG of these Taylor coefficients can be generated by mc::FFGraph as follows:

\code
  mc::FFGraph DAG;
  mc::FFVar T( &DAG );
  mc::FFVar X( &DAG ); 
  mc::FFVar F = X;

  const unsigned NTE = 10;
  const mc::FFVar* F_TAD = DAG.TAD( NTE, 1, &F, 1, &X, &T );
  std::cout << DAG;

  DAG.output( DAG.subgraph( NTE+1, F_TAD ), " F_TAD" );
  std::ofstream o_TAD( "F_TAD.dot", std::ios_base::out );
  DAG.dot_script( NTE+1, F_TAD, o_TAD );
  o_TAD.close();

  delete[] F_TAD;
\endcode

The resulting DAG of the Taylor coefficients up to order 10 is shown below.

<CENTER><TABLE BORDER=0>
<TR>
<TD>
\verbatim
DAG VARIABLES:
  X0	 => { }
  X1	 => { Z1 }

DAG INTERMEDIATES:
  Z1	<=  X1 * Z0		 => { Z3 }
  Z3	<=  Z1 * Z2		 => { Z5 }
  Z5	<=  Z3 * Z4		 => { Z7 }
  Z7	<=  Z5 * Z6		 => { Z9 }
  Z9	<=  Z7 * Z8		 => { Z11 }
  Z11	<=  Z9 * Z10	 => { Z13 }
  Z13	<=  Z11 * Z12	 => { Z15 }
  Z15	<=  Z13 * Z14	 => { Z17 }
  Z17	<=  Z15 * Z16	 => { }
  Z16	<=  0.1(D)		 => { Z17 }
  Z14	<=  0.111111(D)	 => { Z15 }
  Z12	<=  0.125(D)	 => { Z13 }
  Z10	<=  0.142857(D)	 => { Z11 }
  Z8	<=  0.166667(D)	 => { Z9 }
  Z6	<=  0.2(D)		 => { Z7 }
  Z4	<=  0.25(D)		 => { Z5 }
  Z2	<=  0.333333(D)	 => { Z3 }
  Z0	<=  0.5(D)		 => { Z1 }

OPERATIONS IN SUBGRAPH F_TAD:
  X1	<=  VARIABLE
  Z0	<=  0.5(D)	
  Z1	<=  X1 * Z0	
  Z2	<=  0.333333(D)	
  Z3	<=  Z1 * Z2	
  Z4	<=  0.25(D)	
  Z5	<=  Z3 * Z4	
  Z6	<=  0.2(D)	
  Z7	<=  Z5 * Z6	
  Z8	<=  0.166667(D)	
  Z9	<=  Z7 * Z8	
  Z10	<=  0.142857(D)	
  Z11	<=  Z9 * Z10	
  Z12	<=  0.125(D)	
  Z13	<=  Z11 * Z12	
  Z14	<=  0.111111(D)	
  Z15	<=  Z13 * Z14	
  Z16	<=  0.1(D)	
  Z17	<=  Z15 * Z16	

DEPENDENTS IN SUBGRAPH F_TAD:
  0:  X1
  1:  X1
  2:  Z1
  3:  Z3
  4:  Z5
  5:  Z7
  6:  Z9
  7:  Z11
  8:  Z13
  9:  Z15
  10:  Z17
\endverbatim
<TD>\image html F_TAD.png "Figure: Graph for file F_TAD.dot"</TD>
</TR>
</TABLE></CENTER>

In turn, the resulting DAG of Taylor coefficients may be differentiated using mc::FFGraph::FAD or mc::FFGraph::BAD, or evaluated in any compatible arithmetic as explained next.


\section sec_FFUNC_eval How do I evaluate the DAG of a factorable function in a given arithmetic?

Having created the DAG of a factorable function or its derivatives, one can evaluate these functions in any arithmetic implemented in MC++ using the method mc::FFGraph::eval.

Coming back to our initial example, suppose that we want to compute interval bounds on the first-order derivatives of the factorable function 
\f{align*}
  {\bf f} = \left(\begin{array}{c} x_2x_3-x_0\\ x_0(\exp(x_2x_3)+3.0)^4)+x_1\end{array}\right)
\f}
in the direction \f$(0,1,1,0)\f$, with \f$x_0\in[0,0.5]\f$, \f$x_1\in[1,2]\f$, \f$x_2\in[-1,-0.8]\f$, and \f$x_3\in[0.5,1]\f$.

For simplicity, the default interval type mc::Interval of MC++ is used here:

\code
      #include "ffunc.hpp"
      #include "interval.hpp"
      typedef mc::Interval I;
\endcode

A DAG of the directional derivatives of \f$f\f$ is constructed first:

\code
      // DAG environment
      mc::FFGraph DAG;

      // Independent variables and derivative direction
      const unsigned int NX = 4;
      mc::FFVar X[NX], D[NX] = { 0., 1., 1., 0. };
      for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );

      // Dependent variables
      const unsigned int NF = 2;
      mc::FFVar F[NF]
        = { X[2]*X[3]-X[0],
            X[0]*pow(exp(X[2]*X[3])+3.,4)+X[1] };

      // DAG of directional derivatives
      const mc::FFVar* dFdXdir = DAG.DFAD( NF, F, NX, X, D );
\endcode

In a second step, the DAG of directional derivatives is evaluated in interval arithmetic as follows:

\code
      // Evaluation in interval arithmetic
      I IX[NX] = { I(0,0.5), I(1,2), I(-1,-0.8), I(0.5,1) }, IdFdXdir[NF];
      DAG.eval( NF, dFdXdir, IdFdXdir, NX, X, IX );

      // Display results
      for( unsigned i=0; i<NF; i++ )
        std::cout << "  dF("<< i << ")dX·D = " << IdFdXdir[i] << std::endl;
\endcode

The DAG evaluation can be carried out in sparse Chebyshev model arithmetic likewise:

\code
      #include "scmodel.hpp"
      typedef mc::SCModel<I,mc::FFVar*,mc::lt_FFVar> SCM;
      typedef mc::SCVar<I,mc::FFVar*,mc::lt_FFVar> SCV;
\endcode

\code
      // Evaluation in 3rd-order Chebyshev model arithmetic
      const unsigned ORD = 3;
      SCM CMenv( ORD );
      SCV CMX[NX], CMdFdXdir[NF];
      for( unsigned i=0; i<NX; i++ ) CMX[i].set( &CMenv, i, IX[i] );
      std::vector<SCV> SCVWK;
      DAG.eval( SCVWK, NF, dFdXdir, CMdFdXdir, NX, X, CMX );

      // Display results
      for( unsigned i=0; i<NF; i++ )
        std::cout << "  dF("<< i << ")dX·D = " << CMdFdXdir[i] << std::endl;
\endcode

Note that for repeated function evaluations, it is recommended to pass a pre-sized working array and a subgraph of the function to mc::FFGraph::eval for efficiency. These evaluations produce the following results:

<h3>Evaluation in Interval Arithmetic</h3>
\verbatim
      dF(0)dX·D = [  5.00000e-01 :  1.00000e+00 ]
      dF(1)dX·D = [  1.00000e+00 :  6.72863e+01 ]
\endverbatim

<h3>Evaluation in Sparse Chebyshev Model Arithmetic</h3>
\verbatim
      dF(0)dX·D = 
         7.5000000e-01   0  1
         2.5000000e-01   1  T1[X3]
         R     =  [ 0.0000000e+00, 0.0000000e+00]
         B     =  [ 5.0000000e-01, 1.0000000e+00]

      dF(1)dX·D = 
         1.7166062e+01   0  1
         1.6166062e+01   1  T1[X0]
         1.7305014e+00   1  T1[X2]
         3.4528311e-01   1  T1[X3]
         1.7305014e+00   2  T1[X0]·T1[X2]
         3.4528311e-01   2  T1[X0]·T1[X3]
         5.6453307e-02   2  T2[X2]
         5.0915026e-01   2  T1[X2]·T1[X3]
        -3.9506502e-01   2  T2[X3]
         5.6453307e-02   3  T1[X0]·T2[X2]
         5.0915026e-01   3  T1[X0]·T1[X2]·T1[X3]
        -3.9506502e-01   3  T1[X0]·T2[X3]
         1.5084111e-03   3  T3[X2]
         3.0422535e-02   3  T2[X2]·T1[X3]
        -4.9916695e-02   3  T1[X2]·T2[X3]
         4.7740483e-02   3  T3[X3]
         R     =  [-1.8750872e-01, 1.8750872e-01]
         B     =  [-5.2770965e+00, 3.9414566e+01]
\endverbatim

Backward propagation is also possible through the DAG, e.g. for contraint propagation. Assuming that interval bounds are known for the directional derivatives, we want to tighten the bounds on the independent variables through reverse DAG propagation:

\code
      // Evaluation in interval arithmetic
      I IX[NX] = { I(0,0.5), I(1,2), I(-1,-0.8), I(0.5,1) };
      I IdFdXdir[NF] = { I(0.6,0.9), I(2.,5.) };
      I IINF = 1e20 * I(-1,1);
      int flag = DAG.reval( NF, dFdXdir, IdFdXdir, NX, X, IX, IINF );
      std::cout << "\nDAG interval evaluation w/ " << flag << " forward/backward passes:\n";

      // Display results
      for( unsigned i=0; i<NX; i++ )
        std::cout << "  X(" << i << ") = " << IX[i] << std::endl;
      for( unsigned i=0; i<NF; i++ )
        std::cout << "  dF("<< i << ")dX·D = " << IdFdXdir[i] << std::endl;
\endcode

This evaluation produces the following results:

<h3>Constraint Propagation in Interval Arithmetic</h3>
\verbatim
      DAG interval evaluation w/ 4 forward/backward passes:
      X(0) = [  9.47264e-03 :  1.03696e-01 ]
      X(1) = [  1.00000e+00 :  2.00000e+00 ]
      X(2) = [ -1.00000e+00 : -8.00000e-01 ]
      X(3) = [  6.00000e-01 :  9.00000e-01 ]
      dF(0)dX·D = [  6.00000e-01 :  9.00000e-01 ]
      dF(1)dX·D = [  2.00000e+00 :  5.00000e+00 ]
\endverbatim

In practice, it is paramount to use reverse propagation of verified types; that is, types that account for round-off errors. Otherwise, the behavior could be unreliable and unpredictable; e.g., a feasible set of constraints might be declared infeasible. 


\section sec_FFUNC_ext How do I add an external operation in the DAG of a factorable function?

For illustration, suppose we want to construct a DAG for the factorable function \f$g:\mathbb{R}^3\to\mathbb{R}\f$ defined by
\f{align*}
  g({\bf x}) = \left\|{\bf x}\right\|_2
\f}
Although \f$\left\|\cdot\right\|_2\f$ is not a default operation in mc::FFOp, we may define this operation externally and pass it as a template argument to the DAG class mc::FFGraph. 

A new external operation mc::FFnorm2 is first derived from the base class of mc::FFOp as follows:

\code
      namespace mc
      {
        class FFnorm2
        : public FFOp
        {
        public:
          // Constructors
          FFnorm2
            ()
            : FFOp( (int)EXTERN )
            {}

          // Functor
          FFVar& operator()
            ( unsigned const nVar, FFVar const* pVar )
            const
            {
              auto dep = FFDep();
              for( unsigned i=0; i<nVar; ++i ) dep += pVar[i].dep();
              dep.update( FFDep::TYPE::N );
              return insert_external_operation( *this, dep, nVar, pVar );
            }

          // Evaluation overloads
          template< typename T > void eval
            ( T& vRes, unsigned const nVar, T const* vVar )
            const
            {
              switch( nVar ){
                case 0: vRes = T( 0. ); break;
                case 1: vRes = vVar[0]; break;
                default: vRes = Op<T>::sqr( vVar[0] );
                         for( unsigned i=1; i<nVar; ++i ) vRes += Op<T>::sqr( vVar[i] );
                         vRes = Op<T>::sqrt( vRes ); break;
              }
            }
          void eval
            ( FFVar& vRes, unsigned const nVar, FFVar const* pVar )
            const
            {
              vRes = operator()( nVar, pVar );
            }

          // Properties
          std::string name
            ()
            const
            { return "NORM2"; }
          //! @brief Return whether or not operation is commutative
          bool commutative
            ()
            const
            { return true; }
        };
      }
\endcode

Notice that the identifier passed to mc::FFOp in the constructor of mc::FFnorm2 must be unique to a given external function, starting from the value mc::FFOp::EXTERN. Insertion of this external operation in the DAG of a function is via the mc::FFnorm2::operator() overload, where the type of dependencies is also defined per mc::FFDep. The mc::FFnorm2::eval overloads can be specialized to particular set arithmetics in MC++. Finally, mc::FFnorm2::name and mc::FFnorm2::commutative are used to set, respectively, the name of the external operation and whether the operation is commutative.

Then, a templated environment mc::FFGraph<FFnorm2> is defined for recording the DAG, and the three participating variables are defined using mc::FFVar::set as previously. 

\code
      mc::FFGraph<FFnorm2> DAG;
      const unsigned int NX = 3;
      mc::FFVar X[NX];
      for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
\endcode

After declaring the external operation, we may use it to define the DAG of our function:

\code
      mc::FFnorm2 norm2;
      mc::FFVar F = norm2( NX, X );
      std::cout << DAG;
\endcode

The last line displays the resulting DAG:

\verbatim
    DAG VARIABLES:
      X0	 => { Z0 }
      X1	 => { Z0 }
      X2	 => { Z0 }

    DAG INTERMEDIATES:
      Z0	<=  NORM2( X0, X1, X2 )		 => { }
\endverbatim

Naturally, like with any other DAG, we may evaluate our function in any arithmetic:

\code
      I IX[NX] = { I(0,0.5), I(1,2), I(-1,-0.8) }, IF;
      DAG.eval( 1, &F, &IF, NX, X, IX );
      std::cout << "  " << F << " = " << IF << std::endl;
\endcode

This evaluation produces the following result:

\verbatim
      Z0 = [  1.28062e+00 :  2.29129e+00 ]
\endverbatim


\section sec_FFUNC_err What errors may I encounter while creating or manipulating the DAG of a factorable function?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, an instance of the class mc::FFBase::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during the creation/manipulation of a DAG, and then make the appropriate changes.  Additional exceptions may be sent by the template argument class in propagating a given arithmetic through the DAG. Should an exception be thrown and not caught by the calling program, the execution will abort.


\section sec_FFUNC_refs References

- Chachuat, B, B. Houska, R. Paulen, N. Peric, J. Rajyaguru, M.E. Villanueva, <A href="http://dx.doi.org/10.1016/j.ifacol.2015.09.097">Set-Theoretic Approaches in Analysis, Estimation and Control of Nonlinear Systems</A>, <I>IFAC-PapersOnLine</I>, <b>48</b>(8):981-995, 2015
- McCormick, G.P., <A href="http://dx.doi.org/10.1007/BF01580665">Computability of global solutions to factorable nonconvex programs: Part I. Convex underestimating problems</A>, <i>Mathematical Programming</i>, <b>10</b>(2):147-175, 1976
- Moore, R.E., Cloud, M.J., Kearfott, R.B., <I><A href="http://books.google.co.uk/books/about/Introduction_to_interval_analysis.html?id=tT7ykKbqfEwC&redir_esc=y">Introduction to Interval Analysis</A></I>, SIAM, 2009
- Naumann, U., <I><A href="http://books.google.co.uk/books/about/The_Art_of_Differentiating_Computer_Prog.html?id=OgQuUR4nLu0C&redir_esc=y">The Art of Differentiating Computer Programs: An Introduction to Algorithmic
Differentiation</A></I>, SIAM, 2009
- Puranik, Y., Sahinidis, N.V. <A href="https://doi.org/10.1007/s10601-016-9267-5">Domain reduction techniques for global NLP and MINLP optimization</a>, <i>Constraints</i>, <b>22</b>(3):338-376, 2017.  
- Rajyaguru, J., Villanueva, M.E., Houska, B., Chachuat, B., <A href="https://doi.org/10.1007/s10898-016-0474-9">Chebyshev model arithmetic for factorable functions</a> <i>Journal of Global Optimization</i>, <b>68</b>(2):413-438, 2017
- Schichl, H., Neumaier, A., <a href="http://dx.doi.org/10.1007/s10898-005-0937-x">Interval Analysis on Directed Acyclic Graphs for Global Optimization</a>, <i>Journal of Global Optimization</i>, <b>33</b>:541-562, 2005
- Wechsung, A., Scott, J.K., Watson, H.A.J., Barton, P.I., <A href="https://doi.org/10.1007/s10898-015-0303-6">Reverse propagation of McCormick relaxations</A>, <i>Journal of Global Optimization</i>, <b>63</b>(1):1-36, 2015 
.
*/

#ifndef MC__FFUNC_HPP
#define MC__FFUNC_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <stdarg.h>
#include <set>
#include <list>
#include <vector>
#include <typeinfo>
#include <utility>
#include <type_traits>
#include <tuple>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <climits>

#include "mcop.hpp"
#include "mcfadbad.hpp"
#include "mcfunc.hpp"
#include "mclapack.hpp"
#include "ffdep.hpp"

//#undef  MC__FFUNC_DEBUG
//#undef  MC__FFUNC_DEBUG_TAD

// For time evaluation
//#undef  MC__FFUNC_CPU_EVAL
#ifdef MC__FFUNC_CPU_EVAL
  #include "mctime.hpp"
#endif

// For block decomposition
#ifdef MC__USE_HSL
extern "C" void mc13d_
  ( const int*, const int*, const int*, const int*, const int*, int*, int*, int*, int* );
extern "C" void mc21a_
  ( const int*, const int*, const int*, const int*, const int*, int*, int*, int* );
extern "C" void mc33ad_
  ( const int*, const int*, const int*, int*, const int*, double*, int*, int*, int*,
    int*, int*, int*, int*, int*, int* );
#endif

namespace mc
{
class FFOp;
class FFBase;

//! @brief Structure defining the numeric field of a factorable program variable
////////////////////////////////////////////////////////////////////////
//! mc::FFNum is a C++ structure defining the numeric field of a
//! variable in a factorable function, which can either a real scalar
//! (double), or an integer scalar (int)
////////////////////////////////////////////////////////////////////////
struct FFNum
{
  //! @brief Enumeration type for numeric variables in factorable function
  enum TYPE{
    INT=0,	//!< Integer value
    REAL	//!< Real value
  };

  //! @brief Variable type
  TYPE t;

  //! @brief Integer/real variable value
  union{
    int n;
    double x;
  }; 

  //! @brief Real value
  const double val() const
    { return( t==REAL? x: n ); }

  //! @brief Constructor for an integer variable
  FFNum( const int i=0 )
    : t(INT), n(i)
    {}

  //! @brief Constructor for a real variable
  FFNum( const double d )
    //: t(REAL), x(d)
    { if( std::floor(d)==d && d>=INT_MIN && d<=INT_MAX){ t = INT; n = d; }
      else{ t = REAL;  x = d; } }

  //! @brief Constructor for an integer scalar
  FFNum& operator=
    ( const int i )
    { t = INT; n = i; return *this; }

  //! @brief Constructor for a real scalar
  FFNum& operator=
    ( const double d )
    { if( std::floor(d)==d && d>=INT_MIN && d<=INT_MAX){ t = INT; n = d; }
      else{ t = REAL;  x = d; }
      //t = REAL;  x = d;
      return *this; }

  //! @brief Copy constructor
  FFNum& operator=
    ( const FFNum&num )
    { t = num.t; t==REAL? x=num.x: n=num.n; return *this; }
};

//! @brief Structure comparing values of scalars in factorable functions for equality
////////////////////////////////////////////////////////////////////////
//! mc::eq_FFNum is a C++ structure comparing the numeric field of a
//! FFVar object in a factorable function for equality.
////////////////////////////////////////////////////////////////////////
struct eq_FFNum
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FFNum*Num1, const FFNum*Num2 ) const
    {
      if( Num1->t != Num2->t ) return false;
      switch( Num1->t ){
        case FFNum::INT:  return Num1->n==Num2->n? true: false;
        case FFNum::REAL: return isequal( Num1->x, Num2->x );
      }
    }
};

//! @brief Structure comparing values of scalars in factorable functions for strict inequality
////////////////////////////////////////////////////////////////////////
//! mc::lt_FFNum is a C++ structure comparing the numeric field of a
//! FFVar object in a factorable function for strict inequality.
////////////////////////////////////////////////////////////////////////
struct lt_FFNum
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FFNum*Num1, const FFNum*Num2 ) const
    {
      if( Num1->t < Num2->t ) return true;
      if( Num1->t > Num2->t ) return false;
      switch( Num1->t ){
        case FFNum::INT:  return Num1->n<Num2->n? true: false;
        case FFNum::REAL: return !isequal( Num1->x, Num2->x ) && Num1->x<Num2->x? true: false;
      }
      return false;
    }
};

//! @brief Class defining variables in a factorable function
////////////////////////////////////////////////////////////////////////
//! mc::FFVar is a C++ class defining variables in the factored form of
//! a factorable function.
////////////////////////////////////////////////////////////////////////
class FFVar
////////////////////////////////////////////////////////////////////////
{
  // friends of this class with other classes/structures/operators
  friend class FFBase;
  friend struct lt_FFVar;
  friend std::ostream& operator<< ( std::ostream&, const FFBase& );
  friend std::ostream& operator<< ( std::ostream&, const FFOp& );

  // friends of this class for operator and function overloading
  friend std::ostream& operator<< ( std::ostream&, const FFVar& );
  friend FFVar operator+ ( const FFVar& );
  friend FFVar operator+ ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator+ ( const V&, const FFVar& );
  template <typename V> friend FFVar operator+ ( const FFVar&, const V& );
  friend FFVar operator- ( const FFVar& );
  friend FFVar operator- ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator- ( const V&, const FFVar& );
  template <typename V> friend FFVar operator- ( const FFVar&, const V& );
  friend FFVar operator* ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator* ( const V&, const FFVar& );
  template <typename V> friend FFVar operator* ( const FFVar&, const V& );
  friend FFVar operator/ ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator/ ( const V&, const FFVar& );
  template <typename V> friend FFVar operator/ ( const FFVar&, const V& );
  friend FFVar sum   ( const unsigned int, const FFVar* );
  friend FFVar prod  ( const unsigned int, const FFVar* );
  friend FFVar monom ( const unsigned int, const FFVar*, const unsigned*, const bool );
  friend FFVar max ( const FFVar&, const FFVar& );
  friend FFVar max ( const unsigned int, const FFVar* );
  friend FFVar max ( const double&, const FFVar& );
  friend FFVar max ( const FFVar&, const double& );
  friend FFVar min ( const FFVar&, const FFVar& );
  friend FFVar min ( const unsigned int, const FFVar* );
  friend FFVar min ( const double&, const FFVar& );
  friend FFVar min ( const FFVar&, const double& );
  friend FFVar inter ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar inter ( const V&, const FFVar& );
  template <typename V> friend FFVar inter ( const FFVar&, const V& );
  friend FFVar inv   ( const FFVar& );
  friend FFVar sqr   ( const FFVar& );
  friend FFVar exp   ( const FFVar& );
  friend FFVar log   ( const FFVar& );
  friend FFVar xlog  ( const FFVar& );
  friend FFVar sqrt  ( const FFVar& );
  friend FFVar cos   ( const FFVar& );
  friend FFVar sin   ( const FFVar& );
  friend FFVar tan   ( const FFVar& );
  friend FFVar acos  ( const FFVar& );
  friend FFVar asin  ( const FFVar& );
  friend FFVar atan  ( const FFVar& );
  friend FFVar cosh  ( const FFVar& );
  friend FFVar sinh  ( const FFVar& );
  friend FFVar tanh  ( const FFVar& );
  friend FFVar erf   ( const FFVar& );
  friend FFVar erfc  ( const FFVar& );
  friend FFVar fabs  ( const FFVar& );
  friend FFVar fstep ( const FFVar& );
  friend FFVar bstep ( const FFVar& );
  friend FFVar pow   ( const FFVar&, const int );
  friend FFVar pow   ( const FFVar&, const double );
  friend FFVar pow   ( const FFVar&, const FFVar& );
  friend FFVar pow   ( const double, const FFVar& );
  friend FFVar cheb  ( const FFVar&, const unsigned );

public:

  // other operator overloadings
  bool operator== ( const FFVar& ) const;
  bool operator!= ( const FFVar& ) const;
  FFVar& operator= ( const FFVar& );
  FFVar& operator= ( const int );
  FFVar& operator= ( const double );
  template <typename V> FFVar& operator+= ( const V& );
  template <typename V> FFVar& operator-= ( const V& );
  template <typename V> FFVar& operator*= ( const V& );
  template <typename V> FFVar& operator/= ( const V& );

  typedef std::list< FFOp* > t_Ops;
  typedef typename std::pair< FFOp*, t_Ops > pt_Ops;

  /** @defgroup FFunc Construction, Manipulation and Evaluation of DAGs for Factorable Functions
   *  @{
   */
  //! @brief Index for unreferenced variables in DAG
  static const long NOREF = -33;
  //! @brief Default name for independent variables in DAG
  static const std::string VARNAME;
  //! @brief Default name for auxiliary variables in DAG
  static const std::string AUXNAME;
  //! @brief Enumeration type for variables in DAG
  enum TYPE{
    VAR=0,	//!< Original variable
    AUX,	//!< Auxiliary variable
    CINT,	//!< Integer constant
    CREAL	//!< Real constant
  };
  //! @brief Typedef for variable identifier in DAG
  typedef std::pair< TYPE, long > pt_idVar;
  /** @} */

private:
  //! @brief Pointer to underlying factorable function DAG - _dag := NULL for variable identifier NOREF
  mutable FFBase *_dag;
  //! @brief Identifier (type and index)
  pt_idVar _id;
  //! @brief Numeric field (integer or real)
  mutable FFNum _num;
  //! @brief Dependence and linearity
  FFDep _dep;
  //! @brief Pointer to value field - has to be const_cast'ed in order to retreive original pointer type
  mutable void *_val;
  //! @brief Movability attribute
  mutable bool _mov;
  //! @brief Constness attribute
  mutable bool _cst;
  //! @brief Pointer to parent (_ops.first) and children (_ops.second) operations - _ops.first := NULL for unreferenced constants
  pt_Ops _ops;
  //! @brief Non-default name
  std::string _nam;

public:

  /** @ingroup FFunc
   *  @{
   */
  //! @brief Constructor for variable in DAG <a>*dag</a>
  FFVar
    ( FFBase* dag, std::string const& name="" );

  //! @brief Constructor for variable in DAG <a>*dag</a> with constant real value <a>d</a>
  FFVar
    ( FFBase* dag, double const d );

  //! @brief Constructor for variable in DAG <a>*dag</a> with constant integer value <a>i</a>
  FFVar
    ( FFBase* dag, int const i );

  //! @brief Attach variable to DAG <a>*dag</a>.
  FFVar& set
    ( FFBase* dag, std::string const& name="" )
    { *this = FFVar( dag, name );
      return *this; }

  //! @brief Attach variable to DAG <a>*dag</a>.
  FFVar& set
    ( FFBase* dag, double const& d )
    { *this = FFVar( dag, d );
      return *this; }

  //! @brief Attach variable to DAG <a>*dag</a>.
  FFVar& set
    ( FFBase* dag, int const i )
    { *this = FFVar( dag, i );
      return *this; }

  //! @brief Constructor for integer constant
  FFVar
    ( int const i=0 )
    : _dag( nullptr ), _id( CINT, NOREF ), _num( i ), _dep( i ), _val( nullptr ),
      _mov( false ), _cst( true ), _nam( "" )
    { _ops.first = nullptr; }

  //! @brief Constructor for real parameter
  FFVar
    ( double const& d )
    : _dag( nullptr ), _id( CREAL, NOREF ), _num( d ), _dep( d ), _val( nullptr ),
      _mov( false ), _cst( true ), _nam( "" )
    { _ops.first = nullptr;
      if( _num.t == FFNum::INT ) _id.first = CINT; }

  //! @brief Copy constructor
  FFVar
    ( FFVar const& Var )
    : _dag( Var._dag ), _id( Var._id ), _num( Var._num ), _dep( Var._dep ),
      _val( Var._val ), _mov( Var._mov ), _cst( Var._cst ), _ops( Var._ops ),
      _nam( Var._nam )
    {}
  /** @} */

private:

  //! @brief Constructor for auxiliary variable in factorable function <a>*dag</a> defined from operation <a>*Op</a>
  FFVar
    ( FFBase* dag, FFDep const& dep, FFOp* op=nullptr );

  //! @brief Constructor for a variable with identifier <a>id</a> in factorable function <a>*dag</a>
  FFVar
    ( FFBase* dag, pt_idVar const& id, std::string const& name="" );
    
public:

  /** @ingroup FFunc
   *  @{
   */
  //! @brief Get variable identifier
  std::pair<TYPE,long> const id
    ()
    const
    { return _id; }

  //! @brief Get reference to variable identifier
  std::pair<TYPE,long>& id
    ()
    { return _id; }

  //! @brief Get const reference to variable numeric field
  const FFNum& num
    ()
    const
    { return _num; }

  //! @brief Get/set const reference to variable numeric field
  FFNum& num
    ()
    { return _num; }

  //! @brief Get const reference to variable dependencies
  FFDep const& dep
    ()
    const
    { return _dep; }

  //! @brief Get/set reference to variable dependencies
  FFDep& dep
    ()
    { return _dep; }

  //! @brief Get const pointer to defining operation
  pt_Ops const ops
    ()
    const
    { return _ops; }

  //! @brief Get/set pointer to defining operation
  pt_Ops& ops
    ()
    { return _ops; }

  //! @brief Get const pointer to factorable function dag
  FFBase* dag
    ()
    const
    { return _dag; }

  //! @brief Get/set pointer to factorable function dag
  FFBase*& dag
    ()
    { return _dag; }

  //! @brief Get/set pointer to value field
  void*& val
    ()
    { return _val; }

  //! @brief Get pointer to value field
  template <typename U>
  void reset_val
    ( U const&  U_dum )
    { if( !_val ) return; delete static_cast<U*>( _val ); _val = 0; }

  //! @brief Get/set movability attribute
  bool& mov
    ()
    { return _mov; }

  //! @brief Get variable name
  std::string name
    ( bool const user=false )
    const
    { return ( user || !_nam.empty() ) ? _nam : _name(_id); }

  //! @brief Get constness
  const bool cst
    ()
    const
    { return _cst; }

  //! @brief Get/set constness
  bool& cst
    ()
    { return _cst; }

  //! @brief Set variable at a constant value
  void set
    ( std::string const& name )
    { _nam = name; }

  //! @brief Set variable at a constant value
  void set
    ( const double d )
    const;

  //! @brief Set variable at a constant value
  void set
    ( const int i )
    const;

  //! @brief Unset variable at a constant value
  void unset
    ()
    const;
  /** @} */

private:

  //! @brief Return string with variable name for identifier <a>id</a>
  static std::string _name
    ( std::pair<TYPE,long> const& id )
    {
      std::ostringstream ovar;
      ovar << (id.first==VAR? VARNAME: AUXNAME) << id.second;
      return ovar.str();
    }
};

inline long const FFVar::NOREF;
inline std::string const FFVar::VARNAME = "V";
inline std::string const FFVar::AUXNAME = "Z";

//! @brief Structure comparing variable identifiers in a factorable function for ordering in set FFBase::_Vars
////////////////////////////////////////////////////////////////////////
//! mc::lt_FFVar is a C++ structure comparing variable identifiers in a
//! factorable function for ordering in set FFBase::_Vars.
////////////////////////////////////////////////////////////////////////
struct lt_FFVar
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( FFVar const* Var1, FFVar const* Var2 )
    const
    {
      if( !Var1 ) return false;
      if( !Var2 ) return true;
      // Order variables/constants w.r.t. their types first
      if( Var1->_id.first < Var2->_id.first ) return true;
      if( Var1->_id.first > Var2->_id.first ) return false;
      // If variables, order w.r.t. their index next
      switch( Var1->_id.first ){
        case FFVar::VAR:
        case FFVar::AUX:
          if( Var1->_id.second < Var2->_id.second ) return true;
          if( Var1->_id.second > Var2->_id.second ) return false;
          break;
        case FFVar::CINT: case FFVar::CREAL:
          lt_FFNum ltNum;
          return ltNum( &Var1->_num, &Var2->_num );
          break;
      }
      return false;
    }
};

//! @brief Class defining operations in a factorable function
////////////////////////////////////////////////////////////////////////
//! mc::FFOp is a C++ class defining operations in the factored form of
//! a factorable function.
////////////////////////////////////////////////////////////////////////
class FFOp
////////////////////////////////////////////////////////////////////////
{
public:

  /** @ingroup FFunc
   *  @{
   */
  //! @brief Enumeration type for unary and binary operations
  enum TYPE{
    CNST=0, VAR,
    PLUS, SHIFT, NEG, MINUS, TIMES, SCALE, DIV, INV,
    PROD, IPOW, DPOW, CHEB, SQR, SQRT, EXP, LOG, XLOG,
    SIN, COS, TAN, ASIN, ACOS, ATAN, COSH, SINH, TANH,
    ERF, FABS, FSTEP, MINF, MAXF, INTER, EXTERN
  };

  //! @brief Constructor for unary (or default) operation
  FFOp
    ( int const top=CNST, FFVar* plop=nullptr, FFVar* pres=nullptr );

  //! @brief Constructor for binary operation
  FFOp
    ( int const top, FFVar* plop, FFVar* prop, FFVar* pres );

  //! @brief Constructor for n-ary operation
  FFOp
    ( int const top, unsigned const nops, FFVar** pops, FFVar* pres );

  //! @brief Destructor
  virtual ~FFOp
    ()
    {}

  //! @brief Type of operation
  int type;
  //! @brief Pointer to operation result
  FFVar* pres;
  //! @brief Vector of operands
  std::vector<FFVar*> pops;
  //! @brief Flag for current operation (during a DAG traversal)
  mutable int iflag;


  //! @brief Set unary operand and result
  FFOp& set
    ( FFVar* lop, FFVar* res );
  //! @brief Set binary operands and result
  FFOp& set
    ( FFVar* lop, FFVar* rop, FFVar* res );
  //! @brief Set n-ary operands and result
  FFOp& set
    ( unsigned const nop, FFVar** ops, FFVar* res );

  //! @brief Propagate subset of operations participating in subgraph
  void propagate_subgraph
    ( std::list< std::pair< FFOp const*, bool > >& Ops )
    const;
  //! @brief Reset mc::FFVar::_val field in subgraph
  template <typename U>
  void reset_val_subgraph
    ( U const&  U_dum )
    const;
  //! @brief Propagate script for DAG using DOT and display to <a>os</a>
  void generate_dot_script
    ( std::ostream&os )
    const;
  //! @brief Append script for current operation using DOT to <a>os</a>
  void append_dot_script
    ( std::ostream&os )
    const;
  //! @brief Append script for factor <a>fname</a> using DOT to <a>os</a>
  void append_dot_script_factor
    ( std::ostream&os, const unsigned int fontsize )
    const;
  //! @brief Append script for variable/contant using DOT to <a>os</a>
  void append_dot_script_variable
    ( std::ostream&os, const unsigned int fontsize )
    const;

  //! @brief Evaluate operation in U arithmetic, putting the result at position <a>itU</a>
  template <typename U>
  void evaluate
    ( typename std::vector<U>::iterator itU, U const* pU_dum, bool const movU )
    const;
  //! @brief Forward operation propagation in U arithmetic
  template <typename U>
  bool tighten_forward
    ( const U* pU_dum )
    const;
  //! @brief Backward operation propagation in U arithmetic
  template <typename U>
  bool tighten_backward
    ( const U* pU_dum )
    const;

  //! @brief Evaluate external operation in U arithmetic, putting the result at position <a>itU</a>
  template <typename U, typename ExtOp, typename... NextOps>
  void evaluate_external
    ( typename std::vector<U>::iterator itU, U const* pU_dum, bool const movU,
      ExtOp op, std::tuple<NextOps...> ops )
    const;
  //! @brief Forward external operation propagation in U arithmetic
  template <typename U, typename ExtOp, typename... NextOps>
  bool tighten_forward_external
    ( U const* pU_dum, ExtOp op, std::tuple<NextOps...> ops )
    const;
  //! @brief Backward external operation propagation in U arithmetic
  template <typename U, typename ExtOp, typename... NextOps>
  bool tighten_backward_external
    ( U const* pU_dum, ExtOp op, std::tuple<NextOps...> ops )
    const;

  //! @brief Insert Unary operation <a>Op</a> with operand <a>Var</a> in DAG
  template <typename ExtOp>
  static FFVar& insert_external_operation
    ( ExtOp const& Op, FFDep const& dep, FFVar const& Var );
  //! @brief Insert n-ary operation <a>Op</a> with operands <a>Var1</a> and <a>Var2</a> in DAG
  template <typename ExtOp>
  static FFVar& insert_external_operation
    ( ExtOp const& Op, FFDep const& dep, FFVar const& Var1, FFVar const& Var2 );
  //! @brief Insert n-ary operation <a>Op</a> with operand array <a>pVar</a> of size <a>nVar</a> in DAG
  template <typename ExtOp>
  static FFVar& insert_external_operation
    ( ExtOp const& Op, FFDep const& dep, unsigned const nVar, FFVar const* pVar );

  //! @brief Emulate virtual templated eval function
  template< typename U >
  void eval
    ( U& vRes, unsigned const nVar, U const* vVar )
    const;
  //! @brief Emulate virtual templated reval function
  template< typename U >
  bool reval
    ( U const& vRes, unsigned const nVar, U* vVar )
    const;

  //! @brief Return whether or not operation is commutative
  virtual bool commutative
    ()
    const;      
  //! @brief Return operation name for external operations)
  virtual std::string name
    ()
    const;
  //! @brief Compare type-info
  virtual bool sameid
    ( std::type_info const& id )
    const;
  /** @} */
};

//! @brief C++ structure for comparing operations in a factorable program for ordering in set FFBase::_Ops
////////////////////////////////////////////////////////////////////////
//! mc::lt_FFOp is a C++ structure for comparing operations in a
//! factorable program based on their types and operands for ordering
//! in set FFBase::_Ops.
////////////////////////////////////////////////////////////////////////
struct lt_FFOp
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( FFOp const* Op1, FFOp const* Op2 )
    const
    {
      // Sort by type of operation first
      if( Op1->type < Op2->type ) return true;
      if( Op1->type > Op2->type ) return false;

      // Sort by number of operands next
      if( Op1->pops.size() < Op2->pops.size() ) return true;
      if( Op1->pops.size() > Op2->pops.size() ) return false;

      // Sort by variable type next
      // Different variable ordering is accounted for CAVEAT: This ignores variable ordering; e.g. X1+X2 different from X2+X1
      lt_FFVar ltVar;
      if( Op1->pops.empty() ) return ltVar( Op1->pres, Op2->pres );
      for( auto it1=Op1->pops.begin(), it2=Op2->pops.begin(); 
           it1!=Op1->pops.end() && it2!=Op2->pops.end(); ++it1, ++it2 ){
        if( ltVar( *it1, *it2 ) ) return true;
        if( ltVar( *it2, *it1 ) ) return false;
      }
      return false;
    }
};

//! @brief C++ structure for comparing operations in a factorable program
////////////////////////////////////////////////////////////////////////
//! mc::lt_FFOp is a C++ structure for comparing operations in a
//! factorable program based on their types only.
////////////////////////////////////////////////////////////////////////
struct range_FFOp
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( FFOp const* Op1, FFOp const* Op2 )
    const
    { return ( Op1->type < Op2->type ); }
};


//! @brief C++ structure for holding a subgraph in a DAG
////////////////////////////////////////////////////////////////////////
//! mc::FFSubgraph is a C++ structure for holding a subgraph comparing operations in a
//! factorable program based on their types only.
////////////////////////////////////////////////////////////////////////
struct FFSubgraph
////////////////////////////////////////////////////////////////////////
{
  //! @brief List of (pointers to) operations in DAG
  std::list< std::pair< FFOp const*, bool > > l_op;

  //! @brief Vector of (pointers to) operations defining dependent variables in DAG and whether they may be treated as r-values
  std::vector< FFOp const* > op_dep;

  //! @brief Clear subgraph
  void clear
    ()
    { l_op.clear(); op_dep.clear(); }

  //! @brief Set dependent
  void set_dep
    ( unsigned const& iflag )
    { assert( iflag );
      auto it = l_op.begin();
      std::advance( it, iflag-1 );
      op_dep.push_back( it->first ); }

  //! @brief Set move attributes
  void set_move
    ()
    { for( auto& [op,mov] : l_op ) mov = op->pres->mov(); }
};

//! @brief C++ class representing the DAG of factorable functions
////////////////////////////////////////////////////////////////////////
//! mc::FFBase is a C++ base class representing the DAG of a factorable
//! function.
////////////////////////////////////////////////////////////////////////
class FFBase
////////////////////////////////////////////////////////////////////////
{
  // friends of this class with other classes/structures/operators
  friend class FFVar;
  friend class FFOp;
  friend FFVar operator+ ( const FFVar& );
  friend FFVar operator+ ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator+ ( const V&, const FFVar& );
  template <typename V> friend FFVar operator+ ( const FFVar&, const V& );
  friend FFVar operator- ( const FFVar& );
  friend FFVar operator- ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator- ( const V&, const FFVar& );
  template <typename V> friend FFVar operator- ( const FFVar&, const V& );
  friend FFVar operator* ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator* ( const V&, const FFVar& );
  template <typename V> friend FFVar operator* ( const FFVar&, const V& );
  friend FFVar operator/ ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator/ ( const V&, const FFVar& );
  template <typename V> friend FFVar operator/ ( const FFVar&, const V& );
  friend FFVar sum   ( const unsigned int, const FFVar* );
  friend FFVar prod  ( const unsigned int, const FFVar* );
  friend FFVar monom ( const unsigned int, const FFVar*, const unsigned*, const bool );
  friend FFVar max   ( const FFVar&, const FFVar& );
  friend FFVar max   ( const unsigned int, const FFVar* );
  friend FFVar max   ( const double&, const FFVar& );
  friend FFVar max   ( const FFVar&, const double& );
  friend FFVar min   ( const FFVar&, const FFVar& );
  friend FFVar min   ( const unsigned int, const FFVar* );
  friend FFVar min   ( const double&, const FFVar& );
  friend FFVar min   ( const FFVar&, const double& );
  friend FFVar inter ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar inter ( const V&, const FFVar& );
  template <typename V> friend FFVar inter ( const FFVar&, const V& );
  friend FFVar inv   ( const FFVar& );
  friend FFVar sqr   ( const FFVar& );
  friend FFVar exp   ( const FFVar& );
  friend FFVar log   ( const FFVar& );
  friend FFVar xlog  ( const FFVar& );
  friend FFVar sqrt  ( const FFVar& );
  friend FFVar cos   ( const FFVar& );
  friend FFVar sin   ( const FFVar& );
  friend FFVar tan   ( const FFVar& );
  friend FFVar acos  ( const FFVar& );
  friend FFVar asin  ( const FFVar& );
  friend FFVar atan  ( const FFVar& );
  friend FFVar cosh  ( const FFVar& );
  friend FFVar sinh  ( const FFVar& );
  friend FFVar tanh  ( const FFVar& );
  friend FFVar erf   ( const FFVar& );
  friend FFVar erfc  ( const FFVar& );
  friend FFVar fabs  ( const FFVar& );
  friend FFVar fstep ( const FFVar& );
  friend FFVar bstep ( const FFVar& );
  friend FFVar pow   ( const FFVar&, const int );
  friend FFVar pow   ( const FFVar&, const double );
  friend FFVar pow   ( const FFVar&, const FFVar& );
  friend FFVar pow   ( const double, const FFVar& );
  friend FFVar cheb  ( const FFVar&, const unsigned );

  // friends of this class for operator and function overloading
  friend std::ostream& operator<< ( std::ostream&, const FFBase& );

public:
  /** @ingroup FFunc
   *  @{
   */
  typedef typename FFVar::pt_idVar pt_idVar;
  typedef std::set< FFVar*, lt_FFVar > t_Vars;
  typedef std::set< FFOp*,  lt_FFOp >  t_Ops;
  typedef typename t_Vars::iterator it_Vars;
  typedef typename t_Vars::const_iterator cit_Vars;
  typedef typename t_Ops::iterator  it_Ops;
  typedef typename t_Ops::const_iterator  cit_Ops;
  /** @} */

protected:
  //! @brief Number of original variables in DAG
  unsigned long _nvar;

  //! @brief Number of auxiliary variables in DAG
  unsigned long _naux;

  //! @brief Set of variables in DAG
  t_Vars _Vars;

  //! @brief Set of operations in DAG
  t_Ops _Ops;

  //! @brief Dummy variable used for variable search
  FFVar _dummyVar;

  //! @brief Pointer to current operation in subtree evaluation
  FFOp const* _curOp;

public:
  /** @ingroup FFunc
   *  @{
   */
  //! @brief Default Constructor
  FFBase():
    _nvar( 0 ), _naux( 0 ), _dummyVar( 0 ), _curOp( nullptr )
    {}

  //! @brief Destructor
  virtual ~FFBase
    ()
    { clear(); }

  //! @brief DAG Exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for exception handling
    enum TYPE{
      INIT = 1,		//!< Invalid DAG in variable initialization
      DAG,		    //!< Operation between variables linked to different DAGs
      INTER, 		//!< Empty intersection between constant variables
      MISSVAR,		//!< Missing independent variable during subgraph evaluation
      EVAL,		    //!< Error during subgraph evaluation
      CONSTVAL,		//!< Error due to overriding a constant variable during subgraph evaluation
      INTERN = -1, 	//!< Internal error
      EXTERN = -2, 	//!< Error in external operation
      UNDEF = -33 	//!< Feature not yet implemented
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr=UNDEF ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case INIT:
        return "Invalid DAG in variable initilization";
      case DAG:
        return "Operation between variables linked to different DAGs";
      case INTER:
        return "Empty intersection between constant variables";
      case MISSVAR:
        return "Missing independent variable during subgraph evaluation";
      case EVAL:
        return "Error during subgraph evaluation";
      case CONSTVAL:
        return "Error due to overriding a constant variable during subgraph evaluation";
      case INTERN:
        return "Internal error";
      case EXTERN:
        return "Error in external operation";
      case UNDEF:
        return "Feature not yet implemented";
      default:
        return "Undocumented error";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief DAG Options
  struct Options
  {
    //! @brief Constructor
    Options():
      DETECTSIGNOM( true ),
      CHEBRECURS( true ),
      USEMOVE( true )
      {}
    //! @brief Whether to detect signomial terms as exp(d.log(x)) and handle them as x^d signomial terms
    bool DETECTSIGNOM;
    //! @brief Whether to intersect Chebyshev variables with their recursive expressions -- this may be used to build redundancy in constructing tighter relaxations
    bool CHEBRECURS;
    //! @brief Whether to enable the move semantic during DAG evaluation
    bool USEMOVE;
  } options;

  //! @brief Number of original variables in DAG
  unsigned long nvar
    ()
    const
    { return _nvar; }
  
  //! @brief Number of auxiliary variables in DAG
  unsigned long naux
    ()
    const
    { return _naux; }

  //! @brief Current operation in DAG evaluation
  FFOp const* curOp
    ()
    const
    { return _curOp; }

  //! @brief Current operation in DAG evaluation
  FFOp const*& curOp
    ()
    { return _curOp; }

  //! @brief Reference to set of (all) variables in factorable function
  t_Vars const& Vars
    ()
    const
    { return _Vars; }

  //! @brief Clear DAG (all variables and operations)
  void clear
    ()
    { _clear_variables(); _clear_operations(); _naux = _nvar = 0; }

  //! @brief Looks for the real constant <a>x</a> and adds it if not found
  FFVar const* add_constant
    ( double const x )
    { return _add_constant( x ); }

  //! @brief Extract subgraph corresponding to <a>nDep</a> dependents in array <a>pDep</a>
  FFSubgraph subgraph
    ( unsigned int const nDep, FFVar const* pDep );

  //! @brief Extract subgraph corresponding to dependents indexed by <a>ndxDep</a> in array <a>pDep</a>
  FFSubgraph subgraph
    ( std::set<unsigned> const& ndxDep, FFVar const* pDep );

  //! @brief Extract subgraph corresponding to dependents <a>vDep</a>
  FFSubgraph subgraph
    ( std::vector<const FFVar*> const& vDep );

  //! @brief Extract list of operations corresponding to dependents <a>vDep</a>
  template< typename V, typename COMP>
  FFSubgraph subgraph
    ( std::map<V,FFVar,COMP> const& mDep );

  //! @brief Output list of nodes in <a>Ops</a> to <a>os</a>
  static void output
    ( FFSubgraph const& Ops, std::string const& header="", 
      std::ostream& os=std::cout );

  //! @brief Generate script for DAG visualization of factors <a>*F</a> using DOT
  void dot_script
    ( unsigned int const nDep, FFVar const* pDep, std::ostream& os=std::cout )
    const;

  //! @brief Generate script for DAG visualization of factors <a>*F</a> using DOT
  void dot_script
    ( std::vector<const FFVar*> const& vDep, std::ostream& os=std::cout )
    const;

  //! @brief Compute (symbolic) sum of vector elements in <a>V</a>, possibly weighted by elements in <a>a</a>
  template< typename U>
  static U sum
    ( unsigned const n, U const* V, double const* a=nullptr );

  //! @brief Compute (symbolic) product of vector elements in <a>V</a>
  template< typename U>
  static U prod
    ( unsigned const n, U const* V );

  //! @brief Compute (symbolic) trace of a square matrix
  template< typename U>
  static U trace
    ( unsigned const n, U const* A );

  //! @brief Compute (symbolic) sum of two matrices
  template< typename U>
  static U* sum
    ( unsigned const m, unsigned const n,
      U const* A, U const* B );
  template< typename U>
  static void sum
    ( unsigned const m, unsigned const n, U const* A,
      U const* B, U* AB );

  //! @brief Compute (symbolic) difference of two matrices
  template< typename U>
  static U* sub
    ( unsigned const m, unsigned const n,
      U const* A, U const* B );
  template< typename U>
  static void sub
    ( unsigned const m, unsigned const n, U const* A,
      U const* B, U* AB );

  //! @brief Compute (symbolic) product of two matrices
  template< typename U>
  static U* prod
    ( unsigned const m, unsigned const n, unsigned const p,
      U const* A, U const* B );
  template< typename U>
  static void prod
    ( unsigned const m, unsigned  const n, unsigned const p,
      U const* A, U const* B, U* AB );

  //! @brief Compute (symbolic) determinant of a square matrix
  template< typename U>
  static U* polchar
    ( unsigned int const n, U const* A );
  template< typename U>
  static U det
    ( unsigned int const n, U const* A );
  template< typename U>
  static U* inv
    ( unsigned int const n, U const* A );
  /** @} */

protected:
  //! @brief Erase operation <a>op</a> in set <a>_Ops</a>
  bool _remove_operation
    ( FFOp*op );

  //! @brief Erase all operations in set <a>_Ops</a>
  void _clear_operations
    ()
    { for( auto&& op : _Ops ) delete op;
      _Ops.clear(); }

  //! @brief Reset all operations in set <a>_Ops</a>
  void _reset_operations
    ()
    const
    { for( auto&& op : _Ops ) op->iflag = 0; }

  //! @brief Looks for the operation of type <a>top</a> with operand <a>op</a> in set <a>_Ops</a> and adds it if not found
  FFOp* _insert_operation
    ( int const top, FFVar* op );

  //! @brief Looks for the operation of type <a>top</a> with left and right operands <a>lop</a>, <a>rop</a> in set <a>_Ops</a> and adds it if not found
  FFOp* _insert_operation
    ( int const top, FFVar* lop, FFVar* rop );

  //! @brief Looks for the operation of type <a>top</a> with <a>nop</a> operands <a>ops</a> in set <a>_Ops</a> and adds it if not found
  FFOp* _insert_operation
    ( int const top, unsigned const nop, FFVar** ops );

  //! @brief Looks for the n-ary operation of type <a>top</a> with operand array <a>pVar</a> of size <a>nVar</a> in set <a>_Ops</a> and adds it if not found; also adds new auxiliary variable in set <a>_Vars</a> and update list of dependencies in all operands in <a>_Vars</a>
  static FFVar& _insert_nary_operation
    ( int const top, FFDep const& dep, unsigned const nVar, FFVar const* pVar );

  //! @brief Looks for the binary operation of type <a>top</a> with left and right operands <a>Var1</a>, <a>Var2</a> in set <a>_Ops</a> and adds it if not found; also adds new auxiliary variable in set <a>_Vars</a> and update list of dependencies in both operands in <a>_Vars</a>
  static FFVar& _insert_binary_operation
    ( int const top, FFDep const& dep, FFVar const& Var1, FFVar const& Var2 );

  //! @brief Looks for the binary operation of type <a>top</a> with left and right operands <a>Cst1</a>, <a>Var2</a> in set <a>_Ops</a> and adds it if not found; also adds new auxiliary variable as well as constant <a>Cst1</a> in set <a>_Vars</a> and update list of dependencies in both operands in <a>_Vars</a>
  template <typename U> static FFVar&
  _insert_binary_operation
    ( int const top, FFDep const& dep, U const& Cst1, FFVar const& Var2 );

  //! @brief Looks for the binary operation of type <a>top</a> with left and right operands <a>Var1</a>, <a>Cst2</a> in set <a>_Ops</a> and adds it if not found; also adds new auxiliary variable as well as constant <a>Cst2</a> in set <a>_Vars</a> and update list of dependencies in both operands in <a>_Vars</a>
  template <typename U> static FFVar&
  _insert_binary_operation
    ( int const top, FFDep const& dep, FFVar const& Var1, U const& Cst2 );

  //! @brief Looks for the unary operation of type <a>top</a> with operand <a>Var</a> in set <a>_Ops</a> and adds it if not found; also adds new auxiliary variable in set <a>_Vars</a> and update list of dependencies in operand in <a>_Vars</a>
  static FFVar& _insert_unary_operation
    ( int const top, FFDep const& dep, FFVar const& Var );

  //! @brief Inserts the external n-ary operation <a>Op</a> with operand array <a>pVar</a> of size <a>nVar</a> in set <a>_Ops</a>, if not alrady present, adds new auxiliary variable in set <a>_Vars</a> and update list of dependencies in all operands in <a>_Vars</a>
  template <typename ExtOp>
  static FFVar& _insert_nary_external_operation
    ( ExtOp const& Op, FFDep const& dep, unsigned const nVar, FFVar const* pVar );

  //! @brief Inserts the external binary operation <a>Op</a> with operands <a>Var1</a> and <a>Var2</a> in set <a>_Ops</a>, if not alrady present, adds new auxiliary variable in set <a>_Vars</a> and update list of dependencies in all operands in <a>_Vars</a>
  template <typename ExtOp>
  static FFVar& _insert_binary_external_operation
    ( ExtOp const& Op, FFDep const& dep, FFVar const& Var1, FFVar const& Var2 );

  //! @brief Inserts the external unary operation <a>Op</a> with operand <a>Var</a> in set <a>_Ops</a>, if not alrady present, adds new auxiliary variable in set <a>_Vars</a> and update list of dependencies in all operands in <a>_Vars</a>
  template <typename ExtOp>
  static FFVar& _insert_unary_external_operation
    ( ExtOp const& Op, FFDep const& dep, FFVar const& Var1 );

  //! @brief Adds the auxiliary variable with dependency <a>dep</a> from operation <a>op</a>
  FFVar* _add_auxiliary
    ( FFDep const& dep, FFOp* pOp );

  //! @brief Looks for the real constant <a>x</a> and adds it if not found
  FFVar* _add_constant
    ( double const x );

  //! @brief Looks for the integer constant <a>n</a> and adds it if not found
  FFVar* _add_constant
    ( int const n );

  //! @brief Sets a variables to a constant and sets its numerical field
  FFVar* _set_constant
    ( FFVar const* pVar, FFNum const& num );

  //! @brief Unsets a variables to a constant
  FFVar* _unset_constant
    ( FFVar const* pVar );

  //! @brief Erase all variables in _Vars
  void _clear_variables
    ()
    { it_Vars itv = _Vars.begin();
      for( ; itv != _Vars.end(); ++itv ) delete *itv;
      _Vars.clear(); }

  //! @brief Appends the auxiliary variable <a>pAux</a> and define it in _Ops with type <a>tOp</a>
  void _append_aux
    ( FFVar* pAux, int const tOp );

  //! @brief Appends new auxiliary variable
  virtual void _append_aux
    ( FFVar* pAux );

  //! @brief Appends new original variable
  virtual void _append_var
    ( FFVar* pVar );

  //! @brief Search for the variable with identify <a>id</a> in <a>_Vars</a>
  FFVar* _find_var
    ( typename FFVar::pt_idVar const& id );

  //! @brief Create the variable with identify <a>id</a> and adds it if not found
  FFVar _create_var
    ( typename FFVar::pt_idVar const& id, std::string const& name );

  //! @brief Get DAG constant for variable <a>pVar</a> and adds it if not found; return false if <a>pVar</a> is not a constant
  bool _get_constant
    ( FFVar const*& pVar );

private:
  //! @brief Private methods to block default compiler methods
  FFBase
    ( FFBase const& );
  FFBase& operator=
    ( FFBase const& );
};

////////////////////////////////////////////////////////////////////////
// META-FUNCTIONS FOR MANIPULATING PARAMETER PACKS AND TUPLES

// Declare primary template
template<int I, typename... ExtOps>
struct nth_type_of
{
    using type = FFOp;
};

// Base step
template<typename ExtOp, typename... NextOps>
struct nth_type_of<0, ExtOp, NextOps...>
{
    using type = ExtOp;
};

// Induction step
template<int I, typename ExtOp, typename... NextOps>
struct nth_type_of<I, ExtOp, NextOps...>
{
    using type = typename nth_type_of<I - 1, NextOps...>::type;
};

// Helper meta-function for retrieving the first type in a parameter pack
template<typename... ExtOps>
struct first_type_of
{
    using type = typename nth_type_of<0, ExtOps...>::type;
};

// Helper meta-function for retrieving the last type in a parameter pack
template<typename... ExtOps>
struct last_type_of
{
    using type = typename nth_type_of<sizeof...(ExtOps) - 1, ExtOps...>::type;
};

// Base step
template <typename ExtOp>
struct remove_first_type
{
  typedef std::tuple<> type;
};

// Helper meta-function for removing the first type in a tuple
template <typename ExtOp, typename... NextOps>
struct remove_first_type< std::tuple<ExtOp, NextOps...> >
{
  typedef std::tuple<NextOps...> type;
};

//! @brief C++ class representing the DAG of factorable functions with external operations
////////////////////////////////////////////////////////////////////////
//! mc::FFGraph is a C++ class derived from mc::FFBase enabling
//! basic manipulations on a DAG, evaluation of a DAG in various
//! arithmetics, and external operations in a DAG.
////////////////////////////////////////////////////////////////////////
template <typename... ExtOps>
class FFGraph
: public FFBase
////////////////////////////////////////////////////////////////////////
{
  typedef std::tuple<ExtOps...> t_ExtOps;
  typedef typename remove_first_type< t_ExtOps >::type t_NextExtOps;
  typedef typename first_type_of< ExtOps... >::type FirstExtOps;

protected:
  //! @brief Work array for forward differentiation of subgraph
  std::vector< fadbad::F<FFVar> > _wkSFAD;

  //! @brief Work array for backward differentiation of subgraph
  std::vector< fadbad::B<FFVar> > _wkSBAD;

  //! @brief Work array for Taylor expansion of subgraph
  std::vector< fadbad::T<FFVar> > _wkTAD;

public:
  /** @ingroup FFunc
   *  @{
   */
  //! @brief Default Constructor
  FFGraph()
    : FFBase()
    {}

  //! @brief Destructor
  virtual ~FFGraph()
    {}

  //! @brief Expand DAG with derivatives of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> using fadbad::F. The return value is an array with entries of the dense Jacobian matrix. The function parameter pack <a>args</a> can be any number of extra pairs {const unsigned nIndep, const FFVar* const pIndep}, as well as a final, optional flag {const bool transp} indicating if the entries in the returned Jacobian matrix are ordered row-wise (transp=false, default) or column-wise (transp=true).
  template <typename... Deps>
  FFVar* FAD
    ( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
      const FFVar* const pIndep, Deps... args );

  //! @brief Expand DAG with derivatives of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> for the direction in array <a>pDir</a> using fadbad::F. The return value is an array with entries of the Jacobian matrix-vector product. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nIndep, const FFVar* const pIndep, const FFVar* const pDir}.
  template <typename... Deps> 
  FFVar* DFAD
    ( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
      const FFVar* const pIndep, const FFVar* const pDir, Deps... args );

  //! @brief Expand DAG with derivatives of dependents <a>vDep</a> with respect to independents <a>vIndep</a> using fadbad::F. The returns a vector with the entries of the Jacobian matrix, ordered either row-wise (transp=false) or column-wise (transp=true)
  std::vector<const FFVar*> FAD
    ( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
      const bool transp );

  //! @brief Expand DAG with derivatives of dependents <a>vDep</a> with respect to independents <a>vIndep</a> using fadbad::F (directional derivatives if <a>vDir</a> is specifed). The returns value is a vector with the entries of the Jacobian matrix, ordered row-wise
  std::vector<const FFVar*> FAD
    ( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
      const std::vector<const FFVar*>&vDir=std::vector<const FFVar*>() );

  //! @brief Expand DAG with derivatives of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> using fadbad::F. The return value is a 4-tuple of size and arrays with the row indices, column indices, and non-zero Jacobian entries. The function parameter pack <a>args</a> can be any number of pairs {const unsigned nIndep, const FFVar* const pIndep}, as well as a final, optional flag {const bool LUopt} indicating to only keep the entries in the lower (LUopt=true) or upper (LUopt=false) triangular part of the Jacobian matrix (e.g. for use in square symmetric matrix of derivatives such as Hessians).
  template <typename... Deps>
  std::tuple< unsigned, unsigned*, unsigned*, FFVar* > SFAD
    ( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
      const FFVar* const pIndep, Deps... args );

  //! @brief Expand DAG with derivatives of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> for the direction in array <a>pDir</a> using fadbad::F. The return value is a 4-tuple of size and arrays with the row indices, column indices, and non-zero Jacobian matrix-vector product entries. The function parameter pack <a>args</a> can be any number of triplets {const unsigned nIndep, const FFVar* const pIndep, const FFVar* const pDir}.
  template <typename... Deps>
  std::tuple< unsigned, unsigned*, unsigned*, FFVar* > SDFAD
    ( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
      const FFVar* const pIndep, const FFVar* const pDir, Deps... args );

  //! @brief Expand DAG with derivatives of dependents <a>vDep</a> with respect to independents <a>vIndep</a> using fadbad::F. The return value is a 3-tuple of vectors with the row indices, column indices, and non-zero Jacobian entries.
  std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> > SFAD
    ( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
      const std::vector<const FFVar*>&vDir=std::vector<const FFVar*>() );

  //! @brief Expand DAG with derivatives of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> using fadbad::B. The return value is an array with entries of the dense Jacobian matrix. The function parameter pack <a>args</a> can be any number of extra pairs {const unsigned nIndep, const FFVar* const pIndep}, as well as a final, optional flag {const bool transp} indicating if the entries in the returned Jacobian matrix are ordered row-wise (transp=false, default) or column-wise (transp=true).
  template <typename... Deps> 
  FFVar* BAD
    ( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
      const FFVar* const pIndep, Deps... args );

  //! @brief Expand DAG with derivatives of dependents <a>vDep</a> with respect to independents <a>vIndep</a> using fadbad::B (directional derivatives if <a>vDir</a> is specifed). The returns value is a vector with the entries of the Jacobian matrix, ordered row-wise
  std::vector<const FFVar*> BAD
    ( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
      const bool transp=false );

  //! @brief Expand DAG with derivatives of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> using fadbad::B. The return value is a 4-tuple of size and arrays with the row indices, column indices, and non-zero Jacobian entries. The function parameter pack <a>args</a> can be any number of extra pairs {const unsigned nIndep, const FFVar* const pIndep}, as well as a final, optional flag {const bool LUopt} indicating to only keep the entries in the lower (LUopt=true) or upper (LUopt=false) triangular part of the Jacobian matrix (e.g. for use in square symmetric matrix of derivatives such as Hessians).
  template <typename... Deps>
  std::tuple< unsigned, unsigned*, unsigned*, FFVar* > SBAD
    ( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
      const FFVar* const pIndep, Deps... args );

  //! @brief Expand DAG with derivatives of dependents <a>vDep</a> with respect to independents <a>vIndep</a> using fadbad::B. The return value is a 3-tuple of vectors with the row indices, column indices, and non-zero Jacobian entries.
  std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> > SBAD
    ( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep );

  //! @brief Expand DAG with Taylor coefficients of dependents <a>vDep</a> with respect to independents <a>vIndep</a> using fadbad::T -- Same number of dependents and independent is required, e.g. for expansion of ODE solutions -- Return a vector with the 0th, 1st, ..., ordermax'th order Taylor coefficients ordered sequentially
  std::vector<const FFVar*> TAD
    ( const unsigned int ordermax, const std::vector<const FFVar*>&vDep,
      const std::vector<const FFVar*>&vVar, const FFVar* const pIndep=0 );

  //! @brief Expand DAG with Taylor coefficients of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> using fadbad::T -- Same number of dependents and independent is required, e.g. for expansion of ODE solutions -- Returns an array with the 0th, 1st, ..., ordermax'th order Taylor coefficients ordered sequentially
  const FFVar* TAD
    ( const unsigned int ordermax, const unsigned nDep, const FFVar* const pDep,
      const unsigned nVar, const FFVar* const pVar, const FFVar* const pIndep=0 );

  //! @brief Insert the dependents <a>pDepIn</a> from the DAG <a>dag</a> into the current DAG with the resulting dependents <a>pDepOut</a>. Participating variables share the same indices in both DAGs
  template <typename DAG>
  void insert
    ( DAG* dag, unsigned const nDep, FFVar const* pDepIn, FFVar* pDepOut );

  //! @brief Compose the dependents in <a>vDepOut</a> with those in <a>vDepIn</a>. This function creates the subgraph for the outer dependent variables internally
  std::vector<const FFVar*> compose
    ( std::vector<const FFVar*>&vDepOut,
      std::vector< std::pair<const FFVar*, const FFVar*> >&vDepIn );

  //! @brief Compose the <a>nDepOut</a> dependents in array <a>pDepOut</a> with the <a>nDepIn</a> dependents in array <a>pDepIn</a> for the variables <a>pVarOut</a>. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nDepIn,      const FFVar*pVarOut, const FFVar*pDepIn}. This function creates the subgraph for the outer dependent variables internally
  template <typename... Deps>
  FFVar* compose
    ( unsigned const nDepOut, FFVar const* pDepOut, unsigned const nDepIn,
      FFVar const* pVarOut, FFVar const* pDepIn, Deps... args );

  //! @brief Compose the dependents indexed by <a>ndxDepOut</a> in array <a>pDepOut</a> with the <a>nDepIn</a> dependents in array <a>pDepIn</a> for the variables <a>pVarOut</a>. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nDepIn,      const FFVar*pVarOut, const FFVar*pDepIn}. This function creates the subgraph for the outer dependent variables internally
  template <typename... Deps>
  FFVar* compose
    ( std::set<unsigned> const& ndxDepOut, FFVar const* pDepOut, unsigned const nDepIn,
      FFVar const* pVarOut, FFVar const* pDepIn, Deps... args );
 
  //! @brief Evaluate the dependents in array <a>pDep</a> indexed by <a>ndxDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and write the result into <a>vDep</a>. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, const U*vVar}. This function allocates memory for intermediate operations internally. It creates the subgraph for the dependent variables internally. 
  template <typename U, typename... Deps> 
  void eval
    ( const std::set<unsigned>&ndxDep, const FFVar*pDep, U*vDep,
      const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args );

  //! @brief Evaluate the dependents in array <a>pDep</a> indexed by <a>ndxDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and write the result into <a>vDep</a>. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, const U*vVar}. This function allocates memory for intermediate operations internally.  It uses / creates the subgraph for the dependent variables passed via <a>sgDep</a> (e.g. to reduce the computational burden in repetitive evaluation of the same dependents/functions).
  template <typename U, typename... Deps> 
  void eval
    ( FFSubgraph&sgDep, const std::set<unsigned>&ndxDep, const FFVar*pDep, U*vDep,
      const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args );

  //! @brief Evaluate the dependents in the map <a>pDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and write the result in the map <a>vDep</a>. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, const U*vVar}. This function allocates memory for intermediate operations internally. It creates the subgraph for the dependent variables internally. 
  template <typename U, typename V, typename COMP, typename... Deps> 
  void eval
    ( const std::map<V,FFVar,COMP>&pDep, std::map<V,U,COMP>&vDep,
      const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args );

  //! @brief Evaluate the dependents in the map <a>pDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and write the result in the map <a>vDep</a>. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, const U*vVar}. This function allocates memory for intermediate operations internally.  It uses / creates the subgraph for the dependent variables passed via <a>sgDep</a> (e.g. to reduce the computational burden in repetitive evaluation of the same dependents/functions).
  template <typename U, typename V, typename COMP, typename... Deps> 
  void eval
    ( FFSubgraph&sgDep, const std::map<V,FFVar,COMP>&pDep, std::map<V,U,COMP>&vDep,
      const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and write the result into <a>vDep</a>. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, const U*vVar}, as well as a final, optional flag {const bool add} indicating if the dependent values are to overwrite (add=false) or be added to (add=true) <a>vDep</a>. This function allocates memory for intermediate operations internally. It creates the subgraph for the dependent variables internally. 
  template <typename U, typename... Deps> 
  void eval
    ( const unsigned nDep, const FFVar*pDep, U*vDep,
      const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and write the result into <a>vDep</a>. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, const U*vVar}, as well as a final, optional flag {const bool add} indicating if the dependent values are to overwrite (add=false) or be added to (add=true) <a>vDep</a>. This function allocates memory for intermediate operations internally. It uses / creates the subgraph for the dependent variables passed via <a>sgDep</a> (e.g. to reduce the computational burden in repetitive evaluation of the same dependents/functions).
  template <typename U, typename... Deps> 
  void eval
    ( FFSubgraph&sgDep, const unsigned nDep, const FFVar*pDep,
      U*vDep, const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> using the arithmetic U for the variable sizes and identifiers in lists <a>l_nVar</a> and <a>l_pVar</a>, whose values are specified in the list <a>l_vVar</a>, and write the result into <a>vDep</a>. The final, optional flag {const bool add} indicates if the dependent values are to overwrite (add=false) or be added to (add=true) those in <a>vDep</a>. This function allocates memory for intermediate operations internally. It uses / creates the subgraph for the dependent variables passed via <a>sgDep</a> (e.g. to reduce the computational burden in repetitive evaluation of the same dependents/functions).
  template <typename U> 
  void eval
    ( FFSubgraph&sgDep, const unsigned nDep, const FFVar*pDep,
      U*vDep, const std::list<unsigned>&nVar, const std::list<const FFVar*>&pVar,
      const std::list<const U*>&vVar, const bool add=false );

  //! @brief Evaluate the dependents in array <a>pDep</a> indexed by <a>ndxDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and write the result into <a>vDep</a>. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, const U*vVar}. This function stores the results of intermediate operations in the vector <a>wkDep</a>, resizing it as necessary. It creates the subgraph for the dependent variables internally. 
  template <typename U, typename... Deps> 
  void eval
    ( std::vector<U>&wkDep, const std::set<unsigned>&ndxDep, const FFVar*pDep, U*vDep,
      const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args );

  //! @brief Evaluate the dependents in array <a>pDep</a> indexed by <a>ndxDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and write the result into <a>vDep</a>. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, const U*vVar}. This function stores the results of intermediate operations in the vector <a>wkDep</a>, resizing it as necessary. It uses / creates the subgraph for the dependent variables passed via <a>sgDep</a> (e.g. to reduce the computational burden in repetitive evaluation of the same dependents/functions).
  template <typename U, typename... Deps> 
  void eval
    ( FFSubgraph&sgDep, std::vector<U>&wkDep, const std::set<unsigned>&ndxDep,
      const FFVar*pDep, U*vDep, const unsigned nVar, const FFVar*pVar,
      const U*vVar, Deps... args );

  //! @brief Evaluate the dependents in the map <a>pDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and write the result in the map <a>vDep</a>. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, const U*vVar}. This function stores the results of intermediate operations in the vector <a>wkDep</a>, resizing it as necessary. It creates the subgraph for the dependent variables internally. 
  template <typename U, typename V, typename COMP, typename... Deps> 
  void eval
    ( std::vector<U>&wkDep, const std::map<V,FFVar,COMP>&pDep, std::map<V,U,COMP>&vDep,
      const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args );

  //! @brief Evaluate the dependents in the map <a>pDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and write the result in the map <a>vDep</a>. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, const U*vVar}. This function stores the results of intermediate operations in the vector <a>wkDep</a>, resizing it as necessary. It uses / creates the subgraph for the dependent variables passed via <a>sgDep</a> (e.g. to reduce the computational burden in repetitive evaluation of the same dependents/functions).
  template <typename U, typename V, typename COMP, typename... Deps> 
  void eval
    ( FFSubgraph&sgDep, std::vector<U>&wkDep, const std::map<V,FFVar,COMP>&pDep,
      std::map<V,U,COMP>&vDep, const unsigned nVar, const FFVar*pVar, const U*vVar,
      Deps... args );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and write the result into <a>vDep</a>. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, const U*vVar}, as well as a final, optional flag {const bool add} indicating if the dependent values are to overwrite (add=false) or be added to (add=true) <a>vDep</a>. This function stores the results of intermediate operations in the vector <a>wkDep</a>, resizing it as necessary. It creates the subgraph for the dependent variables internally. 
  template <typename U, typename... Deps> 
  void eval
    ( std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
      U*vDep, const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and write the result into <a>vDep</a>. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, const U*vVar}, as well as a final, optional flag {const bool add} indicating if the dependent values are to overwrite (add=false) or be added to (add=true) <a>vDep</a>. This function stores the results of intermediate operations in the vector <a>wkDep</a>, resizing it as necessary. It uses / creates the subgraph for the dependent variables passed via <a>sgDep</a> (e.g. to reduce the computational burden in repetitive evaluation of the same dependents/functions).
  template <typename U, typename... Deps> 
  void eval
    ( FFSubgraph&sgDep, std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
      U*vDep, const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> using the arithmetic U for the variable sizes and identifiers in lists <a>l_nVar</a> and <a>l_pVar</a>, whose values are specified in the list <a>l_vVar</a>, and write the result into <a>vDep</a>. The final, optional flag {const bool add} indicates if the dependent values are to overwrite (add=false) or be added to (add=true) those in <a>vDep</a>. This function stores the results of intermediate operations in the vector <a>wkDep</a>, resizing it as necessary. It uses / creates the subgraph for the dependent variables passed via <a>sgDep</a> (e.g. to reduce the computational burden in repetitive evaluation of the same dependents/functions).
  template <typename U> 
  void eval
    ( FFSubgraph&sgDep, std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
      U*vDep, const std::list<unsigned>&l_nVar, const std::list<const FFVar*>&l_pVar,
      const std::list<const U*>&l_vVar, const bool add=false );

  //! @brief Evaluate the dependents in array <a>pDep</a> indexed by <a>ndxDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and use a priori information about the dependents in <a>vDep</a> to refine the variables in <a>pVar</a> based on forward/backard propagation. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, U*vVar}, as well as optional flags {const unsigned MAXPASS, const double THRESPASS} indicating the maximum number of forward/backward passes (default: 5) and minimum relative range reduction threshold (default: 0). This function allocates memory for intermediate operations internally. It also creates the subgraph for the dependent variables internally. The return value is the number of forward/backward passes, negative if the contraction leads to an empty intersection.
  template <typename U, typename... Deps> 
  int reval
    ( const std::set<unsigned>&ndxDep, const FFVar*pDep, U*vDep, 
      const unsigned nVar, const FFVar*pVar, U*vVar, Deps... args );

  //! @brief Evaluate the dependents in array <a>pDep</a> indexed by <a>ndxDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and use a priori information about the dependents in <a>vDep</a> to refine the variables in <a>pVar</a> based on forward/backard propagation. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, U*vVar}, as well as optional flags {const unsigned MAXPASS, const double THRESPASS} indicating the maximum number of forward/backward passes (default: 5) and minimum relative range reduction threshold (default: 0). This function stores the results of intermediate operations in the vector <a>wkDep</a>, resizing it as necessary. It creates the subgraph for the dependent variables internally. The return value is the number of forward/backward passes, negative if the contraction leads to an empty intersection.
  template <typename U, typename... Deps> 
  int reval
    ( std::vector<U>&wkDep, const std::set<unsigned>&ndxDep, const FFVar*pDep,
      U*vDep, const unsigned nVar, const FFVar*pVar, U*vVar, Deps... args );

  //! @brief Evaluate the dependents in array <a>pDep</a> indexed by <a>ndxDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and use a priori information about the dependents in <a>vDep</a> to refine the variables in <a>pVar</a> based on forward/backard propagation. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, U*vVar}, including optional flags {const unsigned MAXPASS, const double THRESPASS} indicating the maximum number of forward/backward passes (default: 5) and minimum relative range reduction threshold (default: 0). This function allocates memory for intermediate operations internally. It uses / creates the subgraph for the dependent variables passed via <a>sgDep</a> (e.g. to reduce the computational burden in repetitive evaluation of the same dependents/functions). The return value is the number of forward/backward passes, negative if the contraction leads to an empty intersection.
  template <typename U, typename... Deps> 
  int reval
    ( FFSubgraph&sgDep, const std::set<unsigned>&ndxDep, const FFVar*pDep, 
      U*vDep, const unsigned nVar, const FFVar*pVar, U*vVar, Deps... args );

  //! @brief Evaluate the dependents in array <a>pDep</a> indexed by <a>ndxDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and use a priori information about the dependents in <a>vDep</a> to refine the variables in <a>pVar</a> based on forward/backard propagation. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, U*vVar}, including optional flags {const unsigned MAXPASS, const double THRESPASS} indicating the maximum number of forward/backward passes (default: 5) and minimum relative range reduction threshold (default: 0). This function stores the results of intermediate operations in the vector <a>wkDep</a>, resizing it as necessary. It uses / creates the subgraph for the dependent variables passed via <a>sgDep</a> (e.g. to reduce the computational burden in repetitive evaluation of the same dependents/functions). The return value is the number of forward/backward passes, negative if the contraction leads to an empty intersection.
  template <typename U, typename... Deps> 
  int reval
    ( FFSubgraph&sgDep, std::vector<U>&wkDep, const std::set<unsigned>&ndxDep,
      const FFVar*pDep, U*vDep, const unsigned nVar, const FFVar*pVar,
      U*vVar, Deps... args );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and use a priori information about the dependents in <a>vDep</a> to refine the variables in <a>pVar</a> based on forward/backard propagation. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, U*vVar}, as well as optional flags {const unsigned MAXPASS, const double THRESPASS} indicating the maximum number of forward/backward passes (default: 5) and minimum relative range reduction threshold (default: 0). This function allocates memory for intermediate operations internally. It also creates the subgraph for the dependent variables internally. The return value is the number of forward/backward passes, negative if the contraction leads to an empty intersection.
  template <typename U, typename... Deps> 
  int reval
    ( const unsigned nDep, const FFVar*pDep, U*vDep, 
      const unsigned nVar, const FFVar*pVar, U*vVar, Deps... args );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and use a priori information about the dependents in <a>vDep</a> to refine the variables in <a>pVar</a> based on forward/backard propagation. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, U*vVar}, as well as optional flags {const unsigned MAXPASS, const double THRESPASS} indicating the maximum number of forward/backward passes (default: 5) and minimum relative range reduction threshold (default: 0). This function stores the results of intermediate operations in the vector <a>wkDep</a>, resizing it as necessary. It creates the subgraph for the dependent variables internally. The return value is the number of forward/backward passes, negative if the contraction leads to an empty intersection.
  template <typename U, typename... Deps> 
  int reval
    ( std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
      U*vDep, const unsigned nVar, const FFVar*pVar, U*vVar, Deps... args );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and use a priori information about the dependents in <a>vDep</a> to refine the variables in <a>pVar</a> based on forward/backard propagation. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, U*vVar}, including optional flags {const unsigned MAXPASS, const double THRESPASS} indicating the maximum number of forward/backward passes (default: 5) and minimum relative range reduction threshold (default: 0). This function allocates memory for intermediate operations internally. It uses / creates the subgraph for the dependent variables passed via <a>sgDep</a> (e.g. to reduce the computational burden in repetitive evaluation of the same dependents/functions). The return value is the number of forward/backward passes, negative if the contraction leads to an empty intersection.
  template <typename U, typename... Deps> 
  int reval
    ( FFSubgraph&sgDep, const unsigned nDep, const FFVar*pDep,
      U*vDep, const unsigned nVar, const FFVar*pVar, U*vVar, Deps... args );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> using the arithmetic U for the <a>nVar</a> variables in array <a>pVar</a>, whose values are specified in <a>vVar</a>, and use a priori information about the dependents in <a>vDep</a> to refine the variables in <a>pVar</a> based on forward/backard propagation. The function parameter pack <a>args</a> can be any number of extra triplets {const unsigned nVar, const FFVar*pVar, U*vVar}, including optional flags {const unsigned MAXPASS, const double THRESPASS} indicating the maximum number of forward/backward passes (default: 5) and minimum relative range reduction threshold (default: 0). This function stores the results of intermediate operations in the vector <a>wkDep</a>, resizing it as necessary. It uses / creates the subgraph for the dependent variables passed via <a>sgDep</a> (e.g. to reduce the computational burden in repetitive evaluation of the same dependents/functions). The return value is the number of forward/backward passes, negative if the contraction leads to an empty intersection.
  template <typename U, typename... Deps> 
  int reval
    ( FFSubgraph&sgDep, std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
      U*vDep, const unsigned nVar, const FFVar*pVar, U*vVar, Deps... args );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> using the arithmetic U for the variable sizes and identifiers in lists <a>l_nVar</a> and <a>l_pVar</a>, whose values are specified in the list <a>l_vVar</a>, and use a priori information about the dependents in <a>vDep</a> to refine the variables in <a>l_vVar</a> based on forward/backard propagation. Should the propagation fail for any operation, the default value <a>InfVal</a> is used instead. The optional flags {const unsigned MAXPASS, const double THRESPASS} indicate the maximum number of forward/backward passes (default: 5) and minimum relative range reduction threshold (default: 0). This function stores the results of intermediate operations in the vector <a>wkDep</a>, resizing it as necessary. It uses / creates the subgraph for the dependent variables passed via <a>sgDep</a> (e.g. to reduce the computational burden in repetitive evaluation of the same dependents/functions). The return value is the number of forward/backward passes, negative if the contraction leads to an empty intersection.
  template <typename U> 
  int reval
    ( FFSubgraph&sgDep, std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
      U*vDep, const std::list<unsigned>&l_nVar, const std::list<const FFVar*>&l_pVar,
      const std::list<U*>&l_vVar, U const& InfVal, const unsigned MAXPASS=5,
      const double THRESPASS=0. );

  //! @brief Extract operand values from work array <a>wkIn</a> corresponding to subgraph <a>sgIn</a> and copy them into work array <a>wkOut</a> corresponding to subgraph <a>sgOut</a>. This extraction assumes that the subgraph <a>sgOut</a> is contained within the subgraph <a>sgIn</a>, otherwise the behavior is undefined.
  template <typename U> 
  void wkextract
    ( const FFSubgraph&sgOut, std::vector<U>&wkOut, const FFSubgraph&sgIn, std::vector<U>&wkIn );

#ifdef MC__USE_HSL
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
#endif
  /** @} */

private:
  //! brief Intermediate function for recursive calls in FAD with a function parameter pack.
  template <typename... Deps> 
  std::vector<const FFVar*> FAD
    ( const std::vector<const FFVar*>&vDep, std::vector<const FFVar*>&vIndep,
      const unsigned nIndep, const FFVar* const pIndep, Deps... args );

  //! brief Intermediate function for recursive calls in directional FAD with a function parameter pack.
  template <typename... Deps> 
  std::vector<const FFVar*> FAD
    ( const std::vector<const FFVar*>&vDep, std::vector<const FFVar*>&vIndep,
      std::vector<const FFVar*>&vDir, const unsigned nIndep,
      const FFVar* const pIndep, const FFVar* const pDir, Deps... args );

  //! brief Intermediate function for recursive calls in sparse FAD with a function parameter pack.
  template <typename... Deps>
  std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> > SFAD
    ( const std::vector<const FFVar*>&vDep, std::vector<const FFVar*>&vIndep,
      const unsigned nIndep, const FFVar* const pIndep, Deps... args );

  //! brief Intermediate function for recursive calls in sparse FAD with a function parameter pack.
  std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> > SFAD
    ( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
      const bool LUopt );

  //! brief Intermediate function for recursive calls in sparse directional FAD with a function parameter pack.
  template <typename... Deps>
  std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> > SFAD
    ( const std::vector<const FFVar*>&vDep, std::vector<const FFVar*>&vIndep,
      std::vector<const FFVar*>&vDir, const unsigned nIndep,
      const FFVar* const pIndep, const FFVar* const pDir, Deps... args );

  //! brief Intermediate function for recursive calls in BAD with a function parameter pack.
  template <typename... Deps> 
  std::vector<const FFVar*> BAD
    ( const std::vector<const FFVar*>&vDep, std::vector<const FFVar*>&vIndep,
      const unsigned nIndep, const FFVar* const pIndep, Deps... args );

  //! brief Intermediate function for recursive calls in sparse BAD with a function parameter pack.
  template <typename... Deps>
  std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> > SBAD
    ( const std::vector<const FFVar*>&vDep, std::vector<const FFVar*>&vIndep,
      const unsigned nIndep, const FFVar* const pIndep, Deps... args );

  //! brief Intermediate function for recursive calls in sparse VBAD with a function parameter pack.
  std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> > SBAD
    ( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
      const bool LUopt );

  //! brief Intermediate function for recursive calls in DAG composition with a function parameter pack.
  template <typename... Deps> 
  std::vector<const FFVar*> compose
    ( std::vector<const FFVar*>&vDepOut,
      std::vector< std::pair<const FFVar*, const FFVar*> >&vDepIn,
      const unsigned nDepIn, const FFVar*pVarOut, const FFVar*pDepIn, Deps... args  );

  //! brief Intermediate function for recursive calls in DAG evaluation with a function parameter pack.
  template <typename U, typename... Deps> 
  void eval
    ( FFSubgraph&sgDep, const unsigned nDep, const FFVar*pDep, U*vDep, 
      std::list<unsigned>&l_nVar, std::list<const FFVar*>&l_pVar,
      std::list<const U*>&l_vVar, const unsigned nVar, const FFVar*pVar,
      const U*vVar, Deps... args );

  //! brief Intermediate function for recursive calls in DAG evaluation with a function parameter pack.
  template <typename U, typename... Deps> 
  void eval
    ( FFSubgraph&sgDep, std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
      U*vDep, std::list<unsigned>&l_nVar, std::list<const FFVar*>&l_pVar,
      std::list<const U*>&l_vVar, const unsigned nVar, const FFVar*pVar,
      const U*vVar, Deps... args );

  //! brief Intermediate function for recursive calls in DAG evaluation with a function parameter pack.
  template <typename U, typename... Deps> 
  int reval
    ( FFSubgraph&sgDep, std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
      U*vDep, std::list<unsigned>&l_nVar, std::list<const FFVar*>&l_pVar,
      std::list<U*>&l_vVar, const unsigned nVar, const FFVar*pVar,
      U*vVar, Deps... args );

private:
  //! @brief Private methods to block default compiler methods
  FFGraph
    ( FFGraph<ExtOps...> const& );
  FFGraph& operator=
    ( FFGraph<ExtOps...> const& );
};

////////////////////////////////// FFNum ///////////////////////////////////////

inline std::ostream&
operator<<
( std::ostream& out, FFNum const& Num )
{
  switch( Num.t ){
    case FFNum::INT:
      out << Num.n << "(I)"; break;
    case FFNum::REAL:
      out << Num.x << "(D)"; break;
  }
  return out;
}

////////////////////////////////// FFVar ///////////////////////////////////////

inline std::ostream&
operator<<
( std::ostream& out, FFVar const& Var )
{
  if( Var.id().second == FFVar::NOREF ) out << Var.num();
  else out << Var.name();
   // << " <= " << std::left << Var._num << "\t(" << Var._dag << ")";
  return out;
}

inline FFVar::FFVar
( FFBase* dag, std::string const& name )
: _dag( dag? dag: throw typename FFBase::Exceptions( FFBase::Exceptions::INIT )),
  _id( VAR, _dag->_nvar++ ), _num( 0./0. ), _dep(), _val( nullptr ),
  _mov( false ), _cst( false ), _nam( name )
{
  // Initialize dependence
  _dep.indep(_id.second);

  // Insert new variable in set FFBase::_Vars and corresponding operation in set FFBase::_Ops
  FFVar* pVar = new FFVar( *this );
  //std::cout << "name of inserted DAG variable: " << pVar->name() << std::endl;
  FFOp* pOp = new FFOp( FFOp::VAR, nullptr, pVar );
  _dag->_Ops.insert( pOp );
  pVar->_ops.first = _ops.first = pOp;
  _dag->_append_var( pVar );
}

inline FFVar::FFVar
( FFBase* dag, double const d )
: _dag( dag? dag: throw typename FFBase::Exceptions( FFBase::Exceptions::INIT )),
  _id( VAR, _dag->_nvar++ ), _num( d ), _dep(), _val( nullptr ),
  _mov( false ), _cst( true ), _nam( "" )
{ 
  // Initialize dependence
  _dep.indep(_id.second);

  // Insert new variable in set FFBase::_Vars and corresponding operation in set FFBase::_Ops
  FFVar* pVar = new FFVar( *this );
  FFOp* pOp = new FFOp( FFOp::VAR, nullptr, pVar );
  _dag->_Ops.insert( pOp );
  pVar->_ops.first = _ops.first = pOp;
  _dag->_append_var( pVar );
}

inline FFVar::FFVar
( FFBase* dag, int const i )
: _dag( dag? dag: throw typename FFBase::Exceptions( FFBase::Exceptions::INIT )),
  _id( VAR, _dag->_nvar++ ), _num( i ), _dep(), _val( nullptr ),
  _mov( false ), _cst( true ), _nam( "" )
{ 
  // Initialize dependence
  _dep.indep(_id.second);

  // Insert new variable in set FFBase::_Vars and corresponding operation in set FFBase::_Ops
  FFVar* pVar = new FFVar( *this );
  FFOp* pOp = new FFOp( FFOp::VAR, nullptr, pVar );
  _dag->_Ops.insert( pOp );
  pVar->_ops.first = _ops.first = pOp;
  _dag->_append_var( pVar );
}

inline FFVar::FFVar
( FFBase* dag, FFDep const& dep, FFOp* op )
: _dag( dag ), _id( AUX, dag->_naux++ ), _num( 0./0. ), _dep( dep ),
  _val ( nullptr ), _mov( false ), _cst( false ), _nam( "" )
{ _ops.first = op; }

inline FFVar::FFVar
( FFBase* dag, pt_idVar const& id, std::string const& name )
: _dag( dag? dag: throw typename FFBase::Exceptions( FFBase::Exceptions::INIT )),
  _id( id ), _num( 0./0. ), _dep(), _val( nullptr ), _mov( false ), _cst( false ),
  _nam( name )
{
  // Initialize dependence
  _dep.indep(_id.second);

  // Insert new variable in set FFBase::_Vars and corresponding operation in set FFBase::_Ops
  FFVar* pVar = new FFVar( *this );
  FFOp* pOp = new FFOp( FFOp::VAR, nullptr, pVar );
  _dag->_Ops.insert( pOp );
  pVar->_ops.first = _ops.first = pOp;
  _dag->_append_var( pVar );
}

inline void FFVar::set
( int const i )
const
{ _num = i; _cst = true; _mov = false;
  if( !_dag ) return;
  _dag->_set_constant( this, _num );
  return;
}

inline void FFVar::set
( double const d )
const
{ _num = d; _cst = true; _mov = false;
  if( !_dag ) return;
  _dag->_set_constant( this, _num );
  return;
}

inline void FFVar::unset
()
const
{ _cst = false; _mov = false;
  if( !_dag ) return;
  _dag->_unset_constant( this );
  return;
}

inline bool
FFVar::operator==
( FFVar const& Var )
const
{
  if( _dag != Var._dag || _id != Var._id || _cst ) return false;
  return true;
}

inline bool
FFVar::operator!=
( FFVar const& Var )
const
{
  if( _dag != Var._dag || _id != Var._id || _cst ) return true;
  return false;
}

inline FFVar&
FFVar::operator=
( FFVar const& Var )
{
  if( this == &Var ) return *this;
  _id   = Var._id;
  _num  = Var._num;
  _dep  = Var._dep;
  _dag  = Var._dag;
  _val  = Var._val;
  _mov  = Var._mov;
  _cst  = Var._cst;
  _ops  = Var._ops;
  _nam  = Var._nam;
  return *this;
}

inline FFVar&
FFVar::operator=
( int const i )
{
  _num  = i;
  _id   = std::make_pair(CINT,NOREF);
  _dep  = i;
  _dag  = nullptr;
  _val  = nullptr;
  _mov  = false;
  _cst  = true;
  _nam  = "";
  _ops.first = nullptr;
  _ops.second.clear();
  return *this;
}

inline FFVar&
FFVar::operator=
( double const x )
{
  _num  = x;
  _id   = std::make_pair((_num.t==FFNum::INT?CINT:CREAL),NOREF);
  _dep  = x;
  _dag  = nullptr;
  _val  = nullptr;
  _mov  = false;
  _cst  = true;
  _nam  = "";
  _ops.first = nullptr;
  _ops.second.clear();
  return *this;
}

inline FFVar
operator+
( FFVar const& Var )
{
  return Var;
}

template <typename U>
inline FFVar&
FFVar::operator+=
( U const& Var )
{
  FFVar VarNew = *this + Var;
  *this = VarNew;
  return *this;
}

inline FFVar
operator+
( FFVar const& Var1, FFVar const& Var2 )
{ 
  if( &Var1 == &Var2 || Var1 == Var2 ) return( 2. * Var1 );

  switch( Var1._id.first ){
  case FFVar::CREAL:
    switch( Var2._id.first ){
    case FFVar::CREAL:
      return( Var1._num.x + Var2._num.x );
    case FFVar::CINT:
      return( Var1._num.x + Var2._num.n );
    default:
      return( Var1._num.x + Var2 );
    }
  case FFVar::CINT:
    switch( Var2._id.first ){
    case FFVar::CREAL:
      return( Var1._num.n + Var2._num.x );
    case FFVar::CINT:
      return( Var1._num.n + Var2._num.n );
    default:
      return( Var1._num.n + Var2 );
    }
  default:
    switch( Var2._id.first ){
    case FFVar::CREAL:
      return( Var1 + Var2._num.x );
    case FFVar::CINT:
      return( Var1 + Var2._num.n );
    default:
      return FFBase::_insert_binary_operation( FFOp::PLUS, Var1._dep+Var2._dep, Var1, Var2 );
    }
  }
}

template <typename U>
inline FFVar
operator+
( U const& Cst, FFVar const& Var )
{
  // Case constant is zero
  if( Cst == U(0) ) return Var;

  switch( Var._id.first ){
  case FFVar::CREAL:
    return( Cst + Var._num.x );
  case FFVar::CINT:
    return( Cst + Var._num.n );
  default:
    return FFBase::_insert_binary_operation( FFOp::SHIFT, Cst+Var._dep, Var, (double)Cst );
  }
}

template <typename U>
inline FFVar
operator+
( FFVar const& Var, U const& Cst )
{
  return( Cst + Var );
}

inline FFVar
operator-
( FFVar const& Var )
{
  // Check if expression of type -(-X)
  if( Var._ops.first && Var._ops.first->type == FFOp::NEG )
    return *Var._ops.first->pops[0];

  switch( Var._id.first ){
  case FFVar::CREAL:
    return( -Var._num.x );
  case FFVar::CINT:
    return( -Var._num.n );
  default:
    return FFBase::_insert_unary_operation( FFOp::NEG, -Var._dep, Var );
  }
}

template <typename U>
inline FFVar&
FFVar::operator-=
( U const& Var )
{
  FFVar VarNew = *this - Var;
  *this = VarNew;
  return *this;
}

inline FFVar
operator-
( FFVar const& Var1, FFVar const& Var2 )
{
  if( &Var1 == &Var2 || Var1 == Var2 ) return 0.;

  switch( Var1._id.first ){
  case FFVar::CREAL:
    switch( Var2._id.first ){
    case FFVar::CREAL:
      return( Var1._num.x - Var2._num.x );
    case FFVar::CINT:
      return( Var1._num.x - Var2._num.n );
    default:
      return( Var1._num.x - Var2 );
    }
  case FFVar::CINT:
    switch( Var2._id.first ){
    case FFVar::CREAL:
      return( Var1._num.n - Var2._num.x );
    case FFVar::CINT:
      return( Var1._num.n - Var2._num.n );
    default:
      return( Var1._num.n - Var2 );
    }
  default:
    switch( Var2._id.first ){
    case FFVar::CREAL:
      return( Var1 - Var2._num.x );
    case FFVar::CINT:
      return( Var1 - Var2._num.n );
    default:
      return FFBase::_insert_binary_operation( FFOp::MINUS, Var1._dep-Var2._dep, Var1, Var2 );
    }
  }
}

template <typename U>
inline FFVar
operator-
( FFVar const& Var, U const& Cst )
{
  // Case constant is zero
  if( Cst == U(0) ) return Var;

  switch( Var._id.first ){
  case FFVar::CREAL:
    return( Var._num.x - Cst );
  case FFVar::CINT:
    return( Var._num.n - Cst );
  default:
    return FFBase::_insert_binary_operation( FFOp::SHIFT, Var._dep-Cst, Var, -(double)Cst );
  }
}

template <typename U>
inline FFVar
operator-
( U const& Cst, FFVar const& Var )
{
  return( Cst + (-Var) );
}

template <typename U>
inline FFVar&
FFVar::operator*=
( U const& Var )
{
  FFVar VarNew = *this * Var;
  *this = VarNew;
  return *this;
}

inline FFVar
operator*
( FFVar const& Var1, FFVar const& Var2 )
{
  if( &Var1 == &Var2 || Var1 == Var2 ) return sqr(Var1);

  switch( Var1._id.first ){
  case FFVar::CREAL:
    switch( Var2._id.first ){
    case FFVar::CREAL:
      return( Var1._num.x * Var2._num.x );
    case FFVar::CINT:
      return( Var1._num.x * Var2._num.n );
    default:
      return( Var1._num.x * Var2 );
    }
  case FFVar::CINT:
    switch( Var2._id.first ){
    case FFVar::CREAL:
      return( Var1._num.n * Var2._num.x );
    case FFVar::CINT:
      return( Var1._num.n * Var2._num.n );
    default:
      return( Var1._num.n * Var2 );
    }
  default:
    switch( Var2._id.first ){
    case FFVar::CREAL:
      return( Var1 * Var2._num.x );
    case FFVar::CINT:
      return( Var1 * Var2._num.n );
    default:{
      return FFBase::_insert_binary_operation( FFOp::TIMES, Var1._dep*Var2._dep, Var1, Var2 );
     }
    }
  }
}

template <typename U>
inline FFVar
operator*
( FFVar const& Var, U const& Cst )
{
  return( Cst * Var );
}

template <typename U>
inline FFVar
operator*
( U const& Cst, FFVar const& Var )
{
  // Case constant is zero
  if( Cst == U(0) ) return 0.;
  // Case constant is one
  if( Cst == U(1) ) return Var;
  // Case constant is negative one
  if( Cst == U(-1) ) return -Var;

  switch( Var._id.first ){
  case FFVar::CREAL:
    return( Cst * Var._num.x );
  case FFVar::CINT:
    return( Cst * Var._num.n );
  default:
    return FFBase::_insert_binary_operation( FFOp::SCALE, Cst*Var._dep, Var, (double)Cst );
  }
}

template <typename U>
inline FFVar&
FFVar::operator/=
( U const& Var )
{
  FFVar VarNew = *this / Var;
  *this = VarNew;
  return *this;
}

inline FFVar
operator/
( FFVar const& Var1, FFVar const& Var2 )
{
  if( &Var1 == &Var2 || Var1 == Var2 ) return 1.;
  // UNSURE WHAT THIS IS DOING???
  if( Var2._id.first == FFVar::CREAL || Var2._id.first == FFVar::CINT )
    return std::numeric_limits<double>::quiet_NaN();

  switch( Var1._id.first ){
  case FFVar::CREAL:
    switch( Var2._id.first ){
    case FFVar::CREAL:
      return( Var1._num.x / Var2._num.x );
    case FFVar::CINT:
      return( Var1._num.x / Var2._num.n );
    default:
      return( Var1._num.x / Var2 );
    }
  case FFVar::CINT:
    switch( Var2._id.first ){
    case FFVar::CREAL:
      return( Var1._num.n / Var2._num.x );
    case FFVar::CINT:
      return( Var1._num.n / Var2._num.n );
    default:
      return( Var1._num.n / Var2 );
    }
  default:
    switch( Var2._id.first ){
    case FFVar::CREAL:
      return( Var1 / Var2._num.x );
    case FFVar::CINT:
      return( Var1 / Var2._num.n );
    default:{
      return FFBase::_insert_binary_operation( FFOp::DIV, Var1._dep/Var2._dep, Var1, Var2 );
     }
    }
  }
}

template <typename U>
inline FFVar
operator/
( FFVar const& Var, U const& Cst )
{
  if( Cst == U(0) ) return std::numeric_limits<double>::quiet_NaN();
  return( ( 1. / Cst ) * Var );
}

template <typename U>
inline FFVar
operator/
( U const& Cst, FFVar const& Var )
{
  // Case constant is zero
  if( Cst == U(0) ) return 0.;

  switch( Var._id.first ){
  case FFVar::CREAL:
    return( Cst / Var._num.x );
  case FFVar::CINT:
    return( Cst / Var._num.n );
  default:{
    return FFBase::_insert_binary_operation( FFOp::INV, Cst/Var._dep, (double)Cst, Var );
   }
  }
}

inline FFVar
inv
( FFVar const& Var )
{
  return( 1. / Var );
}

inline FFVar
sqr
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( mc::sqr( Var._num.n ) );
      case FFNum::REAL:  return( mc::sqr( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::SQR, sqr(Var._dep), Var );
}

inline FFVar
sqrt
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::sqrt( Var._num.n ) );
      case FFNum::REAL:  return( std::sqrt( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::SQRT, sqrt(Var._dep), Var );
}

inline FFVar
pow
( FFVar const& Var, int const iExp )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::pow( Var._num.n, iExp ) );
      case FFNum::REAL:  return( std::pow( Var._num.x, iExp ) );
    }
  }

  // Case integer exponent is 0,1, or 2
  if( iExp == 0 )  return( 1. );
  if( iExp == 1 )  return Var;
  if( iExp == 2 )  return sqr(Var);
  if( iExp == -1 ) return 1./Var;

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant iExp if not defined
  return FFBase::_insert_binary_operation( FFOp::IPOW, pow(Var._dep,iExp), Var, iExp );
}

inline FFVar
pow
( FFVar const& Var, double const dExp )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::pow( Var._num.n, dExp ) );
      case FFNum::REAL:  return( std::pow( Var._num.x, dExp ) );
    }
  }

  // Case exponent is negative: compute 1./(Var^dExp)
  if( dExp<0. )  return 1./( pow(Var,-dExp) );

  // Case integer exponent is 0,1, or 2
  if( dExp==0. )  return( 1. );
  if( isequal(dExp,1.) )  return Var;
  if( isequal(dExp,2.) )  return sqr(Var);
  if( isequal(dExp,0.5) ) return sqrt(Var);
  if( std::floor(dExp)==dExp && dExp>=INT_MIN && dExp<=INT_MAX )
    return pow( Var, (int)dExp );

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant dExp if not defined
  return FFBase::_insert_binary_operation( FFOp::DPOW, pow(Var._dep,dExp), Var, dExp );
  //return exp( dExp * log( Var ) );
}

inline FFVar
pow
( FFVar const& Var1, FFVar const& Var2 )
{
  // Case exponent is integer constant
  if( Var2._id.second == FFVar::NOREF && Var2._num.t == FFNum::INT )
    return pow( Var1, Var2._num.n );

  return exp( Var2 * log( Var1 ) );
}

inline FFVar
pow
( double const Cst1, FFVar const& Var2 )
{
  // Case exponent is integer constant
  if( Var2._id.second == FFVar::NOREF && Var2._num.t == FFNum::INT )
    return std::pow( Cst1, Var2._num.n );

  return exp( Var2 * std::log( Cst1 ) );
}

inline FFVar
cheb
( FFVar const& Var, unsigned const iOrd )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( mc::cheb( Var._num.n, iOrd ) );
      case FFNum::REAL:  return( mc::cheb( Var._num.x, iOrd ) );
    }
  }

  // Case integer exponent is 0,1, or 2
  switch( iOrd ){
    case 0: return( 1. );
    case 1: return Var;
    case 2: return 2.*sqr(Var)-1.;
    default: break;
  }
  
  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant iOrd if not defined
  if( !Var._dag || !Var._dag->options.CHEBRECURS )
    return FFBase::_insert_binary_operation( FFOp::CHEB, cheb(Var._dep,iOrd), Var, (int)iOrd );

  if( iOrd==3 ){
    FFVar VarRecu = (4.*sqr(Var)-3.)*Var;
    FFVar VarCheb = FFBase::_insert_binary_operation( FFOp::CHEB, cheb(Var._dep,iOrd), Var, (int)iOrd );
    return inter( VarCheb, VarRecu );
  }
    
  FFVar VarRecu = inter( 2.*Var*cheb(Var,iOrd-1)-cheb(Var,iOrd-2),
                         iOrd%2? 2.*cheb(Var,iOrd/2)*cheb(Var,iOrd/2+1)-Var: 2.*sqr(cheb(Var,iOrd/2))-1. );
  FFVar VarCheb = FFBase::_insert_binary_operation( FFOp::CHEB, cheb(Var._dep,iOrd), Var, (int)iOrd );
  return inter( VarCheb, VarRecu );
}

inline FFVar
exp
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::exp( Var._num.n ) );
      case FFNum::REAL:  return( std::exp( Var._num.x ) );
    }
  }

  // Case exponent is negative: compute 1./(Var^dExp)
  if( Var._dag->options.DETECTSIGNOM
   && Var._ops.first
   && Var._ops.first->type == FFOp::SCALE
   && Var._ops.first->pops[0]->_ops.first
   && Var._ops.first->pops[0]->_ops.first->type == FFOp::LOG ){
    auto const* Aux = Var._ops.first->pops[0]->_ops.first->pops[0];
    auto const* Cst = Var._ops.first->pops[1];
#ifdef MC__FFUNC_DEBUG_SIGNOM
    std::cout << "Signomial term detected: exp(" << Cst->name() << ".log(" << Aux->name() << ")" << std::endl;
#endif
    switch( Cst->_num.t ){
      case FFNum::INT:   return( pow( *Aux, Cst->_num.n ) );
      case FFNum::REAL:  return( pow( *Aux, Cst->_num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::EXP, exp(Var._dep), Var );
}

inline FFVar
log
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::log( Var._num.n ) );
      case FFNum::REAL:  return( std::log( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::LOG, log(Var._dep), Var );
}

inline FFVar
xlog
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( Var._num.n*std::log( Var._num.n ) );
      case FFNum::REAL:  return( Var._num.x*std::log( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::XLOG, xlog(Var._dep), Var );
}

inline FFVar
cos
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::cos( Var._num.n ) );
      case FFNum::REAL:  return( std::cos( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::COS, cos(Var._dep), Var );
}

inline FFVar
sin
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::sin( Var._num.n ) );
      case FFNum::REAL:  return( std::sin( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::SIN, sin(Var._dep), Var );
}

inline FFVar
tan
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::tan( Var._num.n ) );
      case FFNum::REAL:  return( std::tan( Var._num.x ) );

    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::TAN, tan(Var._dep), Var );
}

inline FFVar
asin
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::asin( Var._num.n ) );
      case FFNum::REAL:  return( std::asin( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::ASIN, asin(Var._dep), Var );
}

inline FFVar
acos
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::acos( Var._num.n ) );
      case FFNum::REAL:  return( std::acos( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::ACOS, acos(Var._dep), Var );
}

inline FFVar
atan
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::atan( Var._num.n ) );
      case FFNum::REAL:  return( std::atan( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::ATAN, atan(Var._dep), Var );
}

inline FFVar
sinh
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::sinh( Var._num.n ) );
      case FFNum::REAL:  return( std::sinh( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::SINH, sinh(Var._dep), Var );
}

inline FFVar
cosh
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::cosh( Var._num.n ) );
      case FFNum::REAL:  return( std::cosh( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::COSH, cosh(Var._dep), Var );
}

inline FFVar
tanh
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::tanh( Var._num.n ) );
      case FFNum::REAL:  return( std::tanh( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::TANH, tanh(Var._dep), Var );
}

inline FFVar
erf
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( ::erf( Var._num.n ) );
      case FFNum::REAL:  return( ::erf( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::ERF, erf(Var._dep), Var );
}

inline FFVar
erfc
( FFVar const& Var )
{
  return ( 1. - erf( Var ) );
}

inline FFVar
fabs
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::fabs( Var._num.n ) );
      case FFNum::REAL:  return( std::fabs( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::FABS, fabs(Var._dep), Var );
}

inline FFVar
fstep
( FFVar const& Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( mc::fstep( Var._num.n ) );
      case FFNum::REAL:  return( mc::fstep( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_unary_operation( FFOp::FSTEP, fstep(Var._dep), Var );
}

inline FFVar
bstep
( FFVar const& Var )
{
  return ( fstep( -Var ) );
}

inline FFVar
max
( FFVar const& Var1, FFVar const& Var2 )
{
  if( &Var1 == &Var2 || Var1 == Var2 ) return Var1;
  //if( Var1 == Var2 ) return Var1;

  // Case either or both operands are numeric constants
  if( Var1._id.second == FFVar::NOREF && Var2._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.n>Var2._num.n? Var1._num.n: Var2._num.n );
        case FFNum::REAL:  return( std::max((double)Var1._num.n,Var2._num.x) );
      }
      case FFNum::REAL:
      switch( Var2._num.t ){
        case FFNum::INT:   return( std::max(Var1._num.x,(double)Var2._num.n) );
        case FFNum::REAL:  return( std::max(Var1._num.x,Var2._num.x) );
      }
    }
  }

  if( Var1._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:   return( max((double)Var1._num.n,Var2) );
      case FFNum::REAL:  return( max(Var1._num.x,Var2) );
    }
  }
  
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( max((double)Var2._num.n,Var1) );
      case FFNum::REAL:  return( max(Var2._num.x,Var1) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_binary_operation( FFOp::MAXF, max(Var1._dep,Var2._dep), Var1, Var2 );
}

inline FFVar
max
( FFVar const& Var1, const double&Cst2 )
{
  return( max( Cst2, Var1 ) );
}

inline FFVar
max
( const double&Cst1, FFVar const& Var2 )
{
  // Case right operand is a numeric constant
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( std::max(Cst1,(double)Var2._num.n) );
      case FFNum::REAL:  return( std::max(Cst1,Var2._num.x) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant Cst1 if not defined
  return FFBase::_insert_binary_operation( FFOp::MAXF, max((double)Cst1,Var2._dep), (double)Cst1, Var2 );
}

inline FFVar
max
( const unsigned nVar, const FFVar*pVar )
{
  if( !nVar || !pVar ) return( 0 );
  
  FFVar VarR = pVar[0];
  for( unsigned int i=1; i<nVar; i++ ) VarR = max( VarR, pVar[i] );
  return( VarR );
}

inline FFVar
min
( FFVar const& Var1, FFVar const& Var2 )
{
  if( &Var1 == &Var2 || Var1 == Var2 ) return Var1;
  //if( Var1 == Var2 ) return Var1;

  // Case either or both operands are numeric constants
  if( Var1._id.second == FFVar::NOREF && Var2._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.n<Var2._num.n? Var1._num.n: Var2._num.n );
        case FFNum::REAL:  return( std::min((double)Var1._num.n,Var2._num.x) );
      }
      case FFNum::REAL:
      switch( Var2._num.t ){
        case FFNum::INT:   return( std::min(Var1._num.x,(double)Var2._num.n) );
        case FFNum::REAL:  return( std::min(Var1._num.x,Var2._num.x) );
      }
    }
  }

  if( Var1._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:   return( min((double)Var1._num.n,Var2) );
      case FFNum::REAL:  return( min(Var1._num.x,Var2) );
    }
  }
  
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( min((double)Var2._num.n,Var1) );
      case FFNum::REAL:  return( min(Var2._num.x,Var1) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_binary_operation( FFOp::MINF, min(Var1._dep,Var2._dep), Var1, Var2 );
}

inline FFVar
min
( FFVar const& Var1, const double&Cst2 )
{
  return( min( Cst2, Var1 ) );
}

inline FFVar
min
( const double&Cst1, FFVar const& Var2 )
{
  // Case right operand is a numeric constant
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( std::min(Cst1,(double)Var2._num.n) );
      case FFNum::REAL:  return( std::min(Cst1,Var2._num.x) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant Cst1 if not defined
  return FFBase::_insert_binary_operation( FFOp::MINF, min((double)Cst1,Var2._dep), (double)Cst1, Var2 );
}

inline FFVar
min
( const unsigned nVar, const FFVar*pVar )
{
  if( !nVar || !pVar ) return( 0 );
  
  FFVar VarR = pVar[0];
  for( unsigned int i=1; i<nVar; i++ ) VarR = min( VarR, pVar[i] );
  return( VarR );
}

inline FFVar
inter
( FFVar const& Var1, FFVar const& Var2 )
{
  if( &Var1 == &Var2 || Var1 == Var2 ) return Var1;
  //if( Var1 == Var2 ) return Var1;

  // Case either or both operands are numeric constants
  if( Var1._id.second == FFVar::NOREF && Var2._id.second == FFVar::NOREF ){
    if( Var1._num.val() != Var2._num.val() )
      throw typename FFBase::Exceptions( FFBase::Exceptions::INTER ); 
    return Var1._num.val();
  }
  if( Var1._id.second == FFVar::NOREF )
    return inter( Var2, Var1._num.val() );
  if( Var2._id.second == FFVar::NOREF )
    return inter( Var1, Var2._num.val() );

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_binary_operation( FFOp::INTER, inter(Var1._dep,Var2._dep), Var1, Var2 );
}

template <typename U> inline FFVar
inter
( U const& Cst1, FFVar const& Var2 )
{
  // Case right operand is a numerhic constant
  if( Var2._id.second == FFVar::NOREF ){
    if( Cst1 != Var2._num.val() )
      throw typename FFBase::Exceptions( FFBase::Exceptions::INTER ); 
    return Cst1;
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFBase::_insert_binary_operation( FFOp::INTER, inter((double)Cst1,Var2._dep), (double)Cst1, Var2 );
}

template <typename U>
inline FFVar
inter
( FFVar const& Var1, U const& Cst2 )
{
  return( inter( Cst2, Var1 ) );
}

inline FFVar
sum
( const unsigned nVar, const FFVar*pVar )
{
  switch( nVar ){
   case 0:  return( 0 );
   case 1:  return( pVar[0] );
   default: return( pVar[0] + sum( nVar-1, pVar+1 ) );
  }
}

inline FFVar
prod
( unsigned const nVar, FFVar const* pVar )
{
  // Case 0 operand
  if( !nVar || !pVar ) return( 1 );
  // Case 1 operand
  if( nVar == 1 ) return( pVar[0] );
  // Case 2 operands
  if( nVar == 2 ) return( pVar[0] * pVar[1] );
  // Case >2 operands
  double Scal = 1.;
  std::vector<FFDep> vDep;
  std::vector<FFVar> vVar;
  for( unsigned i=0; i<nVar; i++ ){ 
    // Constant operands are removed from multi-linear term
    if( pVar[i]._id.first == FFVar::CINT || pVar[i]._id.first == FFVar::CREAL ){
      Scal *= pVar[i]._num.val();
      continue;
    }
    vDep.push_back( pVar[i]._dep );
    vVar.push_back( pVar[i] );
  }

  switch( vVar.size() ){
    case 0:  return( Scal );
    case 1:  return( Scal * vVar[0] );
    case 2:  return( Scal * ( vVar[0] * vVar[1] ) );
    default: break;
  }
  FFVar VarProd = FFBase::_insert_nary_operation( FFOp::PROD, prod( vDep.size(), vDep.data() ), vVar.size(), vVar.data() );
  return( vVar.size() < nVar? Scal * VarProd: VarProd );
}

inline FFVar
monom
( unsigned const nVar, FFVar const* pVar, unsigned const* pExp, bool const cmon=false )
{
  // Case 0 operand
  if( !nVar || !pVar ) return( 1 );
  // Case 1 operand
  if( nVar == 1 ) return( cmon? cheb( pVar[0], pExp[0] ): pow( pVar[0], (int)pExp[0] ) );
  // Case 2 operands
  if( nVar == 2 ) return( cmon? cheb( pVar[0], pExp[0] ) * cheb( pVar[1], pExp[1] ):
                                pow( pVar[0], (int)pExp[0] ) * pow( pVar[1], (int)pExp[1] ) );
  // Case >2 operands
  double Scal = 1.;
  std::vector<FFDep> vDep;
  std::vector<FFVar> vVar;
  for( unsigned i=0; i<nVar; i++ ){ 
    // Constant operands are removed from monomial term
    if( pVar[i]._id.first == FFVar::CINT || pVar[i]._id.first == FFVar::CREAL ){
      Scal *= ( cmon? cheb( pVar[i]._num.val(), pExp[i] ): std::pow( pVar[i]._num.val(), (int)pExp[i] ) );
      continue;
    }
    FFVar term = ( cmon? cheb( pVar[i], pExp[i] ): pow( pVar[i], (int)pExp[i] ) );
    vDep.push_back( term._dep );
    vVar.push_back( term );
  }
  switch( vVar.size() ){
    case 0:  return( Scal );
    case 1:  return( Scal * vVar[0] );
    case 2:  return( Scal * ( vVar[0] * vVar[1] ) );
    default: break;
  }
  FFVar VarProd = FFBase::_insert_nary_operation( FFOp::PROD, prod( vDep.size(), vDep.data() ), vVar.size(), vVar.data() );
  return( vVar.size() < nVar? Scal * VarProd: VarProd );
}

////////////////////////////////// FFOp ////////////////////////////////////////

inline std::ostream&
operator <<
( std::ostream& out, FFOp const& Op)
{
  switch( Op.type ){
    case FFOp::CNST:  out << Op.pres->num() << "\t";
                      break;
    case FFOp::VAR:   if( Op.pres->cst() ) out << Op.pres->num(); else out << "VARIABLE";
                      break;
    case FFOp::NEG:   out << "- " << Op.pops[0]->name() << "\t" ;
                      break;
    case FFOp::PLUS:
    case FFOp::SHIFT: 
    case FFOp::MINUS: 
    case FFOp::TIMES:
    case FFOp::SCALE: 
    case FFOp::DIV:
    case FFOp::INV:   out << Op.pops[0]->name() << Op.name() << Op.pops[1]->name() << "\t";
                      break;
    case FFOp::IPOW:
    case FFOp::DPOW:
    case FFOp::CHEB:
    case FFOp::MINF:
    case FFOp::MAXF:
    case FFOp::INTER: out << Op.name() << "( " << Op.pops[0]->name() << ", "
                                               << Op.pops[1]->name() << " )";
                      break;
    case FFOp::SQR:
    case FFOp::SQRT:
    case FFOp::EXP:
    case FFOp::LOG:
    case FFOp::XLOG:
    case FFOp::COS:
    case FFOp::SIN:
    case FFOp::TAN:
    case FFOp::ASIN:
    case FFOp::ACOS:
    case FFOp::ATAN:
    case FFOp::SINH:
    case FFOp::COSH:
    case FFOp::TANH:
    case FFOp::ERF:
    case FFOp::FABS:
    case FFOp::FSTEP: out << Op.name() << "( " << Op.pops[0]->name() << " )\t";
                      break;
    case FFOp::PROD:  out << Op.name() << "( ";
                      for( unsigned i=0; i<Op.pops.size()-1; i++ ) out << Op.pops[i]->name() << ", ";
                      out << Op.pops[Op.pops.size()-1]->name() << " )\t";
                      break;
    default: 
      if( Op.type >= FFOp::EXTERN ){
                      out << Op.name() << "( ";
                      for( unsigned i=0; i<Op.pops.size()-1; i++ ) out << Op.pops[i]->name() << ", ";
                      out << Op.pops[Op.pops.size()-1]->name() << " )\t";
                      break;
      }
      // Should not reach this point
      throw typename FFBase::Exceptions( FFBase::Exceptions::INTERN );
  } 
  return out;
}

inline
FFOp::FFOp
( int const top, FFVar* lop, FFVar* res ):
  type( top ), pres( res ), pops(), iflag(0)
{
  if( lop ) pops.push_back( lop );
}

inline
FFOp::FFOp
( int const top, FFVar* lop, FFVar* rop, FFVar* res ):
  type( top ), pres( res ), pops(), iflag(0)
{
  // Reorder operands in commutative operations
  if( lop && commutative() && lt_FFVar()( rop, lop ) )
    { pops.push_back( rop ); pops.push_back( lop ); }
  else
    { pops.push_back( lop ); pops.push_back( rop ); }
}

inline
FFOp::FFOp
( int const top, unsigned const nop, FFVar** ops, FFVar* res ):
  type( top ), pres( res ), pops( ops, ops+nop ), iflag(0)
{
  if( nop < 2 || !commutative() ) return;
  std::sort( pops.begin(), pops.end(), lt_FFVar() );
}

inline FFOp&
FFOp::set
( FFVar* lop, FFVar* res )
{
  pres = res;
  pops.assign( { lop } );
  return *this;
}

inline FFOp&
FFOp::set
( FFVar* lop, FFVar* rop, FFVar* res )
{
  pres = res;
  // Reorder operands in commutative operations
  if( commutative() && lt_FFVar()( rop, lop ) )
    pops.assign( { rop, lop } );
  else
    pops.assign( { lop, rop } );
  return *this;
}

inline FFOp&
FFOp::set
( unsigned const nop, FFVar** ops, FFVar* res )
{
  pres = res;
  pops.assign( ops, ops+nop );
  if( nop < 2 || !commutative() ) return *this;
  std::sort( pops.begin(), pops.end(), lt_FFVar() );
  return *this;
}

template <typename ExtOp>
inline FFVar&
FFOp::insert_external_operation
( ExtOp const& Op, FFDep const& dep, unsigned const nVar, FFVar const* pVar )
{
  return FFBase::_insert_nary_external_operation( Op, dep, nVar, pVar );
}

template <typename ExtOp>
inline FFVar&
FFOp::insert_external_operation
( ExtOp const& Op, FFDep const& dep, FFVar const& Var1, FFVar const& Var2 )
{
  return FFBase::_insert_binary_external_operation( Op, dep, Var1, Var2 );
}

template <typename ExtOp>
inline FFVar&
FFOp::insert_external_operation
( ExtOp const& Op, FFDep const& dep, FFVar const& Var )
{
  return FFBase::_insert_unary_external_operation( Op, dep, Var );
}

inline void
FFOp::propagate_subgraph
( std::list< std::pair< FFOp const*, bool > >& Ops )
const
{
  if( iflag ){
    // indicate non-movability of current operation result 
    pres->mov() = false;
    return;
  }

  for( auto it=pops.begin(); it!=pops.end(); ++it ){
    if( !(*it) || !(*it)->ops().first ) continue;
    (*it)->ops().first->propagate_subgraph( Ops );
  }

  Ops.push_back( std::make_pair( this, false ) );
  pres->mov() = true;
  iflag = Ops.size();
}

template <typename U> inline void
FFOp::reset_val_subgraph
( U const&  U_dum )
const
{
  if( iflag ) return;
  iflag = 1;

  for( auto it=pops.begin(); it!=pops.end(); ++it ){
    if( !(*it) || !(*it)->ops().first ) continue;
    (*it)->ops().first->reset_val_subgraph( U_dum );
  }
  if( !pres ) return;
  if( pres->val() ) pres->reset_val( U_dum );
  pres->mov() = false;
}

template <typename U> inline void
FFOp::evaluate
( typename std::vector<U>::iterator itU, U const* pU_dum, bool const movU )
const
{
  pres->val()  = &(*itU);
  pres->mov() = movU;
  //std::cout << "evaluation of " << *pres << ": ";

  switch( type ){
   case FFOp::VAR:
    if( !pres->cst() ) break; // do not override constant value if set
    //if( itU ) throw typename FFBase::Exceptions( FFBase::Exceptions::CONSTVAL );

   case FFOp::CNST:
    switch( pres->num().t ){
      case FFNum::INT:  *itU = pres->num().n; break;
      case FFNum::REAL: *itU = pres->num().x; break;
    }
    break;

   case FFOp::SHIFT:
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SHIFT w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = std::move( val += pops[1]->num().val() );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SHIFT w/o move for " << *pops[0] << std::endl;
#endif
      *itU = *static_cast<U*>( pops[0]->val() ) + pops[1]->num().val();
    }
    break;

   case FFOp::PLUS:
    if( &pops[0] == &pops[1] && pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SCALE w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = std::move( val *= 2. );
    }
    else if( &pops[0] == &pops[1] ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SCALE w/o move for " << *pops[0] << std::endl;
#endif
      *itU = *static_cast<U*>( pops[0]->val() ) * 2.;
    }
    else if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling PLUS w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = std::move( val += *static_cast<U*>( pops[1]->val() ) );
    }
    else if( pops[1]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling PLUS w/ move for " << *pops[1] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[1]->val() );
      *itU = std::move( val += *static_cast<U*>( pops[0]->val() ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling PLUS w/o move for " << *pops[0] << "," << *pops[1] << std::endl;
#endif
      *itU = *static_cast<U*>( pops[0]->val() ) + *static_cast<U*>( pops[1]->val() );
    }
    break;

   case FFOp::NEG:
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling NEG w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = - std::move( val );
    }
    else{
      *itU = - *static_cast<U*>( pops[0]->val() );
    }
    break;

   case FFOp::MINUS:
    if( &pops[0] == &pops[1] )
      *itU = 0.;
    else if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling MINUS w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = std::move( val -= *static_cast<U*>( pops[1]->val() ) );
    }
    //else if( pops[1]->mov() ){
    //  std::cout << "calling MINUS w/ move for " << *pops[1] << std::endl;
    //  U& val  = *static_cast<U*>( pops[1]->val() );
    //  U&& neg = operator-( std::move( val ) );
    //  *itU = std::move( neg += *static_cast<U*>( pops[0]->val() ) );
    //}
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling MINUS w/ move for " << *pops[0] << "," << *pops[1] << std::endl;
#endif
      *itU = *static_cast<U*>( pops[0]->val() ) - *static_cast<U*>( pops[1]->val() );
    }
    break;

   case FFOp::SCALE:
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SCALE w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = std::move( val *= pops[1]->num().val() );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SCALE w/o move for " << *pops[0] << std::endl;
#endif
      *itU = *static_cast<U*>( pops[0]->val() ) * pops[1]->num().val();
    }
    break;

   case FFOp::TIMES:
    if( &pops[0] == &pops[1] && pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SQR w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::sqr( std::move( val ) );
    }
    else if( &pops[0] == &pops[1] ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SQR w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::sqr( *static_cast<U*>( pops[0]->val() ) );
    }
    else if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling TIMES w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = std::move( val *= *static_cast<U*>( pops[1]->val() ) );
    }
    else if( pops[1]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling TIMES w/ move for " << *pops[1] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[1]->val() );
      *itU = std::move( val *= *static_cast<U*>( pops[0]->val() ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling TIMES w/o move for " << *pops[0] << "," << *pops[1] << std::endl;
#endif
      *itU = *static_cast<U*>( pops[0]->val() ) * *static_cast<U*>( pops[1]->val() );
    }
    break;

   case FFOp::INV:
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling INV w/ move for " << *pops[1] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[1]->val() );
      *itU = pops[0]->num().val() / std::move( val );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling INV w/o move for " << *pops[1] << std::endl;
#endif
      *itU = pops[0]->num().val() / *static_cast<U*>( pops[1]->val() );
    }
    break;

   case FFOp::DIV:  
    if( &pops[0] == &pops[1] )
      *itU = 1.;
    else if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling DIV w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = std::move( val /= *static_cast<U*>( pops[1]->val() ) );
    }
    //else if( pops[1]->mov() ){
    //  std::cout << "calling DIV w/ move for " << *pops[1] << std::endl;
    //  U& val  = *static_cast<U*>( pops[1]->val() );
    //  U&& inv = Op<U>::inv( std::move( val ) );
    //  *itU = std::move( inv *= *static_cast<U*>( pops[0]->val() ) );
    //}
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling DIV w/o move for " << *pops[0] << "," << *pops[1] << std::endl;
#endif
      *itU = *static_cast<U*>( pops[0]->val() ) / *static_cast<U*>( pops[1]->val() );
    }
    break;

   case FFOp::IPOW:
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling IPOW w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::pow( std::move( val ), pops[1]->num().n );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling IPOW w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::pow( *static_cast<U*>( pops[0]->val() ), pops[1]->num().n );
    }
    break;

   case FFOp::DPOW:
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling DPOW w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::pow( std::move( val ), pops[1]->num().x );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling DPOW w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::pow( *static_cast<U*>( pops[0]->val() ), pops[1]->num().x );
    }
    break;

   case FFOp::CHEB:
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling CHEB w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::cheb( std::move( val ), pops[1]->num().n );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling CHEB w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::cheb( *static_cast<U*>( pops[0]->val() ), pops[1]->num().n );
    }
    break;

   case FFOp::SQR:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SQR w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::sqr( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SQR w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::sqr( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::SQRT: 
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SQRT w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::sqrt( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SQRT w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::sqrt( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::EXP:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling EXP w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::exp( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling EXP w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::exp( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::LOG:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling LOG w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::log( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling LOG w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::log( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::XLOG:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling XLOG w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::xlog( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling XLOG w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::xlog( *static_cast<U*>( pops[0]->val() ) );
    }
    break; 

   case FFOp::COS:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling COS w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::cos( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling COS w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::cos( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::SIN:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SIN w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::sin( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SIN w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::sin( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::TAN:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling TAN w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::tan( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling TAN w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::tan( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::ACOS:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling ACOS w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::acos( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling ACOS w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::acos( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::ASIN:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling ASIN w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::asin( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling ASIN w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::asin( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::ATAN:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling ATAN w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::atan( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling ATAN w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::atan( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::COSH:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling COSH w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::cosh( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling COSH w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::cosh( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::SINH:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SINH w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::sinh( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling SINH w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::sinh( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::TANH:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling TANH w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::tanh( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling TANH w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::tanh( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::ERF:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling ERF w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::erf( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling ERF w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::erf( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::FABS:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling FABS w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::fabs( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling FABS w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::fabs( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::FSTEP:  
    if( pops[0]->mov() ){
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling FSTEP w/ move for " << *pops[0] << std::endl;
#endif
      U& val = *static_cast<U*>( pops[0]->val() );
      *itU = Op<U>::fstep( std::move( val ) );
    }
    else{
#ifdef MC__FFUNC_EVAL_MOVE
      std::cout << "calling FSTEP w/o move for " << *pops[0] << std::endl;
#endif
      *itU = Op<U>::fstep( *static_cast<U*>( pops[0]->val() ) );
    }
    break;

   case FFOp::MINF: 
    if( &pops[0] == &pops[1] && pops[0]->mov() )
      *itU = std::move( *static_cast<U*>( pops[0]->val() ) );
    else if( &pops[0] == &pops[1] )
      *itU = *static_cast<U*>( pops[0]->val() );
    else if( pops[0]->mov() )
      *itU = Op<U>::min( std::move( *static_cast<U*>( pops[0]->val() ) ),
                         *static_cast<U*>( pops[1]->val() ) );
    else if( pops[1]->mov() )
      *itU = Op<U>::min( std::move( *static_cast<U*>( pops[1]->val() ) ),
                         *static_cast<U*>( pops[0]->val() ) );
    else
      *itU = Op<U>::min( *static_cast<U*>( pops[0]->val() ),
                         *static_cast<U*>( pops[1]->val() ) );
    break;

   case FFOp::MAXF: 
    if( &pops[0] == &pops[1] && pops[0]->mov() )
      *itU = std::move( *static_cast<U*>( pops[0]->val() ) );
    else if( &pops[0] == &pops[1] )
      *itU = *static_cast<U*>( pops[0]->val() );
    else if( pops[0]->mov() )
      *itU = Op<U>::max( std::move( *static_cast<U*>( pops[0]->val() ) ),
                         *static_cast<U*>( pops[1]->val() ) );
    else if( pops[1]->mov() )
      *itU = Op<U>::max( std::move( *static_cast<U*>( pops[1]->val() ) ),
                         *static_cast<U*>( pops[0]->val() ) );
    else
      *itU = Op<U>::max( *static_cast<U*>( pops[0]->val() ),
                         *static_cast<U*>( pops[1]->val() ) );
    break;

   case FFOp::INTER: 
    if( &pops[0] == &pops[1] && pops[0]->mov() )
      *itU = std::move( *static_cast<U*>( pops[0]->val() ) );
    else if( &pops[0] == &pops[1] )
      *itU = *static_cast<U*>( pops[0]->val() );
    else if( !Op<U>::inter( *itU, *static_cast<U*>( pops[0]->val() ),
                                  *static_cast<U*>( pops[1]->val() ) ) )
      throw typename FFBase::Exceptions( FFBase::Exceptions::INTER );
    break;

   case FFOp::PROD:{
    std::vector<U> ops_val; ops_val.reserve( pops.size() );
    for( auto it=pops.begin(); it!=pops.end(); ++it ){
      if( (*it)->mov() )
        ops_val.push_back( *static_cast<U*>( std::move( (*it)->val() ) ) );
      else
        ops_val.push_back( *static_cast<U*>( (*it)->val() ) );
    }
    *itU = std::move( Op<U>::prod( ops_val.size(), ops_val.data() ) );
    break;
   }

   default:
    throw typename FFBase::Exceptions( FFBase::Exceptions::INTERN );
  }

  //pres->val() = &(*itU);
  //std::cout << "evaluation of " << *pres << ": " << &(*itU) << std::endl;
  //std::cout << &(*itU) << std::endl;
  return;
}

template <typename U, typename ExtOp, typename... NextOps>
inline void
FFOp::evaluate_external
( typename std::vector<U>::iterator itU, U const* pU_dum, bool const movU,
  ExtOp op, std::tuple<NextOps...> ops )
const
{
  if( type < FFOp::EXTERN )
    throw typename FFBase::Exceptions( FFBase::Exceptions::INTERN );

  // Current operation matches ExtOp type
  if( sameid( typeid(op) ) ){
    pres->val() = &(*itU);
    pres->mov() = movU;
    if( pops.empty() )
      throw typename FFBase::Exceptions( FFBase::Exceptions::EXTERN );
    else if( pops.size() == 1 ){
      op.eval( *itU, 1, static_cast<U*>( pops[0]->val() ) );
//      if( pops[0]->mov() )
//        *itU = op.eval( std::move( *static_cast<U*>( pops[0]->val() ) ) );
//      else
//        *itU = std::move( op.eval( *static_cast<U*>( pops[0]->val() ) ) );
//    }
//    else if( pops.size() == 2 ){
//      if( pops[0]->mov() )
//        *itU = op.eval( std::move( *static_cast<U*>( pops[0]->val() ) ), *static_cast<U*>( pops[1]->val() ) );
//      else if( pops[1]->mov() )
//        *itU = op.eval( std::move( *static_cast<U*>( pops[1]->val() ) ), *static_cast<U*>( pops[0]->val() ) );
//      else
//        *itU = std::move( op.eval( *static_cast<U*>( pops[0]->val() ), *static_cast<U*>( pops[0]->val() ) ) );
    }
    else{
      std::vector<U> ops_val; ops_val.reserve( pops.size() );
      for( auto const& pvar : pops ){
        if( pvar->mov() )
          ops_val.push_back( std::move( *static_cast<U*>( pvar->val() ) ) );
        else
          ops_val.push_back( *static_cast<U*>( pvar->val() ) );
      }
      op.eval( *itU, ops_val.size(), ops_val.data() );
    }
    return;
  }
  
  // No more external operations to peel off parameter pack
  if( !sizeof...( NextOps ) )
    throw typename FFBase::Exceptions( FFBase::Exceptions::EXTERN );

  // Separate first element off parameter pack and recursive call
  typedef std::tuple<NextOps...> t_NextOps;
  typedef typename remove_first_type< t_NextOps >::type t_NextNextOps;
  typedef typename first_type_of< NextOps... >::type FirstNextOps;
  evaluate_external( itU, pU_dum, movU, FirstNextOps(), t_NextNextOps() );
  return;
}

//template <typename U>
//inline U&&
//FFOp::eval
//( U&& Var )
//const
//{
//  return eval( std::forward( Var ) );
//}

//template <typename U>
//inline U
//FFOp::eval
//( U const& Var )
//const
//{
//  return eval( 1, &Var );
//}

//template <typename U>
//inline U&&
//FFOp::eval
//( U&& Var1, U const& Var2 )
//const
//{
//  return eval( std::forward( Var1 ), Var2 );
//}

//template <typename U>
//inline U
//FFOp::eval
//( U const& Var1, U const& Var2 )
//const
//{
//  std::vector<U> vVar; vVar.reserve( 2 );
//  ops_val.push_back( std::move( *static_cast<U*>( Var1.val() ) ) );
//  ops_val.push_back( *static_cast<U*>( Var2.val() ) );
//  return eval( 1, &vVar );
//}

template <typename U>
inline void
FFOp::eval
( U& vRes, unsigned const nVar, U const* vVar )
const
{
  throw typename FFBase::Exceptions( FFBase::Exceptions::EXTERN );
}

template <typename U> inline bool
FFOp::tighten_forward
( const U* pU_dum )
const
{
  switch( type ){
   case FFOp::VAR:
   case FFOp::CNST:
    break;

   case FFOp::SHIFT:
    //*itU = *static_cast<U*>( pops[0]->val() ) + pops[1]->num().val();
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       *static_cast<U*>( pops[0]->val() ) + pops[1]->num().val(),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::PLUS:
    if( &pops[0] == &pops[1] ){
      //*itU = *static_cast<U*>( pops[0]->val() ) * 2.;
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         *static_cast<U*>( pops[0]->val() ) * 2.,
                         *static_cast<U*>( pres->val() ) ) ) return false;
    }
    else{
      //*itU = *static_cast<U*>( pops[0]->val() ) + *static_cast<U*>( pops[1]->val() );
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         *static_cast<U*>( pops[0]->val() ) + *static_cast<U*>( pops[1]->val() ),
                         *static_cast<U*>( pres->val() ) ) ) return false;
    }
    break;

   case FFOp::NEG:
    //*itU = - *static_cast<U*>( pops[0]->val() );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       - *static_cast<U*>( pops[0]->val() ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::MINUS:
    if( &pops[0] == &pops[1] ){
      //*itU = 0.;
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         0.,
                         *static_cast<U*>( pres->val() ) ) ) return false;
    }
    else{
      //*itU = *static_cast<U*>( pops[0]->val() ) - *static_cast<U*>( pops[1]->val() );
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         *static_cast<U*>( pops[0]->val() ) - *static_cast<U*>( pops[1]->val() ),
                         *static_cast<U*>( pres->val() ) ) ) return false;
    }
    break;

   case FFOp::SCALE:
    //*itU = *static_cast<U*>( pops[0]->val() ) * pops[1]->num().val();
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       *static_cast<U*>( pops[0]->val() ) * pops[1]->num().val(),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::TIMES:
    if( &pops[0] == &pops[1] ){
      //*itU = Op<U>::sqr( *static_cast<U*>( pops[0]->val() ) );
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         Op<U>::sqr( *static_cast<U*>( pops[0]->val() ) ),
                         *static_cast<U*>( pres->val() ) ) ) return false;
    }
    else{
      //*itU = *static_cast<U*>( pops[0]->val() ) * *static_cast<U*>( pops[1]->val() );
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         *static_cast<U*>( pops[0]->val() ) * *static_cast<U*>( pops[1]->val() ),
                         *static_cast<U*>( pres->val() ) ) ) return false;
    }
    break;

   case FFOp::INV:
    //*itU = pops[0]->num().val() / *static_cast<U*>( pops[1]->val() );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       pops[0]->num().val() / *static_cast<U*>( pops[1]->val() ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::DIV:  
    if( &pops[0] == &pops[1] ){
      //*itU = 1.;
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         1.,
                         *static_cast<U*>( pres->val() ) ) ) return false;
    }
    else{
      //*itU = *static_cast<U*>( pops[0]->val() ) / *static_cast<U*>( pops[1]->val() );
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         *static_cast<U*>( pops[0]->val() ) / *static_cast<U*>( pops[1]->val() ),
                         *static_cast<U*>( pres->val() ) ) ) return false;
    }
    break;

   case FFOp::IPOW:
    //*itU = Op<U>::pow( *static_cast<U*>( pops[0]->val() ), pops[1]->num().n );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::pow( *static_cast<U*>( pops[0]->val() ), pops[1]->num().n ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::DPOW:
    //*itU = Op<U>::pow( *static_cast<U*>( pops[0]->val() ), pops[1]->num().x );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::pow( *static_cast<U*>( pops[0]->val() ), pops[1]->num().x ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::CHEB:
    //*itU = Op<U>::cheb( *static_cast<U*>( pops[0]->val() ), pops[1]->num().n );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::cheb( *static_cast<U*>( pops[0]->val() ), pops[1]->num().n ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::SQR:  
    //*itU = Op<U>::sqr( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::sqr( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::SQRT: 
    //*itU = Op<U>::sqrt( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::sqrt( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::EXP:  
    //*itU = Op<U>::exp( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::exp( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::LOG:  
    //*itU = Op<U>::log( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::log( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::XLOG:  
    //*itU = Op<U>::xlog( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::xlog( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break; 

   case FFOp::COS:  
    //*itU = Op<U>::cos( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::cos( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::SIN:  
    //*itU = Op<U>::sin( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::sin( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::TAN:  
    //*itU = Op<U>::tan( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::tan( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::ACOS: 
    //*itU = Op<U>::acos( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::acos( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::ASIN: 
    //*itU = Op<U>::asin( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::asin( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::ATAN: 
    //*itU = Op<U>::atan( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::atan( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::COSH: 
    //*itU = Op<U>::cosh( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::cosh( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::SINH: 
    //*itU = Op<U>::sinh( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::sinh( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::TANH: 
    //*itU = Op<U>::tanh( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::tanh( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::ERF:  
    //*itU = Op<U>::erf( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::erf( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::FABS: 
    //*itU = Op<U>::fabs( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::fabs( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::FSTEP:
    //*itU = Op<U>::fstep( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::fstep( *static_cast<U*>( pops[0]->val() ) ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;

   case FFOp::MINF: 
    if( &pops[0] == &pops[1] ){
      //*itU = *static_cast<U*>( pops[0]->val() );
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         *static_cast<U*>( pops[0]->val() ),
                         *static_cast<U*>( pres->val() ) ) ) return false;
    }
    else{
      //*itU = Op<U>::min( *static_cast<U*>( pops[0]->val() ), *static_cast<U*>( pops[1]->val() ) );
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         Op<U>::min( *static_cast<U*>( pops[0]->val() ),
                                     *static_cast<U*>( pops[1]->val() ) ),
                         *static_cast<U*>( pres->val() ) ) ) return false;
    }
    break;

   case FFOp::MAXF: 
    if( &pops[0] == &pops[1] ){
      //*itU = *static_cast<U*>( pops[0]->val() );
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         *static_cast<U*>( pops[0]->val() ),
                         *static_cast<U*>( pres->val() ) ) ) return false;
    }
    else{
      //*itU = Op<U>::max( *static_cast<U*>( pops[0]->val() ), *static_cast<U*>( pops[1]->val() ) );
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         Op<U>::max( *static_cast<U*>( pops[0]->val() ),
                                     *static_cast<U*>( pops[1]->val() ) ),
                         *static_cast<U*>( pres->val() ) ) ) return false;
    }
    break;

   case FFOp::INTER: 
    if( &pops[0] == &pops[1] ){
      //*itU = *static_cast<U*>( pops[0]->val() );
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         *static_cast<U*>( pops[0]->val() ),
                         *static_cast<U*>( pres->val() ) ) ) return false;
    }
    else{
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         *static_cast<U*>( pops[0]->val() ),
                         *static_cast<U*>( pres->val() ) ) ) return false;
      if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                         *static_cast<U*>( pops[1]->val() ),
                         *static_cast<U*>( pres->val() ) ) ) return false;
    }
    break;

   case FFOp::PROD:{
    std::vector<U> ops_val; ops_val.reserve( pops.size() );
    for( auto it=pops.begin(); it!=pops.end(); ++it )
      ops_val.push_back( *static_cast<U*>( (*it)->val() ) );
    //*itU = Op<U>::prod( ops_val.size(), ops_val.data() );
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ),
                       Op<U>::prod( ops_val.size(), ops_val.data() ),
                       *static_cast<U*>( pres->val() ) ) ) return false;
    break;
   }

   default:
    throw typename FFBase::Exceptions( FFBase::Exceptions::INTERN );
  }

  //pres->val() = &(*itU);
  //std::cout << "evaluation of " << *pres << ":\n" << *itU;
  return true;
}

template <typename U, typename ExtOp, typename... NextOps>
inline bool
FFOp::tighten_forward_external
( U const* pU_dum, ExtOp op, std::tuple<NextOps...> ops )
const
{
  if( type < FFOp::EXTERN )
    throw typename FFBase::Exceptions( FFBase::Exceptions::INTERN );

  // Current operation matches ExtOp type
  if( sameid( typeid(op) ) ){
    U vres;
    if( pops.empty() )
      throw typename FFBase::Exceptions( FFBase::Exceptions::EXTERN );
    else if( pops.size() == 1 )
      op.eval( vres, 1, static_cast<U*>( pops[0]->val() ) );
     //else if( pops.size() == 2 )
    //  op.eval( vres, static_cast<U*>( pops[0]->val() ), static_cast<U*>( pops[1]->val() ), false );
    else{
      std::vector<U> ops_val; ops_val.reserve( pops.size() );
      for( auto it=pops.begin(); it!=pops.end(); ++it )
        ops_val.push_back( *static_cast<U*>( (*it)->val() ) );
      op.eval( vres, ops_val.size(), ops_val.data() );
    }
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ), vres,
                       *static_cast<U*>( pres->val() ) ) ) return false;
    return true;
  }
  
  // No more external operations to peel off parameter pack
  if( !sizeof...( NextOps ) )
    throw typename FFBase::Exceptions( FFBase::Exceptions::EXTERN );

  // Separate first element off parameter pack and recursive call
  typedef std::tuple<NextOps...> t_NextOps;
  typedef typename remove_first_type< t_NextOps >::type t_NextNextOps;
  typedef typename first_type_of< NextOps... >::type FirstNextOps;
  return tighten_forward_external( pU_dum, FirstNextOps(), t_NextNextOps() );
}   

template <typename U> inline bool
FFOp::tighten_backward
( const U* pU_dum ) const
{
  switch( type ){
   case FFOp::VAR:
   case FFOp::CNST:
    break;

   case FFOp::SHIFT:
    //*itU = *static_cast<U*>( pops[0]->val() ) + pops[1]->num().val();
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       *static_cast<U*>( pres->val() ) - pops[1]->num().val(),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;
    break;

   case FFOp::PLUS:
    if( &pops[0] == &pops[1] ){
      //*itU = *static_cast<U*>( pops[0]->val() ) * 2.;
      if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                         *static_cast<U*>( pres->val() ) / 2.,
                         *static_cast<U*>( pops[0]->val() ) ) ) return false;
    }
    else{
      //*itU = *static_cast<U*>( pops[0]->val() ) + *static_cast<U*>( pops[1]->val() );
      if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                         *static_cast<U*>( pres->val() ) - *static_cast<U*>( pops[1]->val() ),
                         *static_cast<U*>( pops[0]->val() ) ) ) return false;
      if( !Op<U>::inter( *static_cast<U*>( pops[1]->val() ),
                         *static_cast<U*>( pres->val() ) - *static_cast<U*>( pops[0]->val() ),
                         *static_cast<U*>( pops[1]->val() ) ) ) return false;
    }
    break;

   case FFOp::NEG:
    //*itU = - *static_cast<U*>( pops[0]->val() );
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       - *static_cast<U*>( pres->val() ),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;
    return true;

   case FFOp::MINUS:
    if( &pops[0] == &pops[1] ) break; // nothing to tighten
    //*itU = *static_cast<U*>( pops[0]->val() ) - *static_cast<U*>( pops[1]->val() );
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       *static_cast<U*>( pres->val() ) + *static_cast<U*>( pops[1]->val() ),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;
    if( !Op<U>::inter( *static_cast<U*>( pops[1]->val() ),
                       *static_cast<U*>( pops[0]->val() ) - *static_cast<U*>( pres->val() ),
                       *static_cast<U*>( pops[1]->val() ) ) ) return false;
    break;

   case FFOp::SCALE:
    //*itU = *static_cast<U*>( pops[0]->val() ) * pops[1]->num().val();
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       *static_cast<U*>( pres->val() ) / pops[1]->num().val(),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;
    break;

   case FFOp::TIMES:
    if( &pops[0] == &pops[1] ){
      //*itU = Op<U>::sqr( *static_cast<U*>( pops[0]->val() ) );
      if( Op<U>::l( *static_cast<U*>( pops[0]->val() ) ) >= 0. ){
        if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                           Op<U>::sqrt( *static_cast<U*>( pres->val() ) ),
                           *static_cast<U*>( pops[0]->val() ) ) ) return false;
      }
      else if( Op<U>::u( *static_cast<U*>( pops[0]->val() ) ) <= 0. ){
        if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                           - Op<U>::sqrt( *static_cast<U*>( pres->val() ) ),
                           *static_cast<U*>( pops[0]->val() ) ) ) return false;
      }
      else{
        if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                           Op<U>::hull( - Op<U>::sqrt( *static_cast<U*>( pres->val() ) ),
                                          Op<U>::sqrt( *static_cast<U*>( pres->val() ) ) ),
                           *static_cast<U*>( pops[0]->val() ) ) ) return false;
      }
    }
    else{
      if( *static_cast<U*>( pops[0]->val() ) == U(0)
       || *static_cast<U*>( pops[1]->val() ) == U(0) ) break;
      //*itU = *static_cast<U*>( pops[0]->val() ) * *static_cast<U*>( pops[1]->val() );
      if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                         *static_cast<U*>( pres->val() ) / *static_cast<U*>( pops[1]->val() ),
                         *static_cast<U*>( pops[0]->val() ) ) ) return false;
      if( !Op<U>::inter( *static_cast<U*>( pops[1]->val() ),
                         *static_cast<U*>( pres->val() ) / *static_cast<U*>( pops[0]->val() ),
                         *static_cast<U*>( pops[1]->val() ) ) ) return false;
    }
    break;

   case FFOp::INV:
    //*itU = pops[0]->num().val() / *static_cast<U*>( pops[1]->val() );
    // should not reach this point since forward propagation would have normally thrown an exception earlier
    if( Op<U>::l( *static_cast<U*>( pops[0]->val() ) ) < 0. 
     && Op<U>::u( *static_cast<U*>( pops[0]->val() ) ) > 0. )
      throw typename FFBase::Exceptions( FFBase::Exceptions::EVAL );
    if( !Op<U>::inter( *static_cast<U*>( pops[1]->val() ),
                       pops[0]->num().val() / *static_cast<U*>( pres->val() ),
                       *static_cast<U*>( pops[1]->val() ) ) ) return false;
    break;

   case FFOp::DIV:  
    if( &pops[0] == &pops[1] ) break; // nothing to tighten
    //*itU = *static_cast<U*>( pops[0]->val() ) / *static_cast<U*>( pops[1]->val() );
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       *static_cast<U*>( pres->val() ) * *static_cast<U*>( pops[1]->val() ),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;
    if( !Op<U>::inter( *static_cast<U*>( pops[1]->val() ),
                       *static_cast<U*>( pops[0]->val() ) / *static_cast<U*>( pres->val() ),
                       *static_cast<U*>( pops[1]->val() ) ) ) return false;
    break;

   case FFOp::IPOW:
    //*itU = Op<U>::pow( *static_cast<U*>( pops[0]->val() ), pops[1]->num().n );
    if( pops[1]->num().n > 0 && pops[1]->num().n % 2 ){ // positive odd exponent
      if( Op<U>::l( *static_cast<U*>( pres->val() ) ) >= 0. ){
        if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                           Op<U>::pow( *static_cast<U*>( pres->val() ), 1./pops[1]->num().n ),
                           *static_cast<U*>( pops[0]->val() ) ) ) return false;
      }
      else if( Op<U>::u( *static_cast<U*>( pres->val() ) ) <= 0. ){
        if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                           - Op<U>::pow( - *static_cast<U*>( pres->val() ), 1./pops[1]->num().n ),
                           *static_cast<U*>( pops[0]->val() ) ) ) return false;
      }
      else{
        U bndL = Op<U>::zeroone() * Op<U>::l( *static_cast<U*>( pres->val() ) ),
          bndU = Op<U>::zeroone() * Op<U>::u( *static_cast<U*>( pres->val() ) );
        if( !Op<U>::inter( bndL, *static_cast<U*>( pres->val() ), bndL )
         || !Op<U>::inter( bndU, *static_cast<U*>( pres->val() ), bndU ) )
          throw typename FFBase::Exceptions( FFBase::Exceptions::EVAL );
        if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                           Op<U>::hull( - Op<U>::pow( - bndL, 1./pops[1]->num().n ),
                                          Op<U>::pow(   bndU, 1./pops[1]->num().n ) ),
                           *static_cast<U*>( pops[0]->val() ) ) ) return false;
      }
    }
    else if( pops[1]->num().n > 0 ){ // positive even exponent
      if( Op<U>::l( *static_cast<U*>( pops[0]->val() ) ) >= 0. ){
        if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                           Op<U>::pow( *static_cast<U*>( pres->val() ), 1./pops[1]->num().n ),
                           *static_cast<U*>( pops[0]->val() ) ) ) return false;
      }
      else if( Op<U>::u( *static_cast<U*>( pops[0]->val() ) ) <= 0. ){
        if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                           - Op<U>::pow( *static_cast<U*>( pres->val() ), 1./pops[1]->num().n ),
                           *static_cast<U*>( pops[0]->val() ) ) ) return false;
      }
      else{
        if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                           Op<U>::hull( - Op<U>::pow( *static_cast<U*>( pres->val() ), 1./pops[1]->num().n ),
                                          Op<U>::pow( *static_cast<U*>( pres->val() ), 1./pops[1]->num().n ) ),
                           *static_cast<U*>( pops[0]->val() ) ) ) return false;
      }
    }
    else if( pops[1]->num().n % 2 ){ // negative odd exponent
      if( Op<U>::l( *static_cast<U*>( pops[0]->val() ) ) > 0. ){
        if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                           Op<U>::pow( *static_cast<U*>( pres->val() ), 1./pops[1]->num().n ),
                           *static_cast<U*>( pops[0]->val() ) ) ) return false;
      }
      else if( Op<U>::u( *static_cast<U*>( pops[0]->val() ) ) < 0. ){
        if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                           - Op<U>::pow( - *static_cast<U*>( pres->val() ), 1./pops[1]->num().n ),
                           *static_cast<U*>( pops[0]->val() ) ) ) return false;
      }
      // should not reach this point since forward propagation would have normally thrown an exception earlier
      throw typename FFBase::Exceptions( FFBase::Exceptions::EVAL );
    }
    else{ // negative even exponent
      if( Op<U>::l( *static_cast<U*>( pops[0]->val() ) ) > 0. ){
        if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                           Op<U>::pow( *static_cast<U*>( pres->val() ), 1./pops[1]->num().n ),
                           *static_cast<U*>( pops[0]->val() ) ) ) return false;
      }
      else if( Op<U>::u( *static_cast<U*>( pops[0]->val() ) ) < 0. ){
        if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                           - Op<U>::pow( *static_cast<U*>( pres->val() ), 1./pops[1]->num().n ),
                           *static_cast<U*>( pops[0]->val() ) ) ) return false;
      }
      // should not reach this point since forward propagation would have normally thrown an exception earlier
      throw typename FFBase::Exceptions( FFBase::Exceptions::EVAL );
    }
    break;

   case FFOp::DPOW:
    //*itU = Op<U>::pow( *static_cast<U*>( pops[0]->val() ), pops[1]->num().x );
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       Op<U>::pow( *static_cast<U*>( pres->val() ), 1./pops[1]->num().x ),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;
    break;

   case FFOp::CHEB:
    //*itU = Op<U>::cheb( *static_cast<U*>( pops[0]->val() ), pops[1]->num().n );
    // TBC
    break;

   case FFOp::SQR:
    //*itU = Op<U>::sqr( *static_cast<U*>( pops[0]->val() ) );
    if( Op<U>::l( *static_cast<U*>( pops[0]->val() ) ) >= 0. ){
      if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                         Op<U>::sqrt( *static_cast<U*>( pres->val() ) ),
                         *static_cast<U*>( pops[0]->val() ) ) ) return false;
    }
    else if( Op<U>::u( *static_cast<U*>( pops[0]->val() ) ) <= 0. ){
      if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                         - Op<U>::sqrt( *static_cast<U*>( pres->val() ) ),
                         *static_cast<U*>( pops[0]->val() ) ) ) return false;
    }
    else{
      if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                         Op<U>::hull( - Op<U>::sqrt( *static_cast<U*>( pres->val() ) ),
                                        Op<U>::sqrt( *static_cast<U*>( pres->val() ) ) ),
                         *static_cast<U*>( pops[0]->val() ) ) ) return false;
    }
    break;

   case FFOp::SQRT: 
    //*itU = Op<U>::sqrt( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       Op<U>::sqr( *static_cast<U*>( pres->val() ) ),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;
    break;

   case FFOp::EXP:  
    //*itU = Op<U>::exp( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       Op<U>::log( *static_cast<U*>( pres->val() ) ),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;
    break;

   case FFOp::LOG:  
    //*itU = Op<U>::log( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       Op<U>::exp( *static_cast<U*>( pres->val() ) ),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;
    break;

   case FFOp::XLOG:  
    //*itU = Op<U>::xlog( *static_cast<U*>( pops[0]->val() ) );
    // TBC
    break; 

   case FFOp::COS:  
    //*itU = Op<U>::cos( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       Op<U>::acos( *static_cast<U*>( pres->val() ) ),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;
    break;

   case FFOp::SIN:  
    //*itU = Op<U>::sin( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       Op<U>::asin( *static_cast<U*>( pres->val() ) ),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;
    break;

   case FFOp::TAN:  
    //*itU = Op<U>::tan( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       Op<U>::atan( *static_cast<U*>( pres->val() ) ),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;
    break;

   case FFOp::ACOS: 
    //*itU = Op<U>::acos( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       Op<U>::cos( *static_cast<U*>( pres->val() ) ),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;
    break;

   case FFOp::ASIN: 
    //*itU = Op<U>::asin( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       Op<U>::sin( *static_cast<U*>( pres->val() ) ),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;
    break;

   case FFOp::ATAN: 
    //*itU = Op<U>::atan( *static_cast<U*>( pops[0]->val() ) );
    if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                       Op<U>::tan( *static_cast<U*>( pres->val() ) ),
                       *static_cast<U*>( pops[0]->val() ) ) ) return false;

    break;

   case FFOp::COSH: 
    //*itU = Op<U>::cosh( *static_cast<U*>( pops[0]->val() ) );
    // TBC
    break;

   case FFOp::SINH: 
    //*itU = Op<U>::sinh( *static_cast<U*>( pops[0]->val() ) );
    // TBC
    break;

   case FFOp::TANH: 
    //*itU = Op<U>::tanh( *static_cast<U*>( pops[0]->val() ) );
    // TBC
    break;

   case FFOp::ERF:  
    //*itU = Op<U>::erf( *static_cast<U*>( pops[0]->val() ) );
    // TBC
    break;

   case FFOp::FABS: 
    //*itU = Op<U>::fabs( *static_cast<U*>( pops[0]->val() ) );
    if( Op<U>::l( *static_cast<U*>( pops[0]->val() ) ) >= 0. ){
      if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                         *static_cast<U*>( pres->val() ),
                         *static_cast<U*>( pops[0]->val() ) ) ) return false;
    }
    else if( Op<U>::u( *static_cast<U*>( pops[0]->val() ) ) <= 0. ){
      if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                         - *static_cast<U*>( pres->val() ),
                         *static_cast<U*>( pops[0]->val() ) ) ) return false;
    }
    else{
      if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                         Op<U>::hull( - *static_cast<U*>( pres->val() ),
                                        *static_cast<U*>( pres->val() ) ),
                         *static_cast<U*>( pops[0]->val() ) ) ) return false;
    }
    break;

   case FFOp::FSTEP:
    //*itU = Op<U>::fstep( *static_cast<U*>( pops[0]->val() ) );
    if( Op<U>::l( *static_cast<U*>( pops[0]->val() ) ) >= 0. ){
      if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                         1.,
                         *static_cast<U*>( pops[0]->val() ) ) ) return false;
    }
    else if( Op<U>::u( *static_cast<U*>( pops[0]->val() ) ) <= 0. ){
      if( !Op<U>::inter( *static_cast<U*>( pops[0]->val() ),
                         0.,
                         *static_cast<U*>( pops[0]->val() ) ) ) return false;
    }
    break;

   case FFOp::MINF: 
    //if( &pops[0] == &pops[1] )
    //  *itU = *static_cast<U*>( pops[0]->val() );
    //else
    //  *itU = Op<U>::min( *static_cast<U*>( pops[0]->val() ), *static_cast<U*>( pops[1]->val() ) );
    // TBC
    break;

   case FFOp::MAXF: 
    //if( &pops[0] == &pops[1] )
    //  *itU = *static_cast<U*>( pops[0]->val() );
    //else
    //  *itU = Op<U>::max( *static_cast<U*>( pops[0]->val() ), *static_cast<U*>( pops[1]->val() ) );
    // TBC
    break;

   case FFOp::INTER: 
    //if( &pops[0] == &pops[1] )
    //  *itU = *static_cast<U*>( pops[0]->val() );
    //else if( !Op<U>::inter( *itU, *static_cast<U*>( pops[0]->val() ), *static_cast<U*>( pops[1]->val() ) ) )
    //  throw typename FFBase::Exceptions( FFBase::Exceptions::INTER );
    // TBC
    break;

   case FFOp::PROD:{
    //std::vector<U> ops_val; ops_val.reserve( pops.size() );
    //for( auto it=pops.begin(); it!=pops.end(); ++it ) ops_val.push_back( *static_cast<U*>( (*it)->val() ) );
    //*itU = Op<U>::prod( ops_val.size(), ops_val.data() );
    // TBC
    break;
   }

   default:
    throw typename FFBase::Exceptions( FFBase::Exceptions::INTERN );
  }

  //pres->val() = &(*itU);
  //std::cout << "evaluation of " << *pres << ":\n" << *itU;
  return true;
}

template <typename U, typename ExtOp, typename... NextOps>
inline bool
FFOp::tighten_backward_external
( U const* pU_dum, ExtOp op, std::tuple<NextOps...> ops )
const
{
  if( type < FFOp::EXTERN )
    throw typename FFBase::Exceptions( FFBase::Exceptions::INTERN );

  // Current operation matches ExtOp type
  if( sameid( typeid(op) ) ){
    U vres;
    if( pops.empty() )
      throw typename FFBase::Exceptions( FFBase::Exceptions::EXTERN );
    else if( pops.size() == 1 ){
      //op.eval( vres, 1, static_cast<U*>( pops[0]->val() ) );
      if( !op.reval( vres, 1, static_cast<U*>( pops[0]->val() ) ) ) return false;
    }
    //else if( pops.size() == 2 )
    //  op.eval( *itU, *static_cast<U*>( pops[0]->val() ), *static_cast<U*>( pops[1]->val() ) );
    //  if( !op.reval( vres, static_cast<U*>( pops[0]->val() ),
    //                       static_cast<U*>( pops[1]->val() ) ) ) return false;
    else{
      std::vector<U> ops_val; ops_val.reserve( pops.size() );
      for( auto it=pops.begin(); it!=pops.end(); ++it ) ops_val.push_back( *static_cast<U*>( (*it)->val() ) );
      //op.eval( vres, ops_val.size(), ops_val.data() );
      if( !op.reval( vres, ops_val.size(), ops_val.data() ) ) return false;
    }
    return true;
  }
  
  // No more external operations to peel off parameter pack
  if( !sizeof...( NextOps ) )
    throw typename FFBase::Exceptions( FFBase::Exceptions::EXTERN );

  // Separate first element off parameter pack and recursive call
  typedef std::tuple<NextOps...> t_NextOps;
  typedef typename remove_first_type< t_NextOps >::type t_NextNextOps;
  typedef typename first_type_of< NextOps... >::type FirstNextOps;
  return tighten_forward_external( pU_dum, FirstNextOps(), t_NextNextOps() );
}   

template <typename U>
inline bool
FFOp::reval
( U const& vRes, unsigned const nVar, U* vVar )
const
{
  return true;
  //throw typename FFBase::Exceptions( FFBase::Exceptions::EXTERN );
}

inline void
FFOp::generate_dot_script
( std::ostream&os ) const
{
  if( iflag ) return;
  iflag = 1;

  for( auto it=pops.begin(); it!=pops.end(); ++it ){
    if( !(*it) || !(*it)->ops().first ) continue;
    (*it)->ops().first->generate_dot_script( os );
  }
  append_dot_script( os );
}

inline void
FFOp::append_dot_script
( std::ostream&os ) const
{
  switch( type ){
   case FFOp::VAR:   
   case FFOp::CNST:   return append_dot_script_variable( os, 14 );

   case FFOp::SHIFT:
   case FFOp::PLUS:  
   case FFOp::NEG:   
   case FFOp::MINUS: 
   case FFOp::SCALE:
   case FFOp::TIMES: 
   case FFOp::INV:
   case FFOp::DIV:
   case FFOp::INTER:
   case FFOp::PROD:
   case FFOp::IPOW:
   case FFOp::DPOW:
   case FFOp::CHEB:
   case FFOp::SQR:
   case FFOp::SQRT:
   case FFOp::EXP:
   case FFOp::LOG:
   case FFOp::XLOG:
   case FFOp::COS:
   case FFOp::SIN:
   case FFOp::TAN:
   case FFOp::ACOS:
   case FFOp::ASIN:
   case FFOp::ATAN:
   case FFOp::COSH:
   case FFOp::SINH:
   case FFOp::TANH:
   case FFOp::ERF:
   case FFOp::FABS:
   case FFOp::FSTEP:
   case FFOp::MINF:
   case FFOp::MAXF:  return append_dot_script_factor( os, 14 );

   default: if( type >= FFOp::EXTERN ) return append_dot_script_factor( os, 14 );
            // Should not reach this point
            throw typename FFBase::Exceptions( FFBase::Exceptions::INTERN );
            os << "/* a factor was not displayed */\n";
  }
}

inline void
FFOp::append_dot_script_factor
( std::ostream&os, const unsigned int fontsize )
const
{
  std::ostringstream op_color; op_color  << "black";

  os << "  " << pres->name() << " [shape=Mrecord,fontname=\"Arial\",color=" << op_color.str().c_str()
     << ",fontsize=" << fontsize << ",label=\"{<f0> " << name().c_str() << "|<f1> " << pres->name() << "}\"];\n";
  unsigned int i=1;
  for( auto it=pops.begin(); it!=pops.end(); ++it, ++i ){
    //std::ostringstream label; label  << ",labeldistance=0.8,labelangle=0,taillabel=<<table bgcolor=\"white\" border=\"0\"><tr><td>" << i << "</td></tr></table>>,fontname=\"Arial\",fontsize=10";
    std::ostringstream label; label  << ",label=\"" << i << "\",fontname=\"Arial\",fontsize=10";
    os << "  " << (*it)->name() << " -> " << pres->name() << " [arrowsize=0.7"
       << ((pops.size()>1 && !commutative())? label.str().c_str(): "") << "];\n";
       //<< (!commutative() && (it!=pops.begin())? ",style=dashed];\n": "];\n");
  }
}

inline void
FFOp::append_dot_script_variable
( std::ostream&os, const unsigned int fontsize ) const
{
  std::ostringstream var_color; var_color << "red";
  std::ostringstream cst_color; cst_color << "blue";

  if( type == FFOp::CNST || pres->cst() )
    os << "  " << pres->name() << " [shape=Mrecord,fontname=\"Arial\",color=" << cst_color.str().c_str()
       << ",fontsize=" << fontsize << ",label=\"{<f0> " << pres->num() << "|<f1> " << pres->name() << "}\"];\n"; 
  else
    os << "  " << pres->name() << " [shape=Mrecord,fontname=\"Arial\",color=" << (var_color.str().c_str()) << "];\n";
    //   << (var_color.str().c_str()) << ",label=\"<f0> " << pres->name() << "\"];\n";
}

inline bool
FFOp::commutative
() const
{
  switch( type ){
   // SHIFT and SCALE operators not considered commutative 
   case FFOp::PLUS:
   case FFOp::TIMES:
   case FFOp::MINF:
   case FFOp::MAXF:
   case FFOp::INTER:
   case FFOp::PROD:   return true;
   default:           return false;
  }
}

inline std::string
FFOp::name
()
const
{
  switch( type ){
   case FFOp::VAR:    return "VAR";;
   case FFOp::CNST:   return "CONST";
   case FFOp::SHIFT:
   case FFOp::PLUS:   return " + ";
   case FFOp::NEG:   
   case FFOp::MINUS:  return " - ";
   case FFOp::SCALE:
   case FFOp::TIMES:  return " x ";
   case FFOp::INV:
   case FFOp::DIV:    return " / ";
   case FFOp::MINF:   return "MIN";
   case FFOp::MAXF:   return "MAX";
   case FFOp::INTER:  return "INTER";
   case FFOp::PROD:   return "PROD";
   case FFOp::IPOW:   return "IPOW";
   case FFOp::DPOW:   return "DPOW";
   case FFOp::CHEB:   return "CHEB";
   case FFOp::SQR:    return "SQR";
   case FFOp::SQRT:   return "SQRT";
   case FFOp::EXP:    return "EXP";
   case FFOp::LOG:    return "LOG";
   case FFOp::XLOG:   return "XLOG";
   case FFOp::FABS:   return "FABS";
   case FFOp::COS:    return "COS";
   case FFOp::SIN:    return "SIN";
   case FFOp::TAN:    return "TAN";
   case FFOp::ACOS:   return "ACOS";
   case FFOp::ASIN:   return "ASIN";
   case FFOp::ATAN:   return "ATAN";
   case FFOp::COSH:   return "COSH";
   case FFOp::SINH:   return "SINH";
   case FFOp::TANH:   return "TANH";
   case FFOp::ERF:    return "ERF";
   case FFOp::FSTEP:  return "FSTEP";

   default: if( type >= FFOp::EXTERN ) return "EXTERN";
            // Should not reach this point
            throw typename FFBase::Exceptions( FFBase::Exceptions::INTERN );
  }
}

inline bool
FFOp::sameid
( std::type_info const& id )
const
{
  return( typeid(*this) == id );
}

///////////////////////////////// FFBase //////////////////////////////////////

inline std::ostream&
operator <<
( std::ostream& out, FFBase const& dag)
{
  typename FFBase::t_Vars Vars = dag._Vars;
  typename FFBase::it_Vars itv = Vars.begin();

  out << ( dag._nvar? "\nDAG VARIABLES:\n": "\nNO DAG VARIABLES\n" );
  for( ; itv!=Vars.end() && (*itv)->_id.first<=FFVar::VAR; ++itv ){
    //out << "  " << **itv << "  (" << *itv << ")";
    out << "  " << **itv;
    out << "\t => {";
    typename FFVar::t_Ops Ops = (*itv)->_ops.second;
    typename FFVar::t_Ops::iterator ito = Ops.begin();
    for( ; ito!=Ops.end(); ++ito ) out << " " << *(*ito)->pres;
    out << " }" << std::endl;
  }

  out << ( dag._naux? "\nDAG INTERMEDIATES:\n": "\nNO DAG INTERMEDIATES\n" );
  for( ; itv!=Vars.end(); ++itv ){
    //out << "  " << **itv << "  (" << *itv << ")";
    out << "  " << **itv;
    if( (*itv)->_ops.first ) out << "\t" << "<=  " << *((*itv)->_ops.first);
    out << "\t => {";
    typename FFVar::t_Ops Ops = (*itv)->_ops.second;
    typename FFVar::t_Ops::iterator ito = Ops.begin();
    for( ; ito!=Ops.end(); ++ito ) out << " " << *(*ito)->pres;
    out << " }" << std::endl;
  }

  return out;
}

template <typename ExtOp>
inline FFVar&
FFBase::_insert_nary_external_operation
( ExtOp const& Op, FFDep const& dep, unsigned const nVar, FFVar const* pVar )
{
  auto dag = pVar[0]._dag;
  for( unsigned i=1; i<nVar && dag==nullptr; i++ )
    if( pVar[i]._dag != nullptr ) dag = pVar[i]._dag;
  if( dag == nullptr ) throw Exceptions( Exceptions::DAG );

  std::vector<FFVar*> vVar; vVar.reserve( nVar );
  for( unsigned i=0; i<nVar; i++ ){
    if( pVar[i]._dag == nullptr && pVar[i]._cst ){
      FFVar* pCst = dag->_add_constant( pVar[i]._num.val() );
      vVar.push_back( pCst->_ops.first->pres );
    }
    else
      vVar.push_back( pVar[i]._ops.first->pres );
  }
  
  FFOp* pOp = new ExtOp();
  pOp->set( nVar, vVar.data(), nullptr );
  auto itOp = dag->_Ops.find( pOp );
  if( itOp != dag->_Ops.end() ){
    delete pOp;
    pOp = *itOp;
  }
  dag->_Ops.insert( pOp );

  if( pOp->pres != nullptr ) return *pOp->pres;
  for( unsigned i=0; i<nVar; i++ ) vVar[i]->_ops.second.push_back( pOp );
  pOp->pres = dag->_add_auxiliary( dep, pOp );
  return *(pOp->pres);
}

template <typename ExtOp>
inline FFVar&
FFBase::_insert_binary_external_operation
( ExtOp const& Op, FFDep const& dep, FFVar const& Var1, FFVar const& Var2 )
{
  if( Var1._dag == nullptr && Var1._dag == nullptr ) throw Exceptions( Exceptions::DAG );
  auto dag = ( Var1._dag == nullptr? Var2._dag: Var1._dag );

  FFVar *pVar1 = nullptr;
  if( Var1._dag == nullptr && Var1._cst ){
    FFVar* pCst1 = dag->_add_constant( Var1._num.val() );
    pVar1 = pCst1->_ops.first->pres;
  }
  else
    pVar1 = Var1._ops.first->pres;  

  FFVar *pVar2 = nullptr;  
  if( Var2._dag == nullptr && Var2._cst ){
    FFVar* pCst2 = dag->_add_constant( Var2._num.val() );
    pVar2 = pCst2->_ops.first->pres;
  }
  else
    pVar2 = Var2._ops.first->pres;  

  FFOp* pOp = new ExtOp();
  pOp->set( pVar1, pVar2, nullptr );
  auto itOp = dag->_Ops.find( pOp );
  if( itOp != dag->_Ops.end() ){
    delete pOp;
    pOp = *itOp;
  }
  dag->_Ops.insert( pOp );

  if( pOp->pres ) return *pOp->pres;
  pVar1->_ops.second.push_back( pOp );
  pVar2->_ops.second.push_back( pOp );
  pOp->pres = dag->_add_auxiliary( dep, pOp );
  return *(pOp->pres);
}

template <typename ExtOp>
inline FFVar&
FFBase::_insert_unary_external_operation
( ExtOp const& Op, FFDep const& dep, FFVar const& Var )
{
  if( Var._dag == nullptr ) throw Exceptions( Exceptions::DAG );
  FFVar *pVar = Var._ops.first->pres;

  FFOp* pOp = new ExtOp();
  pOp->set( pVar, nullptr );
  auto itOp = Var._dag->_Ops.find( pOp );
  if( itOp != Var._dag->_Ops.end() ){
    delete pOp;
    pOp = *itOp;
  }
  Var._dag->_Ops.insert( pOp );

  if( pOp->pres ) return *pOp->pres;
  pVar->_ops.second.push_back( pOp );
  pOp->pres = Var._dag->_add_auxiliary( dep, pOp );
  return *(pOp->pres);
}

inline FFVar&
FFBase::_insert_nary_operation
( int const top, FFDep const& dep, unsigned const nVar, FFVar const* pVar )
{
  for( unsigned i=1; i<nVar; i++ ) if( pVar[0]._dag != pVar[i]._dag ) throw Exceptions( Exceptions::DAG );
  FFBase *pdag = pVar[0]._dag;
  std::vector<FFVar*> vVar;
  for( unsigned i=0; i<nVar; i++ ) vVar.push_back( pVar[i]._ops.first->pres );

  FFOp *pOp = pdag->_insert_operation( top, nVar, vVar.data() );
  if( pOp->pres ) return *pOp->pres;
  for( unsigned i=0; i<nVar; i++ ) vVar[i]->_ops.second.push_back( pOp );
  pOp->pres = pdag->_add_auxiliary( dep, pOp );
  return *pOp->pres;
}

inline FFVar&
FFBase::_insert_binary_operation
( int const top, FFDep const& dep, FFVar const& Var1, FFVar const& Var2 )
{
  if( Var1._dag != Var2._dag ) throw Exceptions( Exceptions::DAG );
  FFBase *pdag = Var1._dag;
  FFVar *pVar1 = Var1._ops.first->pres, *pVar2 = Var2._ops.first->pres;

  FFOp *pOp = pdag->_insert_operation( top, pVar1, pVar2 );
  if( pOp->pres ) return *pOp->pres;
  pVar1->_ops.second.push_back( pOp );
  pVar2->_ops.second.push_back( pOp );
  pOp->pres = pdag->_add_auxiliary( dep, pOp );
  return *pOp->pres;
}

template <typename U> inline FFVar&
FFBase::_insert_binary_operation
( int const top, FFDep const& dep, U const& Cst1, FFVar const& Var2 )
{
  FFBase *pdag = Var2._dag;
  FFVar *pVar2 = Var2._ops.first->pres;
  FFVar* pCst1 = pdag->_add_constant( Cst1 );
  FFVar *pVar1 = pCst1->_ops.first->pres;

  FFOp *pOp = pdag->_insert_operation( top, pVar1, pVar2 );
  if( pOp->pres ) return *pOp->pres;
  pVar1->_ops.second.push_back( pOp );
  pVar2->_ops.second.push_back( pOp );
  pOp->pres = pdag->_add_auxiliary( dep, pOp );
  return *pOp->pres;
}

template <typename U> inline FFVar&
FFBase::_insert_binary_operation
( int const top, FFDep const& dep, FFVar const& Var1, U const& Cst2 )
{
  FFBase *pdag = Var1._dag;
  FFVar *pVar1 = Var1._ops.first->pres;
  FFVar* pCst2 = pdag->_add_constant( Cst2 );
  FFVar *pVar2 = pCst2->_ops.first->pres;

  FFOp *pOp = pdag->_insert_operation( top, pVar1, pVar2 );
  if( pOp->pres ) return *pOp->pres;
  pVar1->_ops.second.push_back( pOp );
  pVar2->_ops.second.push_back( pOp );
  pOp->pres = pdag->_add_auxiliary( dep, pOp );
  return *pOp->pres;
}

inline FFVar&
FFBase::_insert_unary_operation
( int const top, FFDep const& dep, FFVar const& Var )
{
  FFBase* pdag = Var._dag;
  FFVar *pVar = Var._ops.first->pres;

  FFOp* pOp = pdag->_insert_operation( top, pVar );
  if( pOp->pres ) return *pOp->pres;
  pVar->_ops.second.push_back( pOp );
  pOp->pres = pdag->_add_auxiliary( dep, pOp );
  return *pOp->pres;
}

inline FFOp*
FFBase::_insert_operation
( int const top, FFVar* lop )
{
  FFOp* op = new FFOp( top, lop, 0 );
  typename FFBase::it_Ops itop = _Ops.find( op );
  if( itop!=_Ops.end() ){ delete op; return *itop; }
  _Ops.insert( op );
  return op;
}

inline FFOp*
FFBase::_insert_operation
( int const top, FFVar* lop, FFVar* rop )
{
  FFOp* op = new FFOp( top, lop, rop, 0 );
  typename FFBase::it_Ops itop = _Ops.find( op );
  if( itop!=_Ops.end() ){ delete op; return *itop; }
  _Ops.insert( op );
  return op;
}

inline FFOp*
FFBase::_insert_operation
( int const top, unsigned const nop, FFVar** ops )
{
  FFOp* op = new FFOp( top, nop, ops, 0 );
  typename FFBase::it_Ops itop = _Ops.find( op );
  if( itop!=_Ops.end() ){ delete op; return *itop; }
  _Ops.insert( op );
  return op;
}

inline bool
FFBase::_remove_operation
( FFOp* op )
{
  typename FFBase::it_Ops itop = _Ops.find( op );
  if( itop == _Ops.end() ) return false;
  delete op;
  _Ops.erase( itop );
  return true;
}
 
inline FFVar*
FFBase::_add_auxiliary
( FFDep const& dep, FFOp* pOp )
{
  pOp->pres = new FFVar( this, dep, pOp );
  _append_aux( pOp->pres );
  return pOp->pres;
}

inline FFVar*
FFBase::_set_constant
( FFVar const* pVar, FFNum const& num )
{
  it_Vars itVar = _Vars.find( const_cast<FFVar*>(pVar) );
  if( itVar!=_Vars.end() ){
    (*itVar)->_num = num;
    (*itVar)->_cst = true;
    (*itVar)->_ops.first->type = FFOp::CNST;
  }
  return *itVar;
}

inline FFVar*
FFBase::_unset_constant
( FFVar const* pVar )
{
  it_Vars itVar = _Vars.find( const_cast<FFVar*>(pVar) );
  if( itVar!=_Vars.end() ){
    (*itVar)->_cst=false;
    (*itVar)->_ops.first->type = FFOp::VAR;
  }
  return *itVar;
}

inline FFVar*
FFBase::_add_constant
( double const x )
{
  // Check if real constant x already defined in _Vars
  FFVar* pAux = new FFVar( x );
  it_Vars iAux = _Vars.find( pAux );
  if( iAux!=_Vars.end() ){ delete pAux; return *iAux; }

  // Otherwise, append constant x
  _append_aux( pAux, FFOp::CNST );
  return pAux;
}

inline FFVar*
FFBase::_add_constant
( int const n )
{
  // Check if integer constant n already defined in _Vars
  FFVar* pAux = new FFVar( n );
  it_Vars iAux = _Vars.find( pAux );
  if( iAux!=_Vars.end() ){ delete pAux; return *iAux; }

  // Otherwise, append constant n
  _append_aux( pAux, FFOp::CNST );
  return pAux;
}

inline void
FFBase::_append_aux
( FFVar* pAux, int const tOp )
{
  FFOp*pOp = new FFOp( tOp, 0, pAux );
  _Ops.insert( pOp );
  pAux->dag() = this;
  pAux->ops().first = pOp;
  pAux->id().second = _naux++;
  _append_aux( pAux );
}

inline void
FFBase::_append_aux
( FFVar* pAux )
{
  _Vars.insert( pAux );
}

inline void
FFBase::_append_var
( FFVar* pVar )
{
  _Vars.insert( pVar );
}   

inline FFVar*
FFBase::_find_var
( typename FFVar::pt_idVar const& id )
{
  if( id.second == FFVar::NOREF ) return nullptr; // Prevents returning a zero constant for unreferenced variables
  _dummyVar.id() = id;
  it_Vars iVar = _Vars.find( &_dummyVar );
  return( iVar == _Vars.end()? nullptr: *iVar );
}

inline FFVar
FFBase::_create_var
( typename FFVar::pt_idVar const& id, std::string const& name )
{
  return FFVar( this, id, name );
}   

inline bool
FFBase::_get_constant
( FFVar const*& pVar )
{
  if( !pVar->cst() ) return false;
  FFNum const& num = pVar->num();
  switch( num.t ){
    case FFNum::INT:  pVar = _add_constant( num.n ); break;
    case FFNum::REAL: pVar = _add_constant( num.x ); break;
  }
  return true;
}

inline FFSubgraph
FFBase::subgraph
( std::vector<FFVar const*> const& vDep )
{
  _reset_operations();
  FFSubgraph sgDep;
  for( auto const& dep : vDep ){
    FFVar const* pVar = dep;
    if( !pVar->ops().first ) assert( _get_constant( pVar ) );
    pVar->ops().first->propagate_subgraph( sgDep.l_op );
    sgDep.set_dep( pVar->ops().first->iflag );
  }
  sgDep.set_move();
  return sgDep;
}

inline FFSubgraph
FFBase::subgraph
( unsigned int const nDep, FFVar const* pDep )
{
  _reset_operations();
  FFSubgraph sgDep;
  for( unsigned int i=0; i<nDep; i++ ){
    FFVar const* pVar = &pDep[i];
    if( !pVar->ops().first ) assert( _get_constant( pVar ) );
    pVar->ops().first->propagate_subgraph( sgDep.l_op );
    sgDep.set_dep( pVar->ops().first->iflag );
  }
  sgDep.set_move();
  return sgDep;
}

inline FFSubgraph
FFBase::subgraph
( std::set<unsigned> const& ndxDep, FFVar const* pDep )
{
  _reset_operations();
  FFSubgraph sgDep;
  for( unsigned const& i : ndxDep ){
    FFVar const* pVar = &pDep[i];
    if( !pVar->ops().first ) assert( _get_constant( pVar ) );
    pVar->ops().first->propagate_subgraph( sgDep.l_op );
    sgDep.set_dep( pVar->ops().first->iflag );
  }
  sgDep.set_move();
}

template< typename V, typename COMP>
inline FFSubgraph
FFBase::subgraph
( std::map<V,FFVar,COMP> const& mDep )
{
  _reset_operations();
  FFSubgraph sgDep;
  for( auto const& iDep : mDep ){
    FFVar const* pVar = &iDep.second;
    if( !pVar->ops().first ) assert( _get_constant( pVar ) );
    pVar->ops().first->propagate_subgraph( sgDep.l_op );
    sgDep.set_dep( pVar->ops().first->iflag );
  }
  sgDep.set_move();
  return sgDep;
}

inline void
FFBase::output
( FFSubgraph const& sgDep, std::string const& header,
  std::ostream&os )
{
  if( sgDep.l_op.empty() ){
    os << "\nEMPTY SUBGRAPH" << header.c_str() << "\n";
    return;
  }
  os << "\nOPERATIONS IN SUBGRAPH" << header.c_str() << ":\n";
  for( auto const& [op,mov] : sgDep.l_op )
    os << "  " << *op->pres << "\t" << (mov?"<<  ":"<-  ") << *op << std::endl;
  os << "DEPENDENTS IN SUBGRAPH" << header.c_str() << ":\n";
  unsigned idep = 0;
  for( auto const& op: sgDep.op_dep )
    os << "  " << idep++ << ":  " << *(op->pres) << std::endl;
}

inline void
FFBase::dot_script
( const unsigned int nDep, const FFVar*pDep, std::ostream&os ) const
{
  _reset_operations();
  os << "\ndigraph G {\n";
  for( unsigned int i=0; i<nDep; i++ ){
    if( !pDep[i].ops().first ) continue;
    pDep[i].ops().first->generate_dot_script( os );
  }
  os << "}\n";
}

inline void
FFBase::dot_script
( const std::vector<const FFVar*>&vDep, std::ostream&os ) const
{
  _reset_operations();
  os << "\ndigraph G {\nnode [shape=record];\n";
  typename std::vector<const FFVar*>::const_iterator it = vDep.begin();
  for( ; it!=vDep.end(); ++it ){
    if( !(*it)->ops().first ) continue;
    (*it)->ops().first->generate_dot_script( os );
  }
  os << "}\n";
}

template< typename U>
inline U
FFBase::sum
( unsigned const n, U const* V, double const* a )
{
  if( !n || !V ) return 0;
  if( !a ){
    U sumV = V[0];
    for( unsigned i=1; i<n; i++ ) sumV += V[i];
    return sumV;
  }
  U sumaV = a[0]*V[0];
  for( unsigned i=1; i<n; i++ ) sumaV += a[i]*V[i];
  return sumaV;
}

template< typename U>
inline U
FFBase::prod
( unsigned const n, U const* V )
{
  if( !n || !V ) return 0;
  U prodV = V[0];
  for( unsigned i=1; i<n; i++ ) prodV *= V[i];
  return prodV;
}

template< typename U>
inline U
FFBase::trace
( unsigned const n, U const* A )
{
  if( !n || !A ) return 0;
  U trA = A[0];
  for( unsigned i=1; i<n; i++ ) trA += A[i*n+i];
  return trA;
}

template< typename U>
inline U*
FFBase::sum
( unsigned const m, unsigned const n,
  U const *A, U const* B )
{
  if( !n || !m || !A || !B ) return 0;
  U* AB = new U[m*n];
  FFBase::sum( m, n, A, B, AB );
  return AB;
}

template< typename U>
inline void
FFBase::sum
( unsigned const m, unsigned const n, U const* A,
  U const* B, U* AB )
{
  if( !m || !n || !A || !B ) return;
  assert( AB );
  for( unsigned i=0; i<m; i++ )
    for( unsigned j=0; j<n; j++ )
      AB[i+j*m] = A[i+j*m] + B[i+j*m];
}

template< typename U>
inline U*
FFBase::sub
( unsigned const m, unsigned const n,
  U const* A, U const* B )
{
  if( !n || !m || !A || !B ) return 0;
  U* AB = new U[m*n];
  FFBase::sub( m, n, A, B, AB );
  return AB;
}

template< typename U>
inline void
FFBase::sub
( unsigned const m, unsigned const n, U const* A,
  U const* B, U* AB )
{
  if( !m || !n || !A || !B ) return;
  assert( AB );
  for( unsigned i=0; i<m; i++ )
    for( unsigned j=0; j<n; j++ )
      AB[i+j*m] = A[i+j*m] - B[i+j*m];
}

template< typename U>
inline U*
FFBase::prod
( unsigned const m, unsigned const n, unsigned const p,
  U const* A, U const* B )
{
  if( !m || !n || !p || !A || !B ) return 0;
  U* AB = new U[m*p];
  FFBase::prod( m, n, p, A, B, AB );
  return AB;
}

template< typename U>
inline void
FFBase::prod
( unsigned const m, unsigned const n, unsigned const p,
  U const* A, U const* B, U* AB )
{
  if( !n || !m || !p || !A || !B ) return;
  assert( AB );
  for( unsigned i=0; i<m; i++ ){
    for( unsigned j=0; j<p; j++ ){
      AB[i+j*m] = 0.;
      for( unsigned k=0; k<n; k++ )
        AB[i+j*m] += A[i+k*m] * B[k+j*n];
    }
  }
}

template< typename U>
inline U*
FFBase::inv
( unsigned const n, U const* A )
{
  if( !n || !A ) return 0;

  // Initialization
  U* Bkm1 = new U[n*n];
  U* Bk   = new U[n*n];
  for( unsigned i=0; i<n*n; ++i ) Bkm1[i] = A[i];
  U* c    = new U[n];
  
  c[0] = FFBase::trace( n, Bkm1 );

  // Main Loop  
  for( unsigned l=1; l<n; ++l, std::swap( Bk, Bkm1 ) ){
    // Auxmat = (B_{k-1} - c_{k_1}*I)  
    for( unsigned i=0; i<n; ++i ) Bkm1[i+i*n] -= c[l-1];
	// B_k = Jacobian*Auxmat
    FFBase::prod( n, n, n, A, Bkm1, Bk );
	// pcoeff_stack[l] = tr(Bk) 
    c[l] = FFBase::trace( n, Bk ) / ( l+1 );
  }
  for( unsigned i=0; i<n*n; ++i ) Bkm1[i] /= c[n-1];

  delete[] Bk;
  delete[] c;
  return Bkm1;
}

template< typename U>
inline U*
FFBase::polchar
( unsigned const n, U const* A )
{
  if( !n || !A ) return 0;

  // Initialization
  std::vector<U> Bkm1( A, A+n*n ), Bk( n*n );
  U* c = new U[n];
  c[0] = FFBase::trace( n, Bkm1.data() );

  // Main Loop  
  for( unsigned l=1; l<n; ++l, Bk.swap( Bkm1 ) ){
    // Auxmat = (B_{k-1} - c_{k_1}*I)  
    for( unsigned i=0; i<n; ++i ) Bkm1[i+i*n] -= c[l-1];
	// B_k = Jacobian*Auxmat
    FFBase::prod( n, n, n, A, Bkm1.data(), Bk.data() );
	// pcoeff_stack[l] = tr(Bk) 
    c[l] = FFBase::trace( n, Bk.data() ) / (double)( l+1 );
  }
  return c;
}

template< typename U>
inline U
FFBase::det
( unsigned const n, U const* A )
{
  if( !n || !A ) return 0.;
  U* c = FFBase::polchar( n, A );
  U d = ( n%2? c[n-1]: -c[n-1] );
  delete[] c;
  return d;
}

///////////////////////////////// FFGraph //////////////////////////////////////

template <typename... ExtOps>
template <typename... Deps> 
inline FFVar*
FFGraph<ExtOps...>::FAD
( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
  const FFVar* const pIndep, Deps... args )
{
  if( !nDep || !nIndep ) return 0;  // Nothing to do!
  assert( pDep && pIndep );

  std::vector<const FFVar*> vDep, vIndep;
  for( unsigned i=0; i<nDep; i++ )   vDep.push_back( pDep+i );
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );
  auto vDep_F = FAD( vDep, vIndep, args... );

  FFVar* pDep_F = new FFVar[ vDep_F.size() ];
  auto it = vDep_F.begin();
  for( unsigned k=0; it!=vDep_F.end(); ++it, k++ ) pDep_F[k] = **it;
  return pDep_F;
}

template <typename... ExtOps>
template <typename... Deps>
inline std::vector<const FFVar*>
FFGraph<ExtOps...>::FAD
( const std::vector<const FFVar*>&vDep, std::vector<const FFVar*>&vIndep,
  const unsigned nIndep, const FFVar* const pIndep, Deps... args )
{
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );
  return FAD( vDep, vIndep, args... );
}

template <typename... ExtOps>
template <typename... Deps> 
inline FFVar*
FFGraph<ExtOps...>::DFAD
( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
  const FFVar* const pIndep, const FFVar* const pDir, Deps... args )
{
  if( !nDep || !nIndep ) return 0;  // Nothing to do!
  assert( pDep && pIndep );

  std::vector<const FFVar*> vDep, vIndep, vDir;
  for( unsigned i=0; i<nDep; i++ ) vDep.push_back( pDep+i );
  for( unsigned i=0; i<nIndep; i++ ){ vIndep.push_back( pIndep+i ); vDir.push_back( pDir+i ); }
  std::vector<const FFVar*> vDep_F = FAD( vDep, vIndep, vDir, args... );

  FFVar* pDep_F = new FFVar[ vDep_F.size() ];
  auto it = vDep_F.begin();
  for( unsigned k=0; it!=vDep_F.end(); ++it, k++ ) pDep_F[k] = **it;
  return pDep_F;
}

template <typename... ExtOps>
template <typename... Deps>
inline std::vector<const FFVar*>
FFGraph<ExtOps...>::FAD
( const std::vector<const FFVar*>&vDep, std::vector<const FFVar*>&vIndep,
  std::vector<const FFVar*>&vDir, const unsigned nIndep,
  const FFVar* const pIndep, const FFVar* const pDir, Deps... args )
{
  for( unsigned i=0; i<nIndep; i++ ){ vIndep.push_back( pIndep+i ); vDir.push_back( pDir+i ); }
  return FAD( vDep, vIndep, vDir, args... );
}

template <typename... ExtOps>
inline std::vector<const FFVar*>
FFGraph<ExtOps...>::FAD
( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
  const std::vector<const FFVar*>&vDir )
{
#ifdef MC__FFUNC_DEBUG_DFAD
  std::cout << "D =";
  for( auto it=vDir.begin(); it!=vDir.end(); ++it ) std::cout << " " << *(*it);
  std::cout << std::endl;
#endif

  auto sDep_F = SFAD( vDep, vIndep, vDir );
  const FFVar*pZero = _add_constant( 0. );
  std::vector<const FFVar*> vDep_F( (vDir.size()?1:vIndep.size())*vDep.size(), pZero );
  for( unsigned ie(0); ie<std::get<2>(sDep_F).size(); ie++ ){
    unsigned pDep_F = std::get<0>(sDep_F)[ie]*(vDir.size()?1:vIndep.size())+std::get<1>(sDep_F)[ie];
    vDep_F[pDep_F] = std::get<2>(sDep_F)[ie];
  }
  return vDep_F;
}

template <typename... ExtOps>
inline std::vector<const FFVar*>
FFGraph<ExtOps...>::FAD
( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
  const bool transp )
{
  auto sDep_F = SFAD( vDep, vIndep );
  const FFVar*pZero = _add_constant( 0. );
  std::vector<const FFVar*> vDep_F( vIndep.size()*vDep.size(), pZero );
  for( unsigned ie(0); ie<std::get<2>(sDep_F).size(); ie++ ){
    unsigned pDep_F = transp?
                      std::get<0>(sDep_F)[ie]+vDep.size()*std::get<1>(sDep_F)[ie]:
                      vIndep.size()*std::get<0>(sDep_F)[ie]+std::get<1>(sDep_F)[ie];
    vDep_F[pDep_F] = std::get<2>(sDep_F)[ie];
  }
  return vDep_F;
}

template <typename... ExtOps>
template <typename... Deps>
inline std::tuple< unsigned, unsigned*, unsigned*, FFVar* >
FFGraph<ExtOps...>::SDFAD
( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
  const FFVar* const pIndep, const FFVar* const pDir, Deps... args )
{
  if( !nDep || !nIndep ) return std::make_tuple(0,(unsigned*)0,(unsigned*)0, (FFVar*)0);  // Nothing to do!
  assert( pDep && pIndep );

  std::vector<const FFVar*> vDep, vIndep, vDir;
  for( unsigned i=0; i<nDep; i++ ) vDep.push_back( pDep+i );
  for( unsigned i=0; i<nIndep; i++ ){ vIndep.push_back( pIndep+i ); vDir.push_back( pDir+i ); }
  auto vDep_F = SFAD( vDep, vIndep, vDir, args... );

  const unsigned nDep_F = std::get<0>(vDep_F).size();
  unsigned* iDep_F = new unsigned[ nDep_F ];
  unsigned* jDep_F = new unsigned[ nDep_F ];
  FFVar* pDep_F = new FFVar[ nDep_F ];
  for( unsigned ie=0; ie<nDep_F; ie++ ){
    iDep_F[ie] = std::get<0>(vDep_F)[ie];
    jDep_F[ie] = std::get<1>(vDep_F)[ie];
    pDep_F[ie] = *std::get<2>(vDep_F)[ie];
  }
  return std::make_tuple( nDep_F, iDep_F, jDep_F, pDep_F );
}

template <typename... ExtOps>
template <typename... Deps>
inline std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> >
FFGraph<ExtOps...>::SFAD
( const std::vector<const FFVar*>&vDep, std::vector<const FFVar*>&vIndep,
  std::vector<const FFVar*>&vDir, const unsigned nIndep,
  const FFVar* const pIndep, const FFVar* const pDir, Deps... args )
{
  for( unsigned i=0; i<nIndep; i++ ){ vIndep.push_back( pIndep+i ); vDir.push_back( pDir+i ); }
  return SFAD( vDep, vIndep, vDir, args... );
}

template <typename... ExtOps>
template <typename... Deps>
inline std::tuple< unsigned, unsigned*, unsigned*, FFVar* >
FFGraph<ExtOps...>::SFAD
( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
  const FFVar* const pIndep, Deps... args )
{
  if( !nDep || !nIndep ) return std::make_tuple(0,(unsigned*)0,(unsigned*)0, (FFVar*)0);  // Nothing to do!
  assert( pDep && pIndep );

  std::vector<const FFVar*> vDep, vIndep;
  for( unsigned i=0; i<nDep; i++ )   vDep.push_back( pDep+i );
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );
  auto vDep_F = SFAD( vDep, vIndep, args... );

  const unsigned nDep_F = std::get<0>(vDep_F).size();
  unsigned* iDep_F = new unsigned[ nDep_F ];
  unsigned* jDep_F = new unsigned[ nDep_F ];
  FFVar* pDep_F = new FFVar[ nDep_F ];
  for( unsigned ie=0; ie<nDep_F; ie++ ){
    iDep_F[ie] = std::get<0>(vDep_F)[ie];
    jDep_F[ie] = std::get<1>(vDep_F)[ie];
    pDep_F[ie] = *std::get<2>(vDep_F)[ie];
  }
  return std::make_tuple( nDep_F, iDep_F, jDep_F, pDep_F );
}

template <typename... ExtOps>
template <typename... Deps>
inline std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> >
FFGraph<ExtOps...>::SFAD
( const std::vector<const FFVar*>&vDep, std::vector<const FFVar*>&vIndep,
  const unsigned nIndep, const FFVar* const pIndep, Deps... args )
{
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );
  return SFAD( vDep, vIndep, args... );
}

template <typename... ExtOps>
inline std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> >
FFGraph<ExtOps...>::SFAD
( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
  const bool LUopt )
{
  auto sDep_F = SFAD( vDep, vIndep );
  std::tuple< std::vector<unsigned>, std::vector<unsigned>,
              std::vector<const FFVar*> > vDep_F;
  for( unsigned iv=0; iv<std::get<0>(sDep_F).size(); iv++ ){
    if( (  LUopt && std::get<0>(sDep_F)[iv] < std::get<1>(sDep_F)[iv] )
     || ( !LUopt && std::get<0>(sDep_F)[iv] > std::get<1>(sDep_F)[iv] ) )
      continue;
    std::get<0>(vDep_F).push_back( std::get<0>(sDep_F)[iv] ); // add row index
    std::get<1>(vDep_F).push_back( std::get<1>(sDep_F)[iv] ); // add column index
    std::get<2>(vDep_F).push_back( std::get<2>(sDep_F)[iv] ); // add Jacobian element
  }
  return vDep_F;
}

template <typename... ExtOps>
inline std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<FFVar const*> >
FFGraph<ExtOps...>::SFAD
( std::vector<FFVar const*> const& vDep, std::vector<FFVar const*> const& vIndep,
  std::vector<FFVar const*> const& vDir )
{
  // Nothing to do!
  if( !vIndep.size() || !vDep.size() ) return std::make_tuple( std::vector<unsigned>(),
    std::vector<unsigned>(), std::vector<FFVar const*>() );
  assert( !vDir.size() || vIndep.size() == vDir.size() );
  
  // Vector holding the results in sparse format
  std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<FFVar const*> > vDep_F;

  // Repeat forward differentiation for each dependent
  for( unsigned int i=0; i<vDep.size(); ++i ){
    if( vDep[i]->cst() ) continue;

    // Populate subgraph if empty
    auto sgDep = subgraph( 1, vDep[i] );
#ifdef MC__FFUNC_SFAD_DEBUG
    output( sgDep );
#endif
    _wkSFAD.clear();
    _wkSFAD.resize( sgDep.l_op.size() );//, fadbad::F<FFVar>() );
    std::map< unsigned, unsigned > mapIndep;

    // Count dependencies
    unsigned nIndep = 0;
    for( auto const& [op,mov] : sgDep.l_op ){
      if( op->type != FFOp::VAR ) continue;
      ++nIndep;
    }
#ifdef MC__FFUNC_SFAD_DEBUG
    std::cerr << "#independents " << nIndep << std::endl;
#endif

    // Propagate values in fadbad::F arithmetic through subgraph
#ifdef MC__FFUNC_CPU_EVAL
    double cputime = -cpuclock();
    std::cerr << "#operations " << sgDep.l_op.size() << std::endl;
#endif
    auto ito = sgDep.l_op.begin();
    auto itw = _wkSFAD.begin();
    for( unsigned int j=0; ito!=sgDep.l_op.end(); ++ito, ++itw ){
      _curOp = ito->first;

      // Initialize variable using values in l_vVar
      if( _curOp->type == FFOp::VAR ){
        *itw = *_curOp->pres;
        auto iti = vIndep.begin();
        auto itd = vDir.begin();
        for( unsigned int ii=0; iti!=vIndep.end(); ++iti, ++itd, ++ii ){
          if( _curOp->pres->id() != (*iti)->id() ) continue;
#ifdef MC__FFUNC_SFAD_DEBUG
          std::cerr << "independent " << j << ": " << itw->val() << std::endl;
#endif
          if( vDir.size() ){
            mapIndep[j] = 0;
            (*itw).diff( 0, 1 ) = **itd;
          }
          else{
            mapIndep[j] = ii;
            (*itw).diff( j++, nIndep );
          }
        }
      }

      // Evaluate current operation
      if( _curOp->type < FFOp::EXTERN )
        _curOp->evaluate( itw, _wkSFAD.data(), false );//ito->second );
      else if( !sizeof...(ExtOps) )
        throw Exceptions( Exceptions::EXTERN );
      else
        _curOp->evaluate_external( itw, _wkSFAD.data(), false, FirstExtOps(), t_NextExtOps() );
    }

    // Copy dependent values in vDep 
    auto const& op = *sgDep.op_dep.begin();
    for( unsigned j=0; j<nIndep; j++ ){
      auto pF_F = static_cast<fadbad::F<FFVar>*>( op->pres->val() );
      auto pdFdX = _find_var( pF_F->deriv(j).id() );
      if( pdFdX == nullptr ){
        auto const& num = pF_F->deriv(j).num();
        switch( num.t ){
          case FFNum::INT:  if( num.n != 0  ) pdFdX = _add_constant( num.n ); break;
          case FFNum::REAL: if( num.x != 0. ) pdFdX = _add_constant( num.x ); break;
        }
      }
      if( pdFdX == nullptr ) continue;
      std::get<0>(vDep_F).push_back( i ); // add row index (dependent)
      std::get<1>(vDep_F).push_back( mapIndep[j] ); // add column index (independent)
      std::get<2>(vDep_F).push_back( pdFdX ); // add Jacobian element
      if( vDir.size() ) break; // interrupt if directional derivatives requested
    }
  }

#ifdef MC__FFUNC_CPU_EVAL
  cputime += cpuclock();
  std::cout << "\nEvaluation time: " << std::fixed << cputime << std::endl;
#endif

  return vDep_F;
}

template <typename... ExtOps>
template <typename... Deps> 
inline FFVar*
FFGraph<ExtOps...>::BAD
( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
  const FFVar* const pIndep, Deps... args )
{
  if( !nDep || !nIndep ) return 0;  // Nothing to do!
  assert( pDep && pIndep );

  std::vector<const FFVar*> vDep, vIndep;
  for( unsigned i=0; i<nDep; i++ )   vDep.push_back( pDep+i );
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );
  auto vDep_B = BAD( vDep, vIndep, args... );

  FFVar* pDep_B = new FFVar[ vDep_B.size() ];
  auto it = vDep_B.begin();
  for( unsigned k=0; it!=vDep_B.end(); ++it, k++ ) pDep_B[k] = **it;
  return pDep_B;
}

template <typename... ExtOps>
template <typename... Deps>
inline std::vector<const FFVar*>
FFGraph<ExtOps...>::BAD
( const std::vector<const FFVar*>&vDep, std::vector<const FFVar*>&vIndep,
  const unsigned nIndep, const FFVar* const pIndep, Deps... args )
{
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );
  return BAD( vDep, vIndep, args... );
}

template <typename... ExtOps>
inline std::vector<const FFVar*>
FFGraph<ExtOps...>::BAD
( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
  const bool transp )
{
  auto sDep_B = SBAD( vDep, vIndep );
  const FFVar*pZero = _add_constant( 0. );
  std::vector<const FFVar*> vDep_B( vIndep.size()*vDep.size(), pZero );
  for( unsigned ie(0); ie<std::get<2>(sDep_B).size(); ie++ ){
    unsigned pDep_B = transp?
                      std::get<0>(sDep_B)[ie]+vDep.size()*std::get<1>(sDep_B)[ie]:
                      vIndep.size()*std::get<0>(sDep_B)[ie]+std::get<1>(sDep_B)[ie];
    vDep_B[pDep_B] = std::get<2>(sDep_B)[ie];
  }
  return vDep_B;
}

template <typename... ExtOps>
template <typename... Deps>
inline std::tuple< unsigned, unsigned*, unsigned*, FFVar* >
FFGraph<ExtOps...>::SBAD
( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
  const FFVar* const pIndep, Deps... args )
{
  if( !nDep || !nIndep ) return std::make_tuple(0,(unsigned*)0,(unsigned*)0, (FFVar*)0);  // Nothing to do!
  assert( pDep && pIndep );

  std::vector<const FFVar*> vDep, vIndep;
  for( unsigned i=0; i<nDep; i++ )   vDep.push_back( pDep+i );
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );
  auto vDep_B = SBAD( vDep, vIndep, args... );

  const unsigned nDep_B = std::get<0>(vDep_B).size();
  unsigned* iDep_B = new unsigned[ nDep_B ];
  unsigned* jDep_B = new unsigned[ nDep_B ];
  FFVar* pDep_B = new FFVar[ nDep_B ];
  for( unsigned ie=0; ie<nDep_B; ie++ ){
    iDep_B[ie] = std::get<0>(vDep_B)[ie];
    jDep_B[ie] = std::get<1>(vDep_B)[ie];
    pDep_B[ie] = *std::get<2>(vDep_B)[ie];
  }
  return std::make_tuple( nDep_B, iDep_B, jDep_B, pDep_B );
}

template <typename... ExtOps>
template <typename... Deps>
inline std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> >
FFGraph<ExtOps...>::SBAD
( const std::vector<const FFVar*>&vDep, std::vector<const FFVar*>&vIndep,
  const unsigned nIndep, const FFVar* const pIndep, Deps... args )
{
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );
  return SBAD( vDep, vIndep, args... );
}

template <typename... ExtOps>
inline std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> >
FFGraph<ExtOps...>::SBAD
( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
  const bool LUopt )
{
  auto sDep_B = SBAD( vDep, vIndep );
  std::tuple< std::vector<unsigned>, std::vector<unsigned>,
              std::vector<const FFVar*> > vDep_B;
  for( unsigned iv=0; iv<std::get<0>(sDep_B).size(); iv++ ){
    if( (  LUopt && std::get<0>(sDep_B)[iv] < std::get<1>(sDep_B)[iv] )
     || ( !LUopt && std::get<0>(sDep_B)[iv] > std::get<1>(sDep_B)[iv] ) )
      continue;
    std::get<0>(vDep_B).push_back( std::get<0>(sDep_B)[iv] ); // add row index
    std::get<1>(vDep_B).push_back( std::get<1>(sDep_B)[iv] ); // add column index
    std::get<2>(vDep_B).push_back( std::get<2>(sDep_B)[iv] ); // add Jacobian element
  }
  return vDep_B;
}

template <typename... ExtOps>
inline std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<FFVar const*> >
FFGraph<ExtOps...>::SBAD
( std::vector<FFVar const*> const& vDep, std::vector<FFVar const*> const& vIndep )
{
  // Nothing to do!
  if( !vIndep.size() || !vDep.size() ) return std::make_tuple( std::vector<unsigned>(),
    std::vector<unsigned>(), std::vector<FFVar const*>() );
  
  // Vector holding the results in sparse format
  std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<FFVar const*> > vDep_B;

  // Repeat backward differentiation for each dependent
  for( unsigned int i=0; i<vDep.size(); ++i ){
    if( vDep[i]->cst() ) continue;

    // Populate subgraph if empty
    auto sgDep = subgraph( 1, vDep[i] );
#ifdef MC__FFUNC_SBAD_DEBUG
    output( sgDep );
#endif
    _wkSBAD.clear();
    _wkSBAD.resize( sgDep.l_op.size() );//, fadbad::B<FFVar>() );
    std::map< unsigned, std::pair< unsigned,typename std::vector< fadbad::B<FFVar> >::iterator > > mapIndep;

    // Count dependencies
    unsigned nIndep = 0;
    for( auto const& [op,mov] : sgDep.l_op ){
      if( op->type != FFOp::VAR ) continue;
      ++nIndep;
    }
#ifdef MC__FFUNC_SBAD_DEBUG
    std::cerr << "#independents " << nIndep << std::endl;
#endif

    // Propagate values in fadbad::F arithmetic through subgraph
#ifdef MC__FFUNC_CPU_EVAL
    double cputime = -cpuclock();
    std::cerr << "#operations " << sgDep.l_op.size() << std::endl;
#endif
    auto ito = sgDep.l_op.begin();
    auto itw = _wkSBAD.begin();
    for( unsigned int j=0; ito!=sgDep.l_op.end(); ++ito, ++itw ){
      _curOp = ito->first;

      // Initialize variable using values in l_vVar
      if( _curOp->type == FFOp::VAR ){
        *itw = *_curOp->pres;
        auto iti = vIndep.begin();
        for( unsigned int ii=0; iti!=vIndep.end(); ++iti, ++ii ){
          if( _curOp->pres->id() != (*iti)->id() ) continue;
#ifdef MC__FFUNC_SBAD_DEBUG
          std::cerr << "independent " << j << ": " << itw->val() << std::endl;
#endif
          //(*itw).diff( j++, nIndep );
          mapIndep[j++] = std::make_pair( ii, itw );
        }
      }
      
      // Evaluate current operation
      if( _curOp->type < FFOp::EXTERN )
        _curOp->evaluate( itw, _wkSBAD.data(), false );//ito->second );
      else if( !sizeof...(ExtOps) )
        throw Exceptions( Exceptions::EXTERN );
      else
        _curOp->evaluate_external( itw, _wkSBAD.data(), false, FirstExtOps(), t_NextExtOps() );
    }

    // Copy values in DepB, IndepB
    fadbad::B<FFVar> DepB = *_wkSBAD.rbegin();
    std::vector<fadbad::B<FFVar>> IndepB( nIndep );
    for( unsigned j=0; j<nIndep; j++ )
      IndepB[j] = *mapIndep[j].second;
    _wkSBAD.clear();
  
    DepB.diff( 0, 1 );
    for( unsigned j=0; j<nIndep; j++ ){
#ifdef MC__FFUNC_SBAD_DEBUG
      std::cerr << "independent " << j << ": " << IndepB[j].val() << std::endl;
#endif
      FFVar dXj = IndepB[j].d(0);
      auto pdFdX = _find_var( dXj.id() );
      if( pdFdX == nullptr ){
        auto const& num = dXj.num();
        switch( num.t ){
          case FFNum::INT:  if( num.n != 0  ) pdFdX = _add_constant( num.n ); break;
          case FFNum::REAL: if( num.x != 0. ) pdFdX = _add_constant( num.x ); break;
        }
      }
      if( pdFdX == nullptr ) continue;
      std::get<0>(vDep_B).push_back( i ); // add row index
      std::get<1>(vDep_B).push_back( mapIndep[j].first ); // add column index
      std::get<2>(vDep_B).push_back( pdFdX ); // add Jacobian element
    }
  }

#ifdef MC__FFUNC_CPU_EVAL
  cputime += cpuclock();
  std::cout << "\nEvaluation time: " << std::fixed << cputime << std::endl;
#endif

  return vDep_B;
}

template <typename... ExtOps>
inline const FFVar*
FFGraph<ExtOps...>::TAD
( const unsigned int ordermax, const unsigned nDep, const FFVar* const pDep,
  const unsigned nVar, const FFVar* const pVar, const FFVar* const pIndep )
{
  if( !nDep || !nVar ) return nullptr;  // Nothing to do!
  assert( pDep && pVar );

  std::vector<const FFVar*> vDep, vVar;
  for( unsigned i=0; i<nDep; i++ ) vDep.push_back( pDep+i );
  for( unsigned i=0; i<nVar; i++ ) vVar.push_back( pVar+i );

  std::vector<const FFVar*> vDep_T = TAD( ordermax, vDep, vVar, pIndep );
  FFVar* pDep_T = new FFVar[ vDep_T.size() ];
  typename std::vector<const FFVar*>::const_iterator it = vDep_T.begin();
  for( unsigned k=0; it!=vDep_T.end(); ++it, k++ ) pDep_T[k] = **it;

  return pDep_T;
}

template <typename... ExtOps>
inline std::vector<FFVar const*>
FFGraph<ExtOps...>::TAD
( unsigned int const ordermax, std::vector<FFVar const*> const& vDep,
  std::vector<FFVar const*> const& vVar, FFVar const* const pIndep )
{
  // Check dependent and independent vector sizes
  if( !vVar.size() || !vDep.size() ) return std::vector<FFVar const*>();
  assert( vVar.size() == vDep.size() );

  // Obtain subgraph
  auto sgDep = subgraph( vDep );
  _wkTAD.clear();
  _wkTAD.resize( sgDep.l_op.size() );

  // Vector holding the results
  std::vector<FFVar const*> vDep_T; // <- vector holding the results
  fadbad::T<FFVar>** pX_T = new fadbad::T<FFVar>*[ vVar.size() ];
  fadbad::T<FFVar>** pF_T = new fadbad::T<FFVar>*[ vDep.size() ];

  // Propagate values in fadbad::T type arithmetic through subgraph
#ifdef MC__FFUNC_CPU_EVAL
  double cputime = -cpuclock();
  std::cerr << "#operations " << sgDep.l_op.size() << std::endl;
#endif
  auto ito = sgDep.l_op.begin();
  auto itw = _wkTAD.begin();
  for( ; ito!=sgDep.l_op.end(); ++ito, ++itw ){
    _curOp = ito->first;

    // Initialize variable
    if( _curOp->type == FFOp::VAR ){
      FFVar* pXi = _curOp->pres;
      *itw = *pXi;
      // Independent variable
      if( pIndep && pXi->id() == pIndep->id() )
        (*itw)[1] = 1.;
      // Dependent variables
      auto itv = vVar.begin();
      for( unsigned int i=0; itv!=vVar.end(); ++itv, i++ ){
        if( pXi->id() != (*itv)->id() ) continue;
        pX_T[i] = &(*itw);
        // Append 0th-order Taylor coefficient of ith-dependent to result vector
        vDep_T.push_back( pXi );
#ifdef MC__FFUNC_DEBUG_TAD
        std::cout << "FFGraph::TAD *** f(" << i << ")[0] = "
                  << *pXi << "  (" << pXi << ")\n";
#endif
      }
      // Attach fadbad::T<FFVar>* variable to corresponding variable
      pXi->val() = &(*itw);
    }
    
    // Evaluate current operation
    if( _curOp->type < FFOp::EXTERN )
      _curOp->evaluate( itw, _wkTAD.data(), false );//ito->second );
    else if( !sizeof...(ExtOps) )
      throw Exceptions( Exceptions::EXTERN );
    else
      _curOp->evaluate_external( itw, _wkTAD.data(), false, FirstExtOps(), t_NextExtOps() );
  }

  // Set pointers to the dependents
  auto itd = vDep.begin();
  for( unsigned j=0; itd!=vDep.end(); ++itd, j++ ){
    FFVar* pFj = _find_var( (*itd)->id() );
    pF_T[j] = ( pFj? static_cast<fadbad::T<FFVar>*>( pFj->val() ): nullptr );
  }

  // Evaluate Taylor coefficients recursively
  for( unsigned int q=0; q<ordermax; q++ ){
    itd = vDep.begin();
    for( unsigned j=0; itd!=vDep.end(); ++itd, j++ ){
      // Case dependent is not a variable
      if( !pF_T[j] ){
        vDep_T.push_back( _add_constant( 0. ) ); continue;
      }
      // Evaluate qth-order Taylor coefficient for jth dependent
      pF_T[j]->eval(q);
      // Set result as (q+1)-th Taylor coefficient for x[i]
      FFVar Xjq = (*pF_T[j])[q] / double(q+1);
      FFVar*pXjq = _find_var( Xjq.id() );
      if( !pXjq ) switch( Xjq.num().t ){
        case FFNum::INT:
          (*pX_T[j])[q+1] = Xjq.num().n;
          vDep_T.push_back( _add_constant( Xjq.num().n ) );
          break;
        case FFNum::REAL:
          (*pX_T[j])[q+1] = Xjq.num().x;
          vDep_T.push_back( _add_constant( Xjq.num().x ) );
          break;
      }
      else{
        (*pX_T[j])[q+1] = *pXjq;
        // Append (q+1)th-order Taylor coefficient of jth-dependent to result vector
        vDep_T.push_back( pXjq );
      }
#ifdef MC__FFUNC_DEBUG_TAD
      std::cout << "FFGraph::TAD *** f(" << j << ")[" << q+1 << "] = "
                << *pXjq << "  (" << pXjq << ")\n";
#endif
    }
  }

  delete[] pX_T;
  delete[] pF_T;

#ifdef MC__FFUNC_CPU_EVAL
  cputime += cpuclock();
  std::cout << "\nEvaluation time: " << std::fixed << cputime << std::endl;
#endif

  return vDep_T;
}

template <typename... ExtOps>
template <typename DAG>
inline void
FFGraph<ExtOps...>::insert
( DAG* dag, unsigned const nDep, FFVar const* pDepIn, FFVar* pDepOut )
{
  if( !nDep ) return;
  assert( pDepOut && pDepIn );

  // Populate variable arrays
  auto sg = dag->subgraph( nDep, pDepIn );
  std::vector<FFVar> vVarIn, vVarOut;
  for( auto const& [op,mov] : sg.l_op ){
    if( op->type != FFOp::VAR ) continue;
    vVarIn.push_back( *op->pres );
    auto iVar = _Vars.find( op->pres );
    if( iVar == _Vars.end() ){
      vVarOut.push_back( _create_var( op->pres->id(), op->pres->name(true) ) );
      if( (long int)_nvar <= op->pres->id().second ) _nvar = op->pres->id().second+1;
    }
    else
      vVarOut.push_back( **iVar );
  }

  // Evaluate dependents in current DAG
  std::vector<FFVar> wk( sg.l_op.size() );
  dag->eval( sg, wk, nDep, pDepIn, pDepOut, vVarIn.size(), vVarIn.data(), vVarOut.data() );
}

template <typename... ExtOps>
template <typename... Deps>
inline FFVar*
FFGraph<ExtOps...>::compose
( std::set<unsigned> const& ndxDepOut, FFVar const* pDepOut, unsigned const nDepIn,
  FFVar const* pVarOut, FFVar const* pDepIn, Deps... args )
{
  if( ndxDepOut.empty() ) return nullptr; // Nothing to do!

  std::vector<FFVar const*> vDepOut; vDepOut.reserve( ndxDepOut.size() );
  for( unsigned const& i : ndxDepOut ) vDepOut.push_back( pDepOut+i );
  std::vector< std::pair<const FFVar*,const FFVar*> > vDepIn; vDepIn.reserve( nDepIn );
  for( unsigned i=0; i<nDepIn; i++ ) vDepIn.push_back( std::make_pair(pVarOut+i,pDepIn+i) );
  std::vector<const FFVar*> vDepComp = compose( vDepOut, vDepIn, args... );

  FFVar* pDepComp = new FFVar[ *ndxDepOut.rbegin()+1 ];
  typename std::vector<FFVar const*>::const_iterator it = vDepComp.begin();
  for( unsigned const& i : ndxDepOut ) pDepComp[i] = **(it++);
  return pDepComp;
}

template <typename... ExtOps>
template <typename... Deps>
inline FFVar*
FFGraph<ExtOps...>::compose
( unsigned const nDepOut, FFVar const* pDepOut, unsigned const nDepIn,
  FFVar const* pVarOut, FFVar const* pDepIn, Deps... args )
{
  if( !nDepOut ) return 0;
  //if( !nDepIn ) return pDepOut;  // Nothing to do!
  assert( pDepOut && pVarOut && pDepIn );

  std::vector<FFVar const*> vDepOut; vDepOut.reserve( nDepOut );
  for( unsigned i=0; i<nDepOut; i++ ) vDepOut.push_back( pDepOut+i );
  std::vector< std::pair<const FFVar*,const FFVar*> > vDepIn; vDepIn.reserve( nDepIn );
  for( unsigned i=0; i<nDepIn; i++ ) vDepIn.push_back( std::make_pair(pVarOut+i,pDepIn+i) );
  std::vector<const FFVar*> vDepComp = compose( vDepOut, vDepIn, args... );

  FFVar* pDepComp = new FFVar[ vDepComp.size() ];
  typename std::vector<FFVar const*>::const_iterator it = vDepComp.begin();
  for( unsigned k=0; it!=vDepComp.end(); ++it, k++ ) pDepComp[k] = **it;
  return pDepComp;
}

template <typename... ExtOps>
template <typename... Deps>
inline std::vector<const FFVar*>
FFGraph<ExtOps...>::compose
( std::vector<const FFVar*>&vDepOut,
  std::vector< std::pair<const FFVar*, const FFVar*> >&vDepIn,
  const unsigned nDepIn, const FFVar*pVarOut, const FFVar*pDepIn, Deps... args  )
{
  for( unsigned i=0; i<nDepIn; i++ ) vDepIn.push_back( std::make_pair(pVarOut+i,pDepIn+i) );
  return compose( vDepOut, vDepIn, args... );
}

template <typename... ExtOps>
inline std::vector<const FFVar*>
FFGraph<ExtOps...>::compose
( std::vector<const FFVar*>&vDepOut,
  std::vector< std::pair<const FFVar*, const FFVar*> >&vDepIn )
{
  // Check dependent and independent vector sizes
  if( !vDepIn.size() || !vDepOut.size() ) return vDepOut;

  // Propagate composition through subgraph
  auto sgDep = subgraph( vDepOut );                     // <- subgraph of current dependents
  std::vector<const FFVar*> vDepComp( vDepOut.size() ); // <- vector to hold new dependents
  std::vector<FFVar> wkDep( sgDep.l_op.size() );        // <- vector to hold intermediates
  auto itWork = wkDep.begin();
  for( auto const& [op,mov] : sgDep.l_op ){
    _curOp = op;

    // Check if _curOp is to be substituted
    bool is_set = false;
    for( auto const& sub : vDepIn ){
      if( sub.first->id() == _curOp->pres->id() ){
        *itWork = *sub.second;
        is_set = true; break;
      }
    }
    if( !is_set && _curOp->type == FFOp::VAR && !_curOp->pres->cst() ){
      *itWork = *_curOp->pres;
      is_set = true;
    }

    // (Re)evaluate current operation
    if( is_set )
      _curOp->pres->val() = &(*itWork);
    else if( _curOp->type < FFOp::EXTERN )
      _curOp->evaluate( itWork, wkDep.data(), false );
    else if( !sizeof...(ExtOps) )
      throw Exceptions( Exceptions::EXTERN );
    else
      _curOp->evaluate_external( itWork, wkDep.data(), false, FirstExtOps(), t_NextExtOps() );

    // Check for a corresponding dependent
    auto itNew = vDepComp.begin();
    for( auto const& dep : vDepOut ){
      if( dep->id() == _curOp->pres->id() ){
        auto pNew = static_cast<const FFVar*>( _curOp->pres->val() );
        *itNew = _find_var( pNew->id() );
        if( !*itNew ){
          assert( pNew->cst() );
          *itNew = _add_constant( pNew->num().val() );
        }
        break;
      }
      ++itNew;
    }
    // Increment iterator in working array wkDep
    ++itWork;
  }

  // Reset FFVar::_val field - IS THIS NECESSARY?
  //for( auto&& op : sgDep.l_op )
  //  if( op->pres->val() ) op->pres->val() = nullptr;

  return vDepComp;
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline void
FFGraph<ExtOps...>::eval
( const std::set<unsigned>&ndxDep, const FFVar*pDep, U*vDep,
  const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args )
{
  if( ndxDep.empty() ) return; // Nothing to do!

  FFSubgraph sgDep;
  return eval( sgDep, ndxDep, pDep, vDep, nVar, pVar, vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline void
FFGraph<ExtOps...>::eval
( std::vector<U>&wkDep, const std::set<unsigned>&ndxDep, const FFVar*pDep, U*vDep,
  const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args )
{
  if( ndxDep.empty() ) return; // Nothing to do!

  FFSubgraph sgDep;
  return eval( sgDep, wkDep, ndxDep, pDep, vDep, nVar, pVar, vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline void
FFGraph<ExtOps...>::eval
( FFSubgraph&sgDep, const std::set<unsigned>&ndxDep, const FFVar*pDep, U*vDep,
  const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args )
{
  if( ndxDep.empty() ) return; // Nothing to do!

  std::vector<U> wkDep; 
  return eval( sgDep, wkDep, ndxDep, pDep, vDep, nVar, pVar, vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline void
FFGraph<ExtOps...>::eval
( FFSubgraph&sgDep, std::vector<U>&wkDep, const std::set<unsigned>&ndxDep,
  const FFVar*pDep, U*vDep, const unsigned nVar, const FFVar*pVar,
  const U*vVar, Deps... args )
{
  if( ndxDep.empty() ) return; // Nothing to do!

  std::vector<FFVar> vpDep( ndxDep.size() );
  std::vector<U> vvDep( ndxDep.size() );
  std::set<unsigned>::const_iterator it = ndxDep.cbegin();
  for( unsigned iDep=0; it != ndxDep.cend(); ++it, iDep++ ) vpDep[iDep] = pDep[*it];

  eval( sgDep, wkDep, vvDep.size(), vpDep.data(), vvDep.data(), nVar, pVar, vVar, args... );

  it = ndxDep.cbegin();
  for( unsigned iDep=0; it != ndxDep.cend(); ++it, iDep++ ) vDep[*it] = vvDep[iDep];
}

template <typename... ExtOps>
template <typename U, typename V, typename COMP, typename... Deps>
 inline void
FFGraph<ExtOps...>::eval
( const std::map<V,FFVar,COMP>&pDep, std::map<V,U,COMP>&vDep,
  const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args )
{
  vDep.clear(); 
  if( pDep.empty() ) return; // Nothing to do!

  FFSubgraph sgDep;
  return eval( sgDep, pDep, vDep, nVar, pVar, vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename V, typename COMP, typename... Deps>
 inline void
FFGraph<ExtOps...>::eval
( std::vector<U>&wkDep, const std::map<V,FFVar,COMP>&pDep, std::map<V,U,COMP>&vDep,
  const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args )
{
  vDep.clear(); 
  if( pDep.empty() ) return; // Nothing to do!

  FFSubgraph sgDep;
  return eval( sgDep, wkDep, pDep, vDep, nVar, pVar, vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename V, typename COMP, typename... Deps>
 inline void
FFGraph<ExtOps...>::eval
( FFSubgraph&sgDep, const std::map<V,FFVar,COMP>&pDep, std::map<V,U,COMP>&vDep,
  const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args )
{
  vDep.clear(); 
  if( pDep.empty() ) return; // Nothing to do!

  std::vector<U> wkDep; 
  return eval( sgDep, wkDep, pDep, vDep, nVar, pVar, vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename V, typename COMP, typename... Deps>
 inline void
FFGraph<ExtOps...>::eval
( FFSubgraph&sgDep, std::vector<U>&wkDep, const std::map<V,FFVar,COMP>&pDep,
  std::map<V,U,COMP>&vDep, const unsigned nVar, const FFVar*pVar, const U*vVar,
  Deps... args )
{
  vDep.clear(); 
  if( pDep.empty() ) return; // Nothing to do!

  std::vector<FFVar> vpDep( pDep.size() );
  std::vector<U> vvDep( pDep.size() );
  auto it = pDep.cbegin();
  for( unsigned iDep=0; it != pDep.cend(); ++it, iDep++ )
    vpDep[iDep] = it->second;

  eval( sgDep, wkDep, vvDep.size(), vpDep.data(), vvDep.data(), nVar, pVar, vVar, args... );

  it = pDep.cbegin();
  for( unsigned iDep=0; it != pDep.cend(); ++it, iDep++ )
    vDep.insert( vDep.end(), std::make_pair( it->first, vvDep[iDep] ) );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline void
FFGraph<ExtOps...>::eval
( const unsigned nDep, const FFVar*pDep, U*vDep, const unsigned nVar,
  const FFVar*pVar, const U*vVar, Deps... args )
{
  // Nothing to do!
  if( !nDep ) return;

  auto sgDep = subgraph( nDep, pDep );
  return eval( sgDep, nDep, pDep, vDep, nVar, pVar, vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline void
FFGraph<ExtOps...>::eval
( std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep, U*vDep,
  const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args )
{
  // Nothing to do!
  if( !nDep ) return;
  
  FFSubgraph sgDep;
  return eval( sgDep, wkDep, nDep, pDep, vDep, nVar, pVar, vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline void
FFGraph<ExtOps...>::eval
( FFSubgraph&sgDep, const unsigned nDep, const FFVar*pDep, U*vDep,
  const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args )
{
  // Nothing to do!
  if( !nDep ) return;
  
  std::vector<U> wkDep; 
  return eval( sgDep, wkDep, nDep, pDep, vDep, nVar, pVar, vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline void
FFGraph<ExtOps...>::eval
( FFSubgraph&sgDep, std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
  U*vDep, const unsigned nVar, const FFVar*pVar, const U*vVar, Deps... args )
{
  std::list<unsigned> l_nVar;     l_nVar.push_back(nVar);
  std::list<const FFVar*> l_pVar; l_pVar.push_back(pVar);
  std::list<const U*> l_vVar;     l_vVar.push_back(vVar);
  return eval( sgDep, wkDep, nDep, pDep, vDep, l_nVar, l_pVar, l_vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline void
FFGraph<ExtOps...>::eval
( FFSubgraph&sgDep, std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
  U*vDep, std::list<unsigned>&l_nVar, std::list<const FFVar*>&l_pVar,
  std::list<const U*>&l_vVar, const unsigned nVar, const FFVar*pVar,
  const U*vVar, Deps... args )
{
  l_nVar.push_back(nVar);
  l_pVar.push_back(pVar);
  l_vVar.push_back(vVar);
  return eval( sgDep, wkDep, nDep, pDep, vDep, l_nVar, l_pVar, l_vVar, args... );
}

template <typename... ExtOps>
template <typename U>
inline void
FFGraph<ExtOps...>::eval
( FFSubgraph&sgDep, std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
  U*vDep, const std::list<unsigned>&l_nVar, const std::list<const FFVar*>&l_pVar,
  const std::list<const U*>&l_vVar, const bool add )
{
#ifdef MC__FFUNC_DEBUG_EVAL
  std::cout << "entering FFGraphExt<>::eval\n";
  std::cout << "Size of parameter pack: " << sizeof...(ExtOps) << std::endl;
#endif

  // Nothing to do!
  if( !nDep ) return;
  assert( pDep && vDep );
  const unsigned nIndep = l_nVar.size();
  assert( l_pVar.size() == nIndep && l_vVar.size() == nIndep );

  // Populate subgraph if empty
  if( sgDep.l_op.empty() ) sgDep = subgraph( nDep, pDep );
  wkDep.resize( sgDep.l_op.size() );

  // Propagate values in U arithmetic through subgraph
#ifdef MC__FFUNC_CPU_EVAL
  double cputime = -cpuclock();
  std::cerr << "#operations " << sgDep.l_op.size() << std::endl;
#endif
  auto ito = sgDep.l_op.begin();
  typename std::vector<U>::iterator itU = wkDep.begin();
  for( ; ito!=sgDep.l_op.end(); ++ito, ++itU ){
    _curOp = ito->first;

    // Initialize variable using values in l_vVar
    if( _curOp->type == FFOp::VAR ){
      FFVar* pF = 0;
      auto itnVar = l_nVar.begin(); auto itpVar = l_pVar.begin(); auto itvVar = l_vVar.begin();
      for( ; !pF && itnVar != l_nVar.end(); ++itnVar, ++itpVar, ++itvVar ){
        for( unsigned i=0; i<(*itnVar); i++ ){
          if( _curOp->pres->id() != (*itpVar)[i].id() ) continue;
          pF = _curOp->pres;
          *itU = (*itvVar)[i];
          break;
        }
      }
      if( !pF ) throw Exceptions( Exceptions::MISSVAR );
    }

    // Evaluate current operation
    if( _curOp->type < FFOp::EXTERN )
      _curOp->evaluate( itU, wkDep.data(), (this->options.USEMOVE? ito->second: false) );
    else if( !sizeof...(ExtOps) )
      throw Exceptions( Exceptions::EXTERN );
    else
      _curOp->evaluate_external( itU, wkDep.data(), (this->options.USEMOVE? ito->second: false),
                                 FirstExtOps(), t_NextExtOps() );
  }

  // Copy dependent values in vDep 
  unsigned int i=0;
  for( auto&& op : sgDep.op_dep ){
    if( !this->options.USEMOVE || !op->pres->mov() ){
      if( !add ) vDep[i++]  = *static_cast<U*>( op->pres->val() );
      else       vDep[i++] += *static_cast<U*>( op->pres->val() );
    }
    else{
      if( !add ) vDep[i++]  = std::move( *static_cast<U*>( op->pres->val() ) );
      else       vDep[i++] += std::move( *static_cast<U*>( op->pres->val() ) );    
    }
  }
  
  //std::cout << "#assigned dependents: " << curdep << std::endl;
#ifdef MC__FFUNC_CPU_EVAL
  cputime += cpuclock();
  std::cout << "\nEvaluation time: " << std::fixed << cputime << std::endl;
#endif

  return;
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline int
FFGraph<ExtOps...>::reval
( const std::set<unsigned>&ndxDep, const FFVar*pDep,
  U*vDep, const unsigned nVar, const FFVar*pVar, U*vVar, Deps... args )
{
  if( ndxDep.empty() ) return 0; // Nothing to do!

  std::vector<U> wkDep;
  return reval( wkDep, ndxDep, pDep, vDep, nVar, pVar, vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline int
FFGraph<ExtOps...>::reval
( std::vector<U>&wkDep, const std::set<unsigned>&ndxDep, const FFVar*pDep,
  U*vDep, const unsigned nVar, const FFVar*pVar, U*vVar, Deps... args )
{
  if( ndxDep.empty() ) return 0; // Nothing to do!

  FFSubgraph sgDep;
  return reval( sgDep, wkDep, ndxDep, pDep, vDep, nVar, pVar, vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline int
FFGraph<ExtOps...>::reval
( FFSubgraph&sgDep, const std::set<unsigned>&ndxDep, const FFVar*pDep,
  U*vDep, const unsigned nVar, const FFVar*pVar, U*vVar, Deps... args )
{
  if( ndxDep.empty() ) return 0; // Nothing to do!

  std::vector<U> wkDep;
  return reval( sgDep, wkDep, ndxDep, pDep, vDep, nVar, pVar, vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline int
FFGraph<ExtOps...>::reval
( FFSubgraph&sgDep, std::vector<U>&wkDep, const std::set<unsigned>&ndxDep,
  const FFVar*pDep, U*vDep, const unsigned nVar, const FFVar*pVar,
  U*vVar, Deps... args )
{
  if( ndxDep.empty() ) return 0; // Nothing to do!

  std::vector<FFVar> vpDep( ndxDep.size() );
  std::vector<U> vvDep( ndxDep.size() );
  std::set<unsigned>::const_iterator it = ndxDep.cbegin();
  for( unsigned iDep=0; it != ndxDep.cend(); ++it, iDep++ ){
    vpDep[iDep] = pDep[*it];
    vvDep[iDep] = vDep[*it];
  }
  int flag = reval( sgDep, wkDep, vvDep.size(), vpDep.data(), vvDep.data(), nVar, pVar, vVar, args... );

  it = ndxDep.cbegin();
  for( unsigned iDep=0; it != ndxDep.cend(); ++it, iDep++ ) vDep[*it] = vvDep[iDep];
  return flag;
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline int
FFGraph<ExtOps...>::reval
( const unsigned nDep, const FFVar*pDep,
  U*vDep, const unsigned nVar, const FFVar*pVar, U*vVar, Deps... args )
{
  auto sgDep = subgraph( nDep, pDep );
  std::vector<U> wkDep;
  std::list<unsigned> l_nVar;     l_nVar.push_back(nVar);
  std::list<const FFVar*> l_pVar; l_pVar.push_back(pVar);
  std::list<U*> l_vVar;           l_vVar.push_back(vVar);
  return reval( sgDep, wkDep, nDep, pDep, vDep, nVar, pVar, vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline int
FFGraph<ExtOps...>::reval
( std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
  U*vDep, const unsigned nVar, const FFVar*pVar, U*vVar, Deps... args )
{
  auto sgDep = subgraph( nDep, pDep );
  std::list<unsigned> l_nVar;     l_nVar.push_back(nVar);
  std::list<const FFVar*> l_pVar; l_pVar.push_back(pVar);
  std::list<U*> l_vVar;           l_vVar.push_back(vVar);
  return reval( sgDep, wkDep, nDep, pDep, vDep, l_nVar, l_pVar, l_vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline int
FFGraph<ExtOps...>::reval
( FFSubgraph&sgDep, const unsigned nDep, const FFVar*pDep,
  U*vDep, const unsigned nVar, const FFVar*pVar, U*vVar, Deps... args )
{
  std::vector<U> wkDep;
  std::list<unsigned> l_nVar;     l_nVar.push_back(nVar);
  std::list<const FFVar*> l_pVar; l_pVar.push_back(pVar);
  std::list<U*> l_vVar;           l_vVar.push_back(vVar);
  return reval( sgDep, wkDep, nDep, pDep, vDep, nVar, pVar, vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline int
FFGraph<ExtOps...>::reval
( FFSubgraph&sgDep, std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
  U*vDep, const unsigned nVar, const FFVar*pVar, U*vVar, Deps... args )
{
  std::list<unsigned> l_nVar;     l_nVar.push_back(nVar);
  std::list<const FFVar*> l_pVar; l_pVar.push_back(pVar);
  std::list<U*> l_vVar;           l_vVar.push_back(vVar);
  return reval( sgDep, wkDep, nDep, pDep, vDep, l_nVar, l_pVar, l_vVar, args... );
}

template <typename... ExtOps>
template <typename U, typename... Deps>
inline int
FFGraph<ExtOps...>::reval
( FFSubgraph&sgDep, std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
  U*vDep, std::list<unsigned>&l_nVar, std::list<const FFVar*>&l_pVar,
  std::list<U*>&l_vVar, const unsigned nVar, const FFVar*pVar,
  U*vVar, Deps... args )
{
  l_nVar.push_back(nVar);
  l_pVar.push_back(pVar);
  l_vVar.push_back(vVar);
  return reval( sgDep, wkDep, nDep, pDep, vDep, l_nVar, l_pVar, l_vVar, args... );
}

template <typename... ExtOps>
template <typename U>
inline int
FFGraph<ExtOps...>::reval
( FFSubgraph&sgDep, std::vector<U>&wkDep, const unsigned nDep, const FFVar*pDep,
  U*vDep, const std::list<unsigned>&l_nVar, const std::list<const FFVar*>&l_pVar,
  const std::list<U*>&l_vVar, U const& InfVal, const unsigned MAXPASS,
  const double THRESPASS )
{
  // Nothing to do!
  if( !nDep ) return 0;
  assert( pDep && vDep );
  const unsigned nIndep = l_nVar.size();
  assert( l_pVar.size() == nIndep && l_vVar.size() == nIndep );

  // Populate subgraph if empty
  if( sgDep.l_op.empty() ) sgDep = subgraph( nDep, pDep );
  auto& opDep = sgDep.l_op;
  //wkDep.resize( 2 * opDep.size() ); // forward results first; backward results second
  wkDep.assign( opDep.size(), InfVal ); // forward results first; backward results second

#ifdef MC__FFUNC_CPU_REVAL
  double cputime = -cpuclock();
  std::cerr << "#operations " << opDep.size() << std::endl;
#endif

  // Initialization of independent variables with values in l_vVar
  std::map<U*,U*> mapVar;// 1st: pointer to wkDep; 2nd: pointer to l_vVar
  auto ito = opDep.begin();
  typename std::vector<U>::iterator itU = wkDep.begin();
  for( ; ito!=opDep.end(); ++ito, ++itU ){
    _curOp = ito->first;
    if( _curOp->type == FFOp::VAR ){
      FFVar* pF = 0;
      auto itnVar = l_nVar.begin(); auto itpVar = l_pVar.begin(); auto itvVar = l_vVar.begin();
      for( ; !pF && itnVar != l_nVar.end(); ++itnVar, ++itpVar, ++itvVar ){
        for( unsigned i=0; i<(*itnVar); i++ ){
          if( _curOp->pres->id() != (*itpVar)[i].id() ) continue;
          pF = _curOp->pres;
          *itU = (*itvVar)[i];
          if( MAXPASS ) mapVar[&*itU] = &(*itvVar)[i];
#ifdef MC__REVAL_DEBUG
          std::cout << "Independent " << *pF << ": " << *itU << std::endl;
#endif
          break;
        }
      }
      if( !pF ) throw Exceptions( Exceptions::MISSVAR );
    }
  }

  // Repeat forward/backward pass
  std::vector<U> curVar; curVar.resize(2); // Variables before backward tightening
  unsigned ipass=0;
  for( ; ipass<MAXPASS; ++ipass ){
#ifdef MC__REVAL_DEBUG
    std::cout << "\nPASS #" << ipass << std::endl;
#endif

    // Forward propagation in U arithmetic through subgraph
    ito = opDep.begin();
    itU = wkDep.begin();
    bool is_feasible = true;
    for( ; ito != opDep.end(); ++ito, ++itU ){
      assert( itU != wkDep.end() );
      // Evaluate current operation
      _curOp = ito->first;
      if( _curOp->type < FFOp::EXTERN ){
        try{
          if( !ipass ) _curOp->evaluate( itU, wkDep.data(), false ); // move disabled for value fields
          else if( !_curOp->tighten_forward( wkDep.data() ) ) is_feasible = false;
        }
        catch(...){ continue; }
      }
      else if( !sizeof...(ExtOps) )
        throw Exceptions( Exceptions::EXTERN );
      else{
        try{
          if( !ipass ) _curOp->evaluate_external( itU, wkDep.data(), false, FirstExtOps(), t_NextExtOps() );
          else if( !_curOp->tighten_forward_external( wkDep.data(), FirstExtOps(), t_NextExtOps() ) ) is_feasible = false;
        }
        catch(...){ continue; }
      }
      if( !is_feasible ) return -ipass-1;
    }

    // Intersection of propagated dependents with vDep 
    for( unsigned i=0; i<nDep; i++ ){
      auto op = sgDep.op_dep[i];
#ifdef MC__REVAL_DEBUG
      std::cout << "Dependent " << *(op->pres) << ": " << *static_cast<U*>( op->pres->val() ) << " ^ " << vDep[i] << std::endl;
#endif
      if( !Op<U>::inter( vDep[i], *static_cast<U*>( op->pres->val() ), vDep[i] ) ){
        output( subgraph( 1, op->pres ) );
        return -ipass-1;
      }
      *static_cast<U*>( op->pres->val() ) = vDep[i];
    }
#ifdef MC__REVAL_DEBUG
    { int dum; std::cout << "PAUSED"; std::cin >> dum; }
#endif

    // Backward propagation in U arithmetic through subgraph
    bool is_tighter = false;
    for( auto rito = opDep.rbegin(); rito!=opDep.rend(); ++rito ){
      _curOp = rito->first;
      // Store current operand variables
      curVar.resize( _curOp->pops.size() );
      unsigned iop = 0;
      for( auto const& op : _curOp->pops )
        curVar[iop++] = *static_cast<U*>( op->val() );
      // Tighten current operation
      try{ if( !_curOp->tighten_backward( wkDep.data() ) ) is_feasible = false; }
      catch(...){ continue; }
#ifdef MC__REVAL_DEBUG
      std::cout << *_curOp->pres << " " << *static_cast<U*>(_curOp->pres->val()) << " = ";
      for( auto const& op : _curOp->pops )        
        std::cout << *op << " " << *static_cast<U*>(op->val()) << " @ ";
      std::cout << std::endl;
#endif
      if( !is_feasible ){
        return -ipass-1;
      }
      // Test improvement of operand variables
      iop = 0;
      for( auto const& op : _curOp->pops ){
        is_tighter = is_tighter
                  || Op<U>::l(*static_cast<U*>( op->val() )) > Op<U>::l(curVar[iop])
                     + THRESPASS*Op<U>::diam(curVar[iop])
                  || Op<U>::u(*static_cast<U*>( op->val() )) < Op<U>::u(curVar[iop])
                     - THRESPASS*Op<U>::diam(curVar[iop]);
                  //|| !Op<U>::ge( *static_cast<U*>( op->val() ), curVar[iop] );
#ifdef MC__REVAL_DEBUG
        std::cout << *static_cast<U*>( op->val() ) << "<" << curVar[iop] << "? "
                  << (Op<U>::l(*static_cast<U*>( op->val() )) > Op<U>::l(curVar[iop])
                     + THRESPASS*Op<U>::diam(curVar[iop]))
                  << (Op<U>::u(*static_cast<U*>( op->val() )) < Op<U>::u(curVar[iop])
                     - THRESPASS*Op<U>::diam(curVar[iop]))
                  //<< !Op<U>::ge( *static_cast<U*>( op->val() ), curVar[iop] )
                  << std::endl;
#endif
        ++iop;
      }
    }

    // Intersection of variable values with those in l_vVar 
    for( auto const& pvar : mapVar ){
      if( !Op<U>::inter( *pvar.second, *pvar.first, *pvar.second ) )
        return -ipass-1;
    }

    // Return when no bound betterment
    if( !is_tighter ) return ipass+1;
  }

  //std::cout << "#assigned dependents: " << curdep << std::endl;
#ifdef MC__FFUNC_CPU_REVAL
  cputime += cpuclock();
  std::cout << "\nEvaluation time: " << std::fixed << cputime << std::endl;
#endif

  return ipass;
}

template <typename... ExtOps>
template <typename U>
inline void
FFGraph<ExtOps...>::wkextract
( FFSubgraph const& sgOut, std::vector<U>& wkOut, FFSubgraph const& sgIn, std::vector<U>& wkIn )
{
  // Anything to do?
  auto const& opOut = sgOut.l_op;
  auto const& opIn  = sgIn.l_op;
  if( opOut.empty() || opIn.empty() || wkIn.size() < opIn.size() ) return;
  wkOut.resize( opOut.size() );

  // Update value fields of input operands
  auto itoi = opIn.begin();
  auto itwi = wkIn.begin();
  for( ; itoi != opIn.end(); ++itoi, ++itwi )
    itoi->first->pres->val() = &(*itwi);

  // Copy value fields to output operands
  auto itoo = opOut.begin();
  auto itwo = wkOut.begin();
  for( ; itoo != opOut.end(); ++itoo, ++itwo )
    *itwo = *static_cast<U*>( itoo->first->pres->val() );
}

#ifdef MC__USE_HSL

template <typename... ExtOps>
inline bool
FFGraph<ExtOps...>::MC13
( unsigned const nDep, FFVar const* pDep, FFVar const* pIndep,
  int* IPERM, int* IOR, int* IB, int& NB, bool const disp, std::ostream& os )
{
  // Get list of operations
  std::vector<FFVar> pVar( pIndep, pIndep+nDep );
  std::vector<FFDep> vVar( nDep );
  for( unsigned i=0; i<nDep; i++ ) vVar[i].indep(i);
  auto sgDep = subgraph( nDep, pDep );
  for( auto const& [op,mov] : sgDep.l_op ){
    // Operation not a variable
    if( op->type != FFOp::VAR ) continue;
    bool isParam = true;
    // Operation is a current independent
    for( unsigned i=0; isParam && i<nDep; i++ )
      if( op->pres->id() == pIndep[i].id() ) isParam = false;
    if( !isParam ) continue;
    // Add dummy variable for parameter
    pVar.push_back( *op->pres );
    vVar.push_back( FFDep() );
  }
  const unsigned nVar = pVar.size();
  std::vector<FFDep> wkDep( sgDep.l_op.size() );
  std::vector<FFDep> vDep( nDep );
  //for( unsigned i=0; i<nVar; i++ )
  //  std::cout << pVar[i] << " = " << vVar[i] << std::endl;
  eval( sgDep, wkDep, nDep, pDep, vDep.data(), nVar, pVar.data(), vVar.data() );

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
  mc21a_( &N, ICN.data(), &LICN, IP.data(), LENR.data(), IPERM,
          &NUMNZ, IW.data() );
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
  mc13d_( &N, ICN.data(), &LICN, IP.data(), LENR.data(), IOR,
          IB, &NB, IW.data() );
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

template <typename... ExtOps>
inline bool
FFGraph<ExtOps...>::MC33
( unsigned const nDep, FFVar const* pDep, unsigned const nIndep,
  FFVar const* pIndep, int* IP, int* IQ, int* IPROF, int* IFLAG, bool const disp,
  std::ostream& os )
{
  // Get list of operations
  std::vector<FFVar> pVar( pIndep, pIndep+nIndep );
  std::vector<FFDep> vVar( nIndep );
  for( unsigned i=0; i<nIndep; i++ ) vVar[i].indep(i);
  auto sgDep = subgraph( nDep, pDep );
  for( auto const& [op,mov] : sgDep.l_op ){
    // Operation not a variable
    if( op->type != FFOp::VAR ) continue;
    bool isParam = true;
    // Operation is a current independent
    for( unsigned i=0; isParam && i<nIndep; i++ )
      if( op->pres->id() == pIndep[i].id() ) isParam = false;
    if( !isParam ) continue;
    // Add dummy variable for parameter
    pVar.push_back( *op->pres );
    vVar.push_back( FFDep() );
  }
  const unsigned nVar = pVar.size();
  std::vector<FFDep> wkDep( sgDep.l_op.size() );
  std::vector<FFDep> vDep( nDep );
  //for( unsigned i=0; i<nVar; i++ )
  //  std::cout << pVar[i] << " = " << vVar[i] << std::endl;
  eval( sgDep, wkDep, nDep, pDep, vDep.data(), nVar, pVar.data(), vVar.data() );

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
#endif

} // namespace mc

namespace mc
{

//! @brief Specialization of the structure mc::Op for use of the type mc::FFVar as a template parameter in other MC++ types
template <> struct Op< mc::FFVar >
{
  typedef mc::FFVar FV;
  static FV point( const double c ) { return FV(c); }
  static FV zeroone() { throw typename FFBase::Exceptions( FFBase::Exceptions::UNDEF ); }
  static void I(FV& x, const FV&y) { x = y; }
  static double l(const FV& x) { throw typename FFBase::Exceptions( FFBase::Exceptions::UNDEF ); }
  static double u(const FV& x) { throw typename FFBase::Exceptions( FFBase::Exceptions::UNDEF ); }
  static double abs (const FV& x) { throw typename FFBase::Exceptions( FFBase::Exceptions::UNDEF );  }
  static double mid (const FV& x) { throw typename FFBase::Exceptions( FFBase::Exceptions::UNDEF );  }
  static double diam(const FV& x) { throw typename FFBase::Exceptions( FFBase::Exceptions::UNDEF ); }
  static FV inv (const FV& x) { return mc::inv(x);  }
  static FV sqr (const FV& x) { return mc::sqr(x);  }
  static FV sqrt(const FV& x) { return mc::sqrt(x); }
  static FV exp (const FV& x) { return mc::exp(x);  }
  static FV log (const FV& x) { return mc::log(x);  }
  static FV xlog(const FV& x) { return mc::xlog(x); }
  static FV sin (const FV& x) { return mc::sin(x);  }
  static FV cos (const FV& x) { return mc::cos(x);  }
  static FV tan (const FV& x) { return mc::tan(x);  }
  static FV asin(const FV& x) { return mc::asin(x); }
  static FV acos(const FV& x) { return mc::acos(x); }
  static FV atan(const FV& x) { return mc::atan(x); }
  static FV sinh(const FV& x) { return mc::sinh(x); }
  static FV cosh(const FV& x) { return mc::cosh(x); }
  static FV tanh(const FV& x) { return mc::tanh(x); }
  static FV erf (const FV& x) { return mc::erf(x);  }
  static FV erfc(const FV& x) { return mc::erfc(x); }
  static FV fabs(const FV& x) { return mc::fabs(x); }
  static FV fstep(const FV& x) { return mc::fstep(x); }
  static FV bstep(const FV& x) { return mc::bstep(x); }
  static FV hull(const FV& x, const FV& y) { throw typename FFBase::Exceptions( FFBase::Exceptions::UNDEF ); }
  static FV min (const FV& x, const FV& y) { return mc::min(x,y);  }
  static FV max (const FV& x, const FV& y) { return mc::max(x,y);  }
  template <typename X, typename Y> static FV pow(const X& x, const Y& y) { return mc::pow(x,y); }
  static FV cheb(const FV& x, const unsigned n) { return mc::cheb(x,n); }
  static FV prod(const unsigned int n, const FV* x) { return mc::prod(n,x); }
  static FV monom(const unsigned int n, const FV* x, const unsigned* k, const bool cheb=false) { return mc::monom(n,x,k,cheb); }
  static bool inter(FV& xIy, const FV& x, const FV& y) { xIy = mc::inter(x,y); return true; }
  static bool eq(const FV& x, const FV& y) { throw typename FFBase::Exceptions( FFBase::Exceptions::UNDEF ); }
  static bool ne(const FV& x, const FV& y) { throw typename FFBase::Exceptions( FFBase::Exceptions::UNDEF ); }
  static bool lt(const FV& x, const FV& y) { throw typename FFBase::Exceptions( FFBase::Exceptions::UNDEF ); }
  static bool le(const FV& x, const FV& y) { throw typename FFBase::Exceptions( FFBase::Exceptions::UNDEF ); }
  static bool gt(const FV& x, const FV& y) { throw typename FFBase::Exceptions( FFBase::Exceptions::UNDEF ); }
  static bool ge(const FV& x, const FV& y) { throw typename FFBase::Exceptions( FFBase::Exceptions::UNDEF ); }
};

} // namespace mc

namespace fadbad
{

//! @brief Specialization of the structure fadbad::Op for use of the type mc::FFVar as a template parameter of the classes fadbad::F, fadbad::B and fadbad::T in FADBAD++
template <> struct Op< mc::FFVar >
{
  typedef double Base;
  typedef mc::FFVar FV;
  static Base myInteger( const int i ) { return Base(i); }
  static Base myZero() { return myInteger(0); }
  static Base myOne() { return myInteger(1);}
  static Base myTwo() { return myInteger(2); }
  static double myPI() { return mc::PI; }
  static FV myPos( const FV& x ) { return  x; }
  static FV myNeg( const FV& x ) { return -x; }
  template <typename U> static FV& myCadd( FV& x, U const&  y ) { return x+=y; }
  template <typename U> static FV& myCsub( FV& x, U const&  y ) { return x-=y; }
  template <typename U> static FV& myCmul( FV& x, U const&  y ) { return x*=y; }
  template <typename U> static FV& myCdiv( FV& x, U const&  y ) { return x/=y; }
  static FV myInv( const FV& x ) { return mc::inv( x ); }
  static FV mySqr( const FV& x ) { return mc::pow( x, 2 ); }
  template <typename X, typename Y> static FV myPow( const X& x, const Y& y ) { return mc::pow( x, y ); }
  //static FV myCheb( const FV& x, const unsigned n ) { return mc::cheb( x, n ); }
  static FV mySqrt( const FV& x ) { return mc::sqrt( x ); }
  static FV myLog( const FV& x ) { return mc::log( x ); }
  static FV myExp( const FV& x ) { return mc::exp( x ); }
  static FV mySin( const FV& x ) { return mc::sin( x ); }
  static FV myCos( const FV& x ) { return mc::cos( x ); }
  static FV myTan( const FV& x ) { return mc::tan( x ); }
  static FV myAsin( const FV& x ) { return mc::asin( x ); }
  static FV myAcos( const FV& x ) { return mc::acos( x ); }
  static FV myAtan( const FV& x ) { return mc::atan( x ); }
  static FV mySinh( const FV& x ) { return mc::sinh( x ); }
  static FV myCosh( const FV& x ) { return mc::cosh( x ); }
  static FV myTanh( const FV& x ) { return mc::tanh( x ); }
  static bool myEq( const FV& x, const FV& y ) { throw std::runtime_error("fadbad::Op<FFVar>::myEq -- operation not permitted"); }
  static bool myNe( const FV& x, const FV& y ) { throw std::runtime_error("fadbad::Op<FFVar>::myNe -- operation not permitted"); }
  static bool myLt( const FV& x, const FV& y ) { throw std::runtime_error("fadbad::Op<FFVar>::myLt -- operation not permitted"); }
  static bool myLe( const FV& x, const FV& y ) { throw std::runtime_error("fadbad::Op<FFVar>::myLe -- operation not permitted"); }
  static bool myGt( const FV& x, const FV& y ) { throw std::runtime_error("fadbad::Op<FFVar>::myGt -- operation not permitted"); }
  static bool myGe( const FV& x, const FV& y ) { throw std::runtime_error("fadbad::Op<FFVar>::myGe -- operation not permitted"); }
};

} // end namespace fadbad

#endif
