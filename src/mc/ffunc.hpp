// Copyright (C) 2013-2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.

/*!
\page page_FFUNC Construction, Manipulation and Evaluation of Factorable Functions
\author Benoit Chachuat & OMEGA Research Group (http://www3.imperial.ac.uk/environmentenergyoptimisation)
\version 1.0
\date 2014
\bug No known bugs.

Originally introduced by McCormick [McCormick, 1976] for the development of a convex/concave relaxation arithmetic, <B>factorable functions</B> cover an extremely inclusive class of functions which can be represented finitely on a computer by means of a code list or a computational graph involving atom operations. These are typically unary and binary operations within a library of atom operators, which can be based for example on the C-code library <tt>math.h</tt>. Besides convex/concave relaxations, factorable functions find applications in automatic differentiation (AD) [Naumann, 2009] as well as in interval analysis [Moore <I>et al.</I>, 2009] and Taylor model arithmetic [Neumaier, 2002].

Factorable functions can be represented using <b>directed acyclic graphs (DAGs)</b>, whose nodes are subexpressions and whose directed edges are computational flows [Schichl & Neumaier, 2005]. Compared to tree-based representations, DAGs offer the essential advantage of more accurately handling the influence of subexpressions shared by several functions during evaluation.

The classes mc::FFGraph, mc::FFVar and mc::FFOp defined in <tt>ffunc.hpp</tt> implement such a DAG construction for factorable functions. They also provide a basis for their manipulation, including differentiation and Taylor expansion, as well as their evaluation, in particular with the types mc::McCormick, mc::Specbnd, mc::TVar and mc::CVar of MC++.


\section sec_FFUNC_dag How Do I Construct the DAG of a Factorable Function?

For illustration, suppose we want to construct a DAG for the factorable function \f${\bf f}:\mathbb{R}^4\to\mathbb{R}^2\f$ defined by
\f{align*}
  {\bf f} = \left(\begin{array}{c} x_2x_3-x_0\\ x_0(\exp(x_2x_3)+3.0)^4)+x_1\end{array}\right)
\f}

The constructions require the header file <tt>ffunc.hpp</tt> to be included:

\code
      #include "ffunc.hpp"
\endcode

An environment <a>mc::FFGraph</a> is first defined for recording the factorable function DAG. All four variables <a>mc::FFVar</a> participating in that function are then defined in the enviornment using the method <a>mc::FFVar::set</a>:

\code
      mc::FFGraph FF;
      const unsigned int NX = 4;
      mc::FFVar X[NX];
      for( unsigned int i=0; i<NX; i++ ) X[i].set( &FF );
\endcode

The two components of the factorable function can be defined next:

\code
      const unsigned int NF = 2;
      mc::FFVar F[NF]
        = { X[2]*X[3]-X[0],
            X[0]*pow(exp(X[2]*X[3])+3.,4)+X[1] };
      std::cout << FF;
\endcode

The last line displays the following information about the factorable function DAG:

\verbatim
    DAG VARIABLES:
      X0     => { Z1 Z7 }
      X1     => { Z8 }
      X2     => { Z0 }
      X3     => { Z0 }

    DAG INTERMEDIATES:
      Z0    <=  X2 * X3              => { Z1 Z2 }
      Z1    <=  Z0 - X0              => { }
      Z2    <=  EXP( Z0 )            => { Z4 }
      Z4    <=  Z2 + Z3              => { Z6 }
      Z6    <=  POW( Z4, Z5 )        => { Z7 }
      Z7    <=  X0 * Z6              => { Z8 }
      Z8    <=  X1 + Z7              => { }
      Z5    <=  4(I)                 => { Z6 }
      Z3    <=  3(D)                 => { Z4 }
\endverbatim

Observe that 9 auxiliary variables, \f$z_0,\ldots,z_8\f$, have been created in the DAG, which correspond to the various unary and binary operations in the factorable function expression, as well as the (integer or real) participating constants. Observe, in particular, that the common sub-expression \f$x_2x_3\f$ is detected here; that is, the intermediate \f$z_0\f$ is reused to obtain both subsequent auxiliary variables \f$z_1\f$ and \f$z_2\f$.

At this point, the member function <a>mc::FFGraph::subgraph</a> can be used to generate the subgraph of a DAG corresponding to a given subset of dependent variables. This function returns a list of const pointers <a>mc::FFOp*</a> to the operations participating in the subgraph in order of appearance, which can then be displayed using the member function <a>mc::FFGraph::output</a>:

\code
      std::list<const mc::FFOp*> op_F  = FF.subgraph( NF, F );    FF.output( op_F );
      std::list<const mc::FFOp*> op_F0 = FF.subgraph( 1, F );     FF.output( op_F0 );
\endcode

The member function <a>FFGraph::subgraph</a> 
Here, the first line generates and displays a subgraph of both components of \f${\bf f}\f$, whereas the second line generates and displays a subgraph of the first component \f$f_0\f$ only:

\verbatim
    FACTORS IN SUBGRAPH:
      X2    <=  VARIABLE
      X3    <=  VARIABLE
      Z0    <=  X2 * X3
      X0    <=  VARIABLE
      Z1    <=  Z0 - X0
      X1    <=  VARIABLE
      Z2    <=  EXP( Z0 )
      Z3    <=  3(D)
      Z4    <=  Z2 + Z3
      Z5    <=  4(I)
      Z6    <=  POW( Z4, Z5 )
      Z7    <=  X0 * Z6
      Z8    <=  X1 + Z7

    FACTORS IN SUBGRAPH:
      X2    <=  VARIABLE
      X3    <=  VARIABLE
      Z0    <=  X2 * X3
      X0    <=  VARIABLE
      Z1    <=  Z0 - X0
\endverbatim

The obtained subgraphs can also be depicted using the (open source) graph plotting program <a href="http://www.graphviz.org/">DOT</a>. The dot files <tt>F.dot</tt> and <tt>F1.dot</tt> can be generated for both subgraphs as follows [which requires the header file <tt>fstream.h</tt>]:

\code
      std::ofstream o_F( "F.dot", std::ios_base::out );
      FF.dot_script( NF, F, o_F );
      o_F.close();

      std::ofstream o_F0( "F0.dot", std::ios_base::out );
      FF.dot_script( 1, F, o_F0 );
      o_F0.close();
\endcode

The graphs can be visualized, e.g., after generating SVG files using the command line as:

\verbatim
    $ dot -Tsvg -O F.dot;  display F.dot.svg
    $ dot -Tsvg -O F0.dot; display F0.dot.svg
\endverbatim

<CENTER><TABLE BORDER=0>
<TR>
<TD><h2>Graph for file <tt>F.dot</tt></h2>\image html F.png</TD>
<TD><h2>Graph for file <tt>F0.dot</tt></h2>\image html F0.png</TD>
</TR>
</TABLE></CENTER>


\section sec_FFUNC_FADBAD How Do I Obtain the DAG of a Factorable Function's Derivatives?

Derivatives of a factorable function in mc::FFGraph can be obtained with the methods mc::FFGraph::FAD and mc::FFGraph::BAD, which implement the forward and reverse mode of automatic differentiation (AD), respectively. It should be noted that mc::FFGraph does <a>not</a> implement these AD methods per se, but uses the classes fadbad::F and fadbad::B as part of <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A>.

In the forward mode of AD, for instance, entries of the Jacobian matrix of the factorable function \f$f\f$ considered in the previous section can be added to the DAG as follows:

\code
      std::vector<const mc::FFVar*> v_X, v_F;
      for( unsigned int i=0; i<NX; i++ ) v_X.push_back( &X[i] );
      for( unsigned int j=0; j<NF; j++ ) v_F.push_back( &F[j] );
      std::vector<const mc::FFVar*> v_dFdX = FF.FAD( v_F, v_X );
      std::cout << FF;
\endcode

The last line displays the following information about the DAG of the factorable function and its Jacobian:

\verbatim
    DAG VARIABLES:
      X0     => { Z1 Z7 Z17 Z18 }
      X1     => { Z8 }
      X2     => { Z0 Z10 }
      X3     => { Z0 Z9 }

    DAG INTERMEDIATES:
      Z0    <=  X2 * X3              => { Z1 Z2 }
      Z1    <=  Z0 - X0              => { }
      Z2    <=  EXP( Z0 )            => { Z4 Z9 Z10 }
      Z4    <=  Z2 + Z3              => { Z6 Z12 }
      Z6    <=  POW( Z4, Z5 )        => { Z7 }
      Z7    <=  X0 * Z6              => { Z8 }
      Z8    <=  X1 + Z7              => { }
      Z9    <=  X3 * Z2              => { Z15 }
      Z10   <=  X2 * Z2              => { Z16 }
      Z12   <=  POW( Z4, Z11 )       => { Z14 }
      Z14   <=  Z12 * Z13            => { Z15 Z16 }
      Z15   <=  Z9 * Z14             => { Z17 }
      Z16   <=  Z10 * Z14            => { Z18 }
      Z17   <=  X0 * Z15             => { }
      Z18   <=  X0 * Z16             => { }
      Z11   <=  3(I)                 => { Z12 }
      Z5    <=  4(I)                 => { Z6 }
      Z19   <=  -1(D)                => { }
      Z20   <=  0(D)                 => { }
      Z21   <=  1(D)                 => { }
      Z3    <=  3(D)                 => { Z4 }
      Z13   <=  4(D)                 => { Z14 }
\endverbatim

Observe that 13 extra auxiliary variables, \f$z_9,\ldots,z_{21}\f$, have been created in the DAG after the application of forward AD. Moreover, the function mc:FFGraph::FAD returns a vector of pointers to the dependent variables representing the entries \f$\frac{\partial f_i}{\partial x_j}\f$ of the Jacobian matrix of \f$f\f$ in the DAG, ordered column-wise as \f$\frac{\partial f_1}{\partial x_1},\ldots,\frac{\partial f_1}{\partial x_n},\frac{\partial f_2}{\partial x_1},\ldots,\frac{\partial f_2}{\partial x_n},\ldots\f$. 

As previously, subgraphs can be constructed for all or part of the derivatives, and dot file can be generated for these subgraphs too, e.g.:


\code
      std::vector<const mc::FFVar*> v_F_dFdX = v_F; v_F_dFdX.insert( v_F_dFdX.end(), v_dFdX.begin(), v_dFdX.end() );
      std::list<const mc::FFOp*> op_F_dFdX = FF.subgraph( v_F_dFdX );
      std::ofstream o_F_dFdX( "F_dFdX.dot", std::ios_base::out );
      FF.dot_script( v_F_dFdX, o_F_dFdX );
      o_F_dFdX.close();
      
      std::vector<const mc::FFVar*> v_dF1dX3; v_dF1dX3.insert( v_dF1dX3.end(), v_dFdX[NX+3] );
      std::list<const mc::FFOp*> op_dF1dX3 = FF.subgraph( v_dF1dX3 );
      std::ofstream o_dF1dX3( "dF1dX3.dot", std::ios_base::out );
      FF.dot_script( v_dF1dX3, o_dF1dX3 );
      o_dF1dX3.close();
\endcode

The first subgraph created above corresponds to both components of the factorable function \f$f\f$ as well as all eight component of its Jacobian matrix \f$\frac{\partial {\bf f}}{\partial {\bf x}}\f$; the second subgraph is for the Jacobian element  \f$\frac{\partial f_1}{\partial x_3}\f$. The corresponding graphs are shown below.

<CENTER><TABLE BORDER=0>
<TR>
<TD><h2>Graph for file <tt>F_dFdX.dot (forward AD)</tt></h2>\image html F_dFdX_FAD.png</TD>
<TD><h2>Graph for file <tt>dF1dX3.dot (forward AD)</tt></h2>\image html dF1dX3_FAD.png</TD>
</TR>
</TABLE></CENTER>

The reserve method of AD can be applied in a likewise manner using the method mc::FFGraph::BAD instead of mc::FFGRAPH::FAD, everything else remaining the same. The corresponding graphs are shown below. Note that, in constrast to forward AD, the reverse mode only requires 11 extra auxiliary variables to construct the DAG of the Jacobian matrix.

<CENTER><TABLE BORDER=0>
<TR>
<TD><h2>Graph for file <tt>F_dFdX.dot (reverse AD)</tt></h2>\image html F_dFdX_BAD.png</TD>
<TD><h2>Graph for file <tt>dF1dX3.dot (reverse AD)</tt></h2>\image html dF1dX3_BAD.png</TD>
</TR>
</TABLE></CENTER>

\section sec_FFUNC_TAD How Do I Obtain the DAG of the Taylor Expansion of ODE solutions?

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
      mc::FFGraph FF;
      mc::FFVar X( &FF ); std::vector<const mc::FFVar*> v_X; v_X.push_back(&X);
      mc::FFVar F = X;    std::vector<const mc::FFVar*> v_F; v_F.push_back(&F);

      const unsigned int TEORDER = 10;
      std::vector<const mc::FFVar*> v_FTE = FF.TAD( TEORDER, v_F, v_X );
      std::cout << FF;

      std::list<const mc::FFOp*> op_FTE = FF.subgraph( v_FTE );
      FF.output( op_FTE );

      std::ofstream o_FTE( "FTE.dot", std::ios_base::out );
      FF.dot_script( v_FTE, o_FTE );
      o_FTE.close();
\endcode

The resulting DAG of the Taylor coefficients up to order 10 is shown below.

<CENTER><TABLE BORDER=0>
<TR>
<TD>
\verbatim
DAG VARIABLES:
  X0     => { Z1 }

DAG INTERMEDIATES:
  Z1    <=  X0 * Z0       => { Z3 }
  Z3    <=  Z1 * Z2       => { Z5 }
  Z5    <=  Z3 * Z4       => { Z7 }
  Z7    <=  Z5 * Z6       => { Z9 }
  Z9    <=  Z7 * Z8       => { Z11 }
  Z11   <=  Z9 * Z10      => { Z13 }
  Z13   <=  Z11 * Z12     => { Z15 }
  Z15   <=  Z13 * Z14     => { Z17 }
  Z17   <=  Z15 * Z16     => { }
  Z16   <=  0.1(D)        => { Z17 }
  Z14   <=  0.111111(D)   => { Z15 }
  Z12   <=  0.125(D)      => { Z13 }
  Z10   <=  0.142857(D)   => { Z11 }
  Z8    <=  0.166667(D)   => { Z9 }
  Z6    <=  0.2(D)        => { Z7 }
  Z4    <=  0.25(D)       => { Z5 }
  Z2    <=  0.333333(D)   => { Z3 }
  Z0    <=  0.5(D)        => { Z1 }

FACTORS IN SUBGRAPH:
  X0    <=  VARIABLE
  Z0    <=  0.5(D)
  Z1    <=  X0 * Z0
  Z2    <=  0.333333(D)
  Z3    <=  Z1 * Z2
  Z4    <=  0.25(D)
  Z5    <=  Z3 * Z4
  Z6    <=  0.2(D)
  Z7    <=  Z5 * Z6
  Z8    <=  0.166667(D)
  Z9    <=  Z7 * Z8
  Z10   <=  0.142857(D)
  Z11   <=  Z9 * Z10
  Z12   <=  0.125(D)
  Z13   <=  Z11 * Z12
  Z14   <=  0.111111(D)
  Z15   <=  Z13 * Z14
  Z16   <=  0.1(D)
  Z17   <=  Z15 * Z16
\endverbatim
<TD><h2>Graph for file <tt>TFE.dot</tt></h2>\image html FTE.png</TD>
</TR>
</TABLE></CENTER>

Naturally, the resulting DAG of Taylor coefficients can be differentiated using mc::FFGraph::FAD or mc::FFGraph::BAD in turn, or evaluated in any compatible arithmetic as explained next.


\section sec_FFUNC_eval How Do I Evaluate the DAG of a Factorable Function in a Given Arithmetic?

Having created the DAG of a factorable function as well as its derivatives or Taylor coefficients, one can evaluate these functions in any compatible arithmetic using the method mc::FFGraph::eval.

Coming back to our initial example, suppose that we want to compute interval bounds on the second derivatives of the factorable function 
\f{align*}
  {\bf f} = \left(\begin{array}{c} x_2x_3-x_0\\ x_0(\exp(x_2x_3)+3.0)^4)+x_1\end{array}\right)
\f}
with \f$x_0\in[0,0.5]\f$, \f$x_1\in[1,2]\f$, \f$x_2\in[-1,-0.8]\f$, and \f$x_3\in[0.5,1]\f$. For simplicity, the default interval type mc::Interval of MC++ is used here:

\code
      #include "ffunc.hpp"
      #include "interval.hpp"
      typedef mc::Interval I;
\endcode

First, a DAG of the second-order derivatives of \f$f\f$ is constructed as explained above---here using forward AD twice:

\code
      // DAG environment
      mc::FFGraph FF;

      // Independent variables
      const unsigned int NX = 4;
      mc::FFVar X[NX];
      std::vector<const mc::FFVar*> v_X;
      for( unsigned int i=0; i<NX; i++ ){
        X[i].set( &FF );
        v_X.push_back(&X[i]);
      }

      // Dependent variables
      const unsigned int NF = 2;
      mc::FFVar F[NF]
        = { X[2]*X[3]-X[0],
            X[0]*pow(exp(X[2]*X[3])+3.,4)+X[1] };
      std::vector<const mc::FFVar*> v_F;
      for( unsigned int j=0; j<NF; j++ )
        v_F.push_back( &F[j] );

      // DAG of second-order derivatives
      std::vector<const mc::FFVar*> v_d2FdX2 = FF.FAD( FF.FAD( v_F, v_X ), v_X );
\endcode

In a second step, the DAG of second-order derivatives is evaluated in real-arithmetic as follows:

\code
      // Evaluation in interval arithmetic
      I IX[NX] = { I(0,0.5), I(1,2), I(-1,-0.8), I(0.5,1) };
      std::vector< std::pair<const mc::FFVar*,I> > v_IX;
      for( unsigned int i=0; i<NX; i++ )
        v_IX.push_back( std::make_pair( &X[i], IX[i] ) );
      std::vector<I> v_Id2FdX2 = FF.eval( v_d2FdX2, v_IX );

      // Display results
      for( unsigned i=0, k=0; i<v_Id2FdX2.size(); i++ ){
        std::cout << "  d2FdX2(" << i << ") = " << v_Id2FdX2[i] << std::endl;
        if( ++k == NX ){ std::cout << std::endl; k = 0; }
      }
\endcode

The DAG evaluation can be carried out in Taylor model arithmetic likewise:

\code
      #include "tmodel.hpp"
      typedef mc::TModel<I> TM;
      typedef mc::TVar<I> TV;
\endcode

\code
      // Evaluation in 5th-order Taylor model arithmetic
      TM TM_env( NX, 5 );
      TV TVX[NX];
      std::vector< std::pair<const mc::FFVar*,TV> > v_TVX;
      for( unsigned int i=0; i<NX; i++ ){
        TVX[i].set( &TM_env, i, IX[i] );
        v_TVX.push_back( std::make_pair( &X[i], TVX[i] ) );
      }
      std::vector<TV> v_TVd2FdX2 = FF.eval( v_d2FdX2, v_TVX );

      // Display results
      for( unsigned i=0, k=0; i<v_TVd2FdX2.size(); i++ ){
        std::cout << "  d2FdX2(" << i << ") = " << v_TVd2FdX2[i].B() << std::endl;
        if( ++k == NX ){ std::cout << std::endl; k = 0; }
\endcode


These evaluations produce the following results:
<CENTER><TABLE BORDER=0>
<TR>
<TD><h2>Evaluation in Interval Arithmetic</h2>
\verbatim
  d2FdX2(0) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(1) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(2) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(3) = [  0.00000e+00 :  0.00000e+00 ]

  d2FdX2(4) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(5) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(6) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(7) = [  0.00000e+00 :  0.00000e+00 ]

  d2FdX2(8) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(9) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(10) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(11) = [  1.00000e+00 :  1.00000e+00 ]

  d2FdX2(12) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(13) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(14) = [  1.00000e+00 :  1.00000e+00 ]
  d2FdX2(15) = [  0.00000e+00 :  0.00000e+00 ]

  d2FdX2(16) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(17) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(18) = [  2.81064e+01 :  1.32573e+02 ]
  d2FdX2(19) = [ -1.32573e+02 : -4.49702e+01 ]

  d2FdX2(20) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(21) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(22) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(23) = [  0.00000e+00 :  0.00000e+00 ]

  d2FdX2(24) = [  2.81064e+01 :  1.32573e+02 ]
  d2FdX2(25) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(26) = [  0.00000e+00 :  1.02604e+02 ]
  d2FdX2(27) = [ -6.62258e+01 :  4.80507e+01 ]

  d2FdX2(28) = [ -1.32573e+02 : -4.49702e+01 ]
  d2FdX2(29) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(30) = [ -6.62258e+01 :  4.80507e+01 ]
  d2FdX2(31) = [  0.00000e+00 :  1.02604e+02 ]
\endverbatim
<TD><h2>Evaluation in Taylor Model Arithmetic</h2>
 \verbatim
  d2FdX2(0) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(1) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(2) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(3) = [  0.00000e+00 :  0.00000e+00 ]

  d2FdX2(4) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(5) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(6) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(7) = [  0.00000e+00 :  0.00000e+00 ]

  d2FdX2(8) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(9) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(10) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(11) = [  1.00000e+00 :  1.00000e+00 ]

  d2FdX2(12) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(13) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(14) = [  1.00000e+00 :  1.00000e+00 ]
  d2FdX2(15) = [  0.00000e+00 :  0.00000e+00 ]

  d2FdX2(16) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(17) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(18) = [  5.19577e+01 :  7.71431e+01 ]
  d2FdX2(19) = [ -1.14333e+02 : -5.32913e+01 ]

  d2FdX2(20) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(21) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(22) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(23) = [  0.00000e+00 :  0.00000e+00 ]

  d2FdX2(24) = [  5.19577e+01 :  7.71431e+01 ]
  d2FdX2(25) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(26) = [ -1.97513e+01 :  5.36111e+01 ]
  d2FdX2(27) = [ -1.92705e+01 :  2.57043e+01 ]

  d2FdX2(28) = [ -1.14333e+02 : -5.32913e+01 ]
  d2FdX2(29) = [  0.00000e+00 :  0.00000e+00 ]
  d2FdX2(30) = [ -1.92681e+01 :  2.57022e+01 ]
  d2FdX2(31) = [ -3.07240e+01 :  8.59907e+01 ]
\endverbatim
</TR>
</TABLE></CENTER>


\section sec_FFUNC_eval How Do I Evaluate the DAG of a Factorable Function in a Given Arithmetic?

\section sec_FFUNC_err What Errors Can I Encounter While Creating or Manipulating the DAG of a Factorable Function?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type mc::FFGraph::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during the creation/manipulation of a DAG, and then make the appropriate changes. Should an exception be thrown and not caught by the calling program, the execution will abort.

Possible errors encountered during the creation/manipulation of a DAG are:

<TABLE border="1">
<CAPTION><EM>Errors during the Creation/Manipulation of a DAG</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH><tt>1</tt> <TD>Invalid mc::FFGraph* pointer in initialization of an mc::FFVar variable
     <TR><TH><tt>2</tt> <TD>Operation between variables linked to different DAGs
     <TR><TH><tt>3</tt> <TD>Subgraph cannot be evaluated because an independent variable is missing
     <TR><TH><tt>-1</tt> <TD>Internal Error
     <TR><TH><tt>-33</tt> <TD>Feature not yet implemented in mc::FFGraph
</TABLE>

Further exceptions can be thrown by the underlying arithmetic used for the evaluation of a DAG.

\section sec_FFUNC_refs References

- McCormick, G. P., <A href="http://dx.doi.org/10.1007/BF01580665">Computability of global solutions to factorable nonconvex programs: Part I. Convex underestimating problems</A>, <i>Mathematical Programming</i>, <b>10</b>(2):147-175, 1976
- Moore, R.E., M.J. Cloud, R.B. Kearfott, <I><A href="http://books.google.co.uk/books/about/Introduction_to_interval_analysis.html?id=tT7ykKbqfEwC&redir_esc=y">"Introduction to Interval Analysis"</A></I>, SIAM, 2009
- Naumann, U., <I><A href="http://books.google.co.uk/books/about/The_Art_of_Differentiating_Computer_Prog.html?id=OgQuUR4nLu0C&redir_esc=y">The Art of Differentiating Computer Programs: An Introduction to Algorithmic
Differentiation</A></I>, SIAM, 2009
- Neumaier, A., <A href="http://dx.doi.org/10.1023/A:1023061927787">Taylor forms--Use and limits</A>, <i>Reliable Computing</i>, <b>9</b>(1):43-79, 2002
- Schichl, H., A. Neumaier, <a href="http://dx.doi.org/10.1007/s10898-005-0937-x">Interval Analysis on Directed Acyclic Graphs for Global Optimization</a>, <i>Journal of Global Optimization</i>, <b>33</b>:541–562, 2005
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
#include <utility>
#include <tuple>
#include <algorithm>
#include <stdexcept>
#include <sys/time.h>

#include "mcfadbad.hpp"
#include "mcfunc.hpp"
#include "mclapack.hpp"
#include "ffdep.hpp"

#undef  MC__FFUNC_DEBUG
#undef  MC__FFUNC_DEBUG_TAD
#undef  MC__FFUNC_DEBUG_DEPMAP
#undef  MC__FFUNC_CPU_EVAL

namespace mc
{
class FFOp;
class FFGraph;

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
  FFNum( const int i=0 ):
    t(INT), n(i)
    {}
  //! @brief Constructor for a real variable
  FFNum( const double d ):
    t(REAL), x(d)
    {}

  //! @brief Constructor for an integer scalar
  FFNum& operator=
    ( const int i )
    { t = INT; n = i; return *this; }
  //! @brief Constructor for a real scalar
  FFNum& operator=
    ( const double d )
    { t = REAL; x = d; return *this; }
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
  friend class FFGraph;
  friend struct lt_FFVar;
  friend std::ostream& operator<< ( std::ostream&, const FFGraph& );
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
  friend FFVar max ( const FFVar&, const FFVar& );
  friend FFVar max ( const unsigned int, const FFVar* );
  template <typename V> friend FFVar max ( const V&, const FFVar& );
  template <typename V> friend FFVar max ( const FFVar&, const V& );
  friend FFVar min ( const FFVar&, const FFVar& );
  friend FFVar min ( const unsigned int, const FFVar* );
  template <typename V> friend FFVar min ( const V&, const FFVar& );
  template <typename V> friend FFVar min ( const FFVar&, const V& );
  friend FFVar inter ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar inter ( const V&, const FFVar& );
  template <typename V> friend FFVar inter ( const FFVar&, const V& );
  friend FFVar inv   ( const FFVar& );
  friend FFVar sqr   ( const FFVar& );
  friend FFVar exp   ( const FFVar& );
  friend FFVar log   ( const FFVar& );
  friend FFVar sqrt  ( const FFVar& );
  friend FFVar fabs  ( const FFVar& );
  friend FFVar cos   ( const FFVar& );
  friend FFVar sin   ( const FFVar& );
  friend FFVar tan   ( const FFVar& );
  friend FFVar acos  ( const FFVar& );
  friend FFVar asin  ( const FFVar& );
  friend FFVar atan  ( const FFVar& );
  friend FFVar erf   ( const FFVar& );
  friend FFVar erfc  ( const FFVar& );
  friend FFVar fstep ( const FFVar& );
  friend FFVar bstep ( const FFVar& );
  friend FFVar pow   ( const FFVar&, const int );
  friend FFVar pow   ( const FFVar&, const double );
  friend FFVar pow   ( const FFVar&, const FFVar& );
  friend FFVar pow   ( const double, const FFVar& );
  friend FFVar cheb  ( const FFVar&, const unsigned );

public:

  // other operator overloadings
  FFVar& operator= ( const FFVar& );
  FFVar& operator= ( const int );
  FFVar& operator= ( const double );
  template <typename V> FFVar& operator+= ( const V& );
  template <typename V> FFVar& operator-= ( const V& );
  template <typename V> FFVar& operator*= ( const V& );
  template <typename V> FFVar& operator/= ( const V& );

  typedef std::list< FFOp* > t_Ops;
  //typedef typename t_Ops::iterator it_Ops;
  //typedef typename t_Ops::const_iterator cit_Ops;
  typedef typename std::pair< FFOp*, t_Ops > pt_Ops;

  /** @defgroup FFunc Construction, Manipulation and Evaluation of DAGs for Factorable Functions
   *  @{
   */
  //! @brief Index for 'free' variables in factorable function
  static const long NOREF = -33;
  //! @brief Enumeration type for variables in factorable function
  enum TYPE{
    VAR=0,	//!< Original variable
    AUX,	//!< Auxiliary variable
    CINT,	//!< Integer constant
    CREAL	//!< Real constant
  };
  //! @brief Typedef for variable identifier in factorable function
  typedef std::pair< TYPE, long > pt_idVar;
  /** @} */

private:
  //! @brief Pointer to underlying factorable function DAG - _dag := NULL for variable identifier NOREF
  mutable FFGraph *_dag;
  //! @brief Identifier (type and index)
  pt_idVar _id;
  //! @brief Numeric field (integer or real)
  mutable FFNum _num;
  //! @brief Dependence and linearity
  FFDep _dep;
  //! @brief Pointer to value field - has to be const_cast'ed in order to retreive original pointer type
  mutable void *_val;
  //! @brief Constness
  mutable bool _cst;
  //! @brief Pointer to parent (_ops.first) and children (_ops.second) operations - _ops.first := NULL for unreferenced constants
  pt_Ops _ops;

public:

  /** @ingroup FFunc
   *  @{
   */
  //! @brief Constructor for variable in DAG <a>*dag</a>
  FFVar
    ( FFGraph*dag );

  //! @brief Constructor for variable in DAG <a>*dag</a> with constant real value <a>d</a>
  FFVar
    ( FFGraph*dag, const double d );

  //! @brief Constructor for variable in DAG <a>*dag</a> with constant integer value <a>i</a>
  FFVar
    ( FFGraph*dag, const int i );

  //! @brief Attach variable to DAG <a>*dag</a>.
  FFVar& set
    ( FFGraph*dag )
    { *this = FFVar( dag );
      return *this; }

  //! @brief Attach variable to DAG <a>*dag</a>.
  FFVar& set
    ( FFGraph*dag, const double d )
    { *this = FFVar( dag, d );
      return *this; }

  //! @brief Attach variable to DAG <a>*dag</a>.
  FFVar& set
    ( FFGraph*dag, const int i )
    { *this = FFVar( dag, i );
      return *this; }

  //! @brief Constructor for integer constant
  FFVar
    ( const int i=0 )
    : _dag( 0 ), _id( CINT, NOREF ), _num( i ), _dep( i ), _val( 0 ), _cst( true )
    { _ops.first = 0; }

  //! @brief Constructor for real parameter
  FFVar
    ( const double d )
    : _dag( 0 ), _id( CREAL, NOREF ), _num( d ), _dep( d ), _val( 0 ), _cst( true )
    { _ops.first = 0; }

  //! @brief Copy constructor
  FFVar
    ( const FFVar&Var )
    : _dag( Var._dag ), _id( Var._id ), _num( Var._num ), _dep( Var._dep ),
      _val( Var._val ), _cst( Var._cst ), _ops( Var._ops )
    {}
  /** @} */

private:

  //! @brief Constructor for auxiliary variable in factorable function <a>*dag</a> defined from operation <a>*Op</a>
  FFVar
    ( FFGraph*dag, const FFDep&dep, FFOp*op=0 );

  //! @brief Constructor for a variable with identifier <a>id</a> in factorable function <a>*dag</a>
  FFVar
    ( FFGraph*dag, const pt_idVar&id );
    
public:

  /** @ingroup FFunc
   *  @{
   */
  //! @brief Get variable identifier
  const std::pair<TYPE,long> id() const
    { return _id; }

  //! @brief Get reference to variable identifier
  std::pair<TYPE,long>& id()
    { return _id; }

  //! @brief Get const reference to variable numeric field
  const FFNum& num() const
    { return _num; }

  //! @brief Get/set const reference to variable numeric field
  FFNum& num()
    { return _num; }

  //! @brief Get const reference to variable dependencies
  const FFDep& dep() const
    { return _dep; }

  //! @brief Get/set reference to variable dependencies
  FFDep& dep()
    { return _dep; }

  //! @brief Get const pointer to defining operation
  const pt_Ops ops() const
    { return _ops; }

  //! @brief Get/set pointer to defining operation
  pt_Ops& ops()
    { return _ops; }

  //! @brief Get const pointer to factorable function dag
  FFGraph* dag() const
    { return _dag; }

  //! @brief Get/set pointer to factorable function dag
  FFGraph*& dag()
    { return _dag; }

  //! @brief Get/set pointer to value field
  void*& val()
    { return _val; }

  //! @brief Get pointer to value field
  template <typename U> void reset_val
    ( const U& U_dum )
    { delete static_cast<U*>( _val ); _val = 0; }

  //! @brief Get variable name
  std::string name() const
    { return _name(_id); }

  //! @brief Get constness
  const bool cst() const
    { return _cst; }

  //! @brief Get/set constness
  bool& cst()
    { return _cst; }

  //! @brief Set variable at a constant value
  void set
    ( const double d ) const;

  ////! @brief Set variable at a constant value
  void set
    ( const int i ) const;

  ////! @brief Unset variable at a constant value
  void unset
    () const;
  /** @} */

private:

  //! @brief Return string with variable name for identifier <a>id</a>
  static std::string _name
    ( const std::pair<TYPE,long> id )
    {
      std::ostringstream ovar;
      id.first==VAR? ovar<<"X": ovar<<"Z";
      ovar<<id.second;
      return ovar.str();
    }
};
const long FFVar::NOREF;

//! @brief Structure comparing variable identifiers in a factorable function for ordering in set FFGraph::_Vars
////////////////////////////////////////////////////////////////////////
//! mc::lt_FFVar is a C++ structure comparing variable identifiers in a
//! factorable function for ordering in set FFGraph::_Vars.
////////////////////////////////////////////////////////////////////////
struct lt_FFVar
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FFVar*Var1, const FFVar*Var2 ) const
    {
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

  /** @ingroup FP
   *  @{
   */
  //! @brief Enumeration type for unary and binary operations
  enum TYPE{
    CNST=0, VAR,
    PLUS, NEG, MINUS, TIMES, SCALE, DIV,
    EXP, LOG, SQRT, SQR, IPOW, POW, SIN, COS, TAN, ASIN, ACOS, ATAN,
    FABS, ERF, FSTEP, MINF, MAXF, INTER, CHEB
  };

  //! @brief Constructor
  FFOp( TYPE top, FFVar*lop=0, FFVar*rop=0, FFVar*res=0 );

  //! @brief Destructor
  ~FFOp()
    {}

  //! @brief Type of operation
  TYPE type;
  //! @brief Pointer to operation result
  FFVar* pres;
  //! @brief Pointer to left operand
  FFVar* plop;
  //! @brief Pointer to right operand
  FFVar* prop;

  //! @brief Propagate subset of operations participating in subgraph
  void propagate_subgraph
    ( std::list<const FFOp*>&Ops ) const;
  //! @brief Reset mc::FFVar::_val field in subgraph
  template <typename U> void reset_val_subgraph
    ( const U& U_dum ) const;
  //! @brief Reset mc::FFVar::_val field in subgraph
  template <typename U> void reset_val_subgraph
    ( const U& U_dum, const std::vector<const FFVar*>&vDep,
      const std::vector<const FFVar*>&vIndep ) const;
  //! @brief Propagate script for DAG using DOT and display to <a>os</a>
  void generate_dot_script
    ( std::ostream&os ) const;
  //! @brief Append script for current operation using DOT to <a>os</a>
  void append_dot_script
    ( std::ostream&os ) const;
  //! @brief Append script for factor <a>fname</a> using DOT to <a>os</a>
  void append_dot_script_factor
    ( std::ostream&os, const std::string&fname, const bool unary,
      const unsigned int fontsize, const bool dotted=false ) const;
  //! @brief Append script for variable/contant using DOT to <a>os</a>
  void append_dot_script_variable
    ( std::ostream&os, const bool constant, const unsigned int fontsize ) const;
  //! @brief Evaluate operation in U arithmetic, dynamically allocating the result
  template <typename U> void evaluate
    ( const U& U_dum ) const;
  //! @brief Evaluate operation in U arithmetic, putting the result in <a>pUres</a>
  template <typename U> void evaluate
    ( U* pUres ) const;

  //! @brief Flag operation as visited or not
  void flag
    ( const bool visited=true ) const
    { _visited = visited; }
  //! @brief Retreive operation status (visited or not)
  const bool stat() const
    { return _visited; }
  //! @brief Retreive/set operation status (visited or not)
  const bool& stat()
    { return _visited; }
  //! @brief Return whether or not operation is univariate
  bool is_univariate() const;
  
  /** @} */

private:
  //! @brief Whether a node has been visited (during a DAG navigation)
  mutable bool _visited;
};

//! @brief C++ structure for comparing operations in a factorable program for ordering in set FFGraph::_Ops
////////////////////////////////////////////////////////////////////////
//! mc::lt_FFOp is a C++ structure for comparing operations in a
//! factorable program based on their types and operands for ordering
//! in set FFGraph::_Ops.
////////////////////////////////////////////////////////////////////////
struct lt_FFOp
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FFOp*Op1, const FFOp*Op2 ) const
    {
      // Sort by type of operation first
      if( Op1->type < Op2->type ) return true;
      if( Op1->type > Op2->type ) return false;

      // Sort by variable type next
      lt_FFVar ltVar;
      if( !Op1->plop ) return ltVar( Op1->pres, Op2->pres );
      if( ltVar( Op1->plop, Op2->plop ) ) return true;
      if( ltVar( Op2->plop, Op1->plop ) ) return false;
      if( Op1->prop ) return ltVar( Op1->prop, Op2->prop );
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
    ( const FFOp*Op1, const FFOp*Op2 ) const
    { return ( Op1->type < Op2->type ); }
};

//! @brief C++ class representing the DAG of factorable functions
////////////////////////////////////////////////////////////////////////
//! mc::FFGraph is a C++ class representing the DAG of a factorable
//! function, enabling basic manipulations on that DAG, and evaluating
//! the DAG in various arithmetics.
////////////////////////////////////////////////////////////////////////
class FFGraph
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
  friend FFVar max ( const FFVar&, const FFVar& );
  friend FFVar max ( const unsigned int, const FFVar* );
  template <typename V> friend FFVar max ( const V&, const FFVar& );
  template <typename V> friend FFVar max ( const FFVar&, const V& );
  friend FFVar min ( const FFVar&, const FFVar& );
  friend FFVar min ( const unsigned int, const FFVar* );
  template <typename V> friend FFVar min ( const V&, const FFVar& );
  template <typename V> friend FFVar min ( const FFVar&, const V& );
  friend FFVar inter ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar inter ( const V&, const FFVar& );
  template <typename V> friend FFVar inter ( const FFVar&, const V& );
  friend FFVar inv   ( const FFVar& );
  friend FFVar sqr   ( const FFVar& );
  friend FFVar exp   ( const FFVar& );
  friend FFVar log   ( const FFVar& );
  friend FFVar sqrt  ( const FFVar& );
  friend FFVar fabs  ( const FFVar& );
  friend FFVar cos   ( const FFVar& );
  friend FFVar sin   ( const FFVar& );
  friend FFVar tan   ( const FFVar& );
  friend FFVar acos  ( const FFVar& );
  friend FFVar asin  ( const FFVar& );
  friend FFVar atan  ( const FFVar& );
  friend FFVar erf   ( const FFVar& );
  friend FFVar erfc  ( const FFVar& );
  friend FFVar fstep ( const FFVar& );
  friend FFVar bstep ( const FFVar& );
  friend FFVar pow   ( const FFVar&, const int );
  friend FFVar pow   ( const FFVar&, const double );
  friend FFVar pow   ( const FFVar&, const FFVar& );
  friend FFVar pow   ( const double, const FFVar& );
  friend FFVar cheb  ( const FFVar&, const unsigned );

  // friends of this class for operator and function overloading
  friend std::ostream& operator<< ( std::ostream&, const FFGraph& );

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

  //! @brief Pointer to current operation in subtree evaluation
  const FFOp* _curOp;

public:
  /** @ingroup FFunc
   *  @{
   */
  //! @brief Default Constructor
  FFGraph():
    _nvar( 0 ), _naux( 0 ), _curOp(0)
    {}

  //! @brief Destructor
  virtual ~FFGraph()
    { clear(); }

  //! @brief Exceptions of mc::FFGraph
  class Exceptions
  {
  public:
    //! @brief Enumeration type for exception handling
    enum TYPE{
      INIT = 1,		//!< Error due to a invalid FFGraph pointer in initialization of FFVar
      DAG,		//!< Error due to an operation between variables linked to different DAGs
      INTER, 		//!< Error due to an empty intersection between constant variables
      MISSVAR,		//!< Error due to a missing independent variable for evaluating a given subgraph using FFGraph::eval
      EVAL,		//!< Error during subgraph evaluation using FFGraph::eval
      CONST,		//!< Error due to trying to attach a value field to a constant FFVAR
      INTERN = -1, 	//!< Internal error
      UNDEF = -33, 	//!< Error due to calling a function/feature not yet implemented in MC++
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr=UNDEF ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case INIT:
        return "Invalid DAG passed for initilization of a variable";
      case DAG:
        return "Operation between variables linked to different DAGs";
      case INTER:
        return "Empty intersection between constant variables";
      case MISSVAR:
        return "Missing independent variable during subgraph evaluation";
      case EVAL:
        return "Exception thrown during subgraph evaluation";
      case CONST:
        return "Trying to override a constant variable during subgraph evaluation";
      case UNDEF:
        return "Feature not yet implemented in mc::FFGraph";
      default:
        return "Undocumented error";
      }
    }
  private:
    TYPE _ierr;
  };
/*
  //! @brief Options of mc::FFGraph
  struct Options
  {
    //! @brief Constructor
    Options():
      DCDECOMPOSE(false)
      {}
    //! @brief Whether or not to append DC decomposition of product/division terms
    bool DCDECOMPOSE;
  } options;
*/
  //! @brief Current operation in DAG evaluation
  const FFOp* curOp() const
    { return _curOp; }

  //! @brief Current operation in DAG evaluation
  const FFOp*& curOp()
    { return _curOp; }

  //! @brief Number of original variables in DAG
  unsigned long nvar() const
    { return _nvar; }
  
  //! @brief Number of auxiliary variables in DAG
  unsigned long naux() const
    { return _naux; }
  
  //! @brief Reference to set of (all) variables in factorable function
  const t_Vars& Vars() const
    { return _Vars; }

  //! @brief Clear DAG (all variables and operations)
  void clear()
    { _clear_variables(); _clear_operations(); _naux=_nvar=0; }

  //! @brief Extract list of operations corresponding to <a>nDep</a> dependents in array <a>pDep</a>
  std::list<const FFOp*> subgraph
    ( const unsigned int nDep, const FFVar*pDep ) const;

  //! @brief Extract list of operations corresponding to dependents <a>vDep</a>
  std::list<const FFOp*> subgraph
    ( const std::vector<const FFVar*>&vDep ) const;

   //! @brief Create dependency map corresponding to <a>nDep</a> dependents in array <a>pDep</a> and <a>nIndep</a> independents in array <a>pIndep</a>
  CPPL::dssmatrix depmap
    ( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
      const FFVar* const pIndep );

   //! @brief Create dependency map corresponding to <a>nDep</a> dependents in array <a>pDep</a> and <a>nIndep</a> independents in array <a>pIndep</a> -- This function uses the subgraph for the dependents given in <a>opDep</a>
  CPPL::dssmatrix depmap
    ( std::list<const FFOp*>&opDep, const unsigned nDep, const FFVar* const pDep,
      const unsigned nIndep, const FFVar* const pIndep );

   //! @brief Create dependency map corresponding to the dependent variables <a>pDep</a> and independent variables <a>pIndep</a>
  CPPL::dssmatrix depmap
    ( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep ) ;

   //! @brief Create dependency map corresponding to the dependent variables <a>pDep</a> and independent variables <a>pIndep</a> -- This function uses the subgraph for the dependents given in <a>opDep</a>
  CPPL::dssmatrix depmap
    ( std::list<const FFOp*>&opDep, const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep ) ;

  //! @brief Output list of nodes in <a>Ops</a> to <a>os</a>
  static void output
    ( const std::list<const FFOp*>&Ops, std::ostream&os=std::cout );

  //! @brief Generate script for DAG visualization of factors <a>*F</a> using DOT
  void dot_script
    ( const unsigned int nDep, const FFVar*pDep, std::ostream&os=std::cout ) const;

  //! @brief Generate script for DAG visualization of factors <a>*F</a> using DOT
  void dot_script
    ( const std::vector<const FFVar*>&vDep, std::ostream&os=std::cout ) const;

  //! @brief Expand DAG with derivatives of dependents <a>vDep</a> with respect to independents <a>vIndep</a> using fadbad::F -- Returns a 3-tuple of vectors with the row index, column index, and actual Jacobian element for each non-zero entry of the Jacobian matrix
  std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> > SFAD
    ( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
      const std::vector<const FFVar*>&vDir=std::vector<const FFVar*>() );

  //! @brief Expand DAG with derivatives of dependents <a>vDep</a> with respect to independents <a>vIndep</a> using fadbad::F (directional derivatives if <a>vDir</a> is specifed) -- Returns a vector with the entries of the Jacobian matrix ordered row-wise (transp=false) or column-wise (transp=true)
  std::vector<const FFVar*> FAD
    ( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
      const bool transp );

  //! @brief Expand DAG with derivatives of dependents <a>vDep</a> with respect to independents <a>vIndep</a> using fadbad::F (directional derivatives if <a>vDir</a> is specifed) -- Returns a vector with the entries of the Jacobian matrix ordered row-wise
  std::vector<const FFVar*> FAD
    ( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
      const std::vector<const FFVar*>&vDir=std::vector<const FFVar*>() );

  //! @brief Expand DAG with derivatives of of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> using fadbad::F (directional derivatives if <a>pDir</a> is specifed) -- Returns a 4-tuple of size and arrays with the row index, column index, and actual Jacobian element for each non-zero entry of the Jacobian matrix
  std::tuple< unsigned, const unsigned*, const unsigned*, const FFVar* > SFAD
    ( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
      const FFVar* const pIndep, const FFVar* const pDir=0 );

  //! @brief Expand DAG with derivatives of of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> using fadbad::F -- Returns a 4-tuple of size and arrays with the row index, column index, and actual Jacobian element for each non-zero entry in the lower (LUopt=true) or upper (LUopt=false) traingular part of the Jacobian matrix
  std::tuple< unsigned, const unsigned*, const unsigned*, const FFVar* > SFAD
    ( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
      const FFVar* const pIndep, const bool LUopt );

  //! @brief Expand DAG with derivatives of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> using fadbad::F -- Returns an array with entries of the Jacobian matrix ordered row-wise (transp=false) or column-wise (transp=true)
  const FFVar* FAD
    ( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
      const FFVar* const pIndep, const bool transp );

  //! @brief Expand DAG with directional derivatives of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> for the direction in array <a>pDir</a> using fadbad::F -- Returns an array with entries of the Jacobian matrix ordered row-wise
  const FFVar* FAD
    ( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
      const FFVar* const pIndep, const FFVar* const pDir=0 );

  //! @brief Expand DAG with derivatives of dependents <a>vDep</a> with respect to independents <a>vIndep</a> using fadbad::B -- Returns a 3-tuple of vectors with the row index, column index, and actual Jacobian element for each non-zero entry of the Jacobian matrix
  std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> > SBAD
    ( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep );

  //! @brief Expand DAG with derivatives of dependents <a>vDep</a> with respect to independents <a>vIndep</a> using fadbad::B -- Returns a vector with the entries of the Jacobian matrix ordered row-wise
  std::vector<const FFVar*> BAD
    ( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
      const bool transp=false );

  //! @brief Expand DAG with derivatives of of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> using fadbad::B (directional derivatives if <a>pDir</a> is specifed) -- Returns a 4-tuple of size and arrays with the row index, column index, and actual Jacobian element for each non-zero entry of the Jacobian matrix
  std::tuple< unsigned, const unsigned*, const unsigned*, const FFVar* > SBAD
    ( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
      const FFVar* const pIndep=0 );

  //! @brief Expand DAG with derivatives of of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> using fadbad::B -- Returns a 4-tuple of size and arrays with the row index, column index, and actual Jacobian element for each non-zero entry in the lower (LUopt=true) or upper (LUopt=false) traingular part of the Jacobian matrix
  std::tuple< unsigned, const unsigned*, const unsigned*, const FFVar* > SBAD
    ( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
      const FFVar* const pIndep, const bool LUopt );

  //! @brief Expand DAG with derivatives of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> using fadbad::B -- Returns an array with entries of the Jacobian matrix ordered row-wise
  const FFVar* BAD
    ( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
      const FFVar* const pIndep, const bool transp=false );

  //! @brief Expand DAG with Taylor coefficients of dependents <a>vDep</a> with respect to independents <a>vIndep</a> using fadbad::T -- Same number of dependents and independent is required, e.g. for expansion of ODE solutions -- Return a vector with the 0th, 1st, ..., ordermax'th order Taylor coefficients ordered sequentially
  std::vector<const FFVar*> TAD
    ( const unsigned int ordermax, const std::vector<const FFVar*>&vDep,
      const std::vector<const FFVar*>&vVar, const FFVar* const pIndep=0 );

  //! @brief Expand DAG with Taylor coefficients of <a>nDep</a> dependents in array <a>pDep</a> with respect to <a>nIndep</a> independents in array <a>pIndep</a> using fadbad::T -- Same number of dependents and independent is required, e.g. for expansion of ODE solutions -- Returns an array with the 0th, 1st, ..., ordermax'th order Taylor coefficients ordered sequentially
  const FFVar* TAD
    ( const unsigned int ordermax, const unsigned nDep, const FFVar* const pDep,
      const unsigned nVar, const FFVar* const pVar, const FFVar* const pIndep=0 );

  //! @brief Compose the dependents in <a>vDepOut</a> with those in <a>vDepIn</a> -- This function creates the subgraph for the outer dependent variables internally
  std::vector<const FFVar*> compose
    ( const std::vector<const FFVar*>&vDepOut,
      const std::vector< std::pair<const FFVar*,const FFVar*> >&vDepIn );

  //! @brief Compose the <a>nDepOut</a> dependents in array <a>pDepOut</a> with the <a>nDepIn</a> dependents in array <a>pDepIn</a> for the variables <a>pVarOut</a> -- This function creates the subgraph for the outer dependent variables internally
  const FFVar* compose
    ( const unsigned nDepOut, const FFVar*pDepOut, const unsigned nDepIn,
      const FFVar*pVarOut,  const FFVar*pDepIn );

  //! @brief Decompose the dependent <a>pDep</a> as a linear combination of: (i) <a>nLin</a> variables <a>pLin</a> with coefficients <a>cLin</a>; and (ii) <a>nNLin</a> nonlinear terms <a>pNLin</a> with coefficients <a>cNLin</a> -- This function does not introduce extra variables in the DAG
  void decompose
    ( const FFVar&pDep, unsigned&nLin, double*cLin, FFVar*pLin,
      unsigned&nNLin, double*cNLin, FFVar*pNLin );

  //! @brief Evaluate the dependents in <a>vDep</a> in U arithmetic for the variable values specified in <a>vVar</a> -- This function creates the subgraph for the dependent variables internally
  template <typename U> std::vector<U> eval
    ( const std::vector<const FFVar*>&vDep,
      const std::vector< std::pair<const FFVar*,U> >&vVar );

  //! @brief Evaluate the dependents in <a>vDep</a> in U arithmetic for the variable values specified in <a>vVar</a> -- This function uses the subgraph for the dependent variables given in <a>opDep</a>
  template <typename U> std::vector<U> eval
    ( std::list<const FFOp*>&opDep, const std::vector<const FFVar*>&vDep,
      const std::vector< std::pair<const FFVar*,U> >&vVar );

  //! @brief Evaluate the dependents in <a>vDep</a> in U arithmetic for the variable values specified in <a>vVar</a> -- This function uses the subgraph for the dependent variables given in <a>opDep</a> as well as the preallocated array <a>opRes</a> of size <a>opDep.size()</a> to store intermediate results during the evaluation
  template <typename U> std::vector<U> eval
    ( std::list<const FFOp*>&opDep, U*opRes, const std::vector<const FFVar*>&vDep,
      const std::vector< std::pair<const FFVar*,U> >&vVar );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> in U arithmetic for the <a>nVar</a> variable in array <a>pVar</a> whose values are specified in <a>vVar</a> and write the result in <a>vDep</a> (or add the result to <a>vDep</a> if <a>add</a>==true) -- This function creates the subgraph for the dependent variables internally
  template <typename U> void eval
    ( const unsigned nDep, const FFVar*pDep, U*vDep,
      const unsigned nVar, const FFVar*pVar, const U*vVar, const bool add=false );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> in U arithmetic for the <a>nVar</a> variable in array <a>pVar</a> whose values are specified in <a>vVar</a> and write the result in <a>vDep</a> (or add the result to <a>vDep</a> if <a>add</a>==true) -- This function uses the subgraph for the dependent variables given in <a>opDep</a>
  template <typename U> void eval
    ( std::list<const FFOp*>&opDep, const unsigned nDep, const FFVar*pDep,
      U*vDep, const unsigned nVar, const FFVar*pVar, const U*vVar, const bool add=false );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> in U arithmetic for the variables in list of arrays <a>pVar</a> with values specified in <a>vVar</a> and write the result in <a>vDep</a> (or add the result to <a>vDep</a> if <a>add</a>==true) -- This function creates the subgraph for the dependent variables internally
  template <typename U> void eval
    ( const unsigned nDep, const FFVar*pDep, U*vDep, const std::list<unsigned>&nVar,
      const std::list<const FFVar*>&pVar, const std::list<const U*>&vVar, const bool add=false );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> in U arithmetic for the variables in list of arrays <a>pVar</a> with values specified in <a>vVar</a> and write the result in <a>vDep</a> (or add the result to <a>vDep</a> if <a>add</a>==true) -- This function uses the subgraph for the dependent variables given in <a>opDep</a>
  template <typename U> void eval
    ( std::list<const FFOp*>&opDep, const unsigned nDep, const FFVar*pDep,
      U*vDep, const std::list<unsigned>&nVar, const std::list<const FFVar*>&pVar,
      const std::list<const U*>&vVar, const bool add=false );

  //! @brief Evaluate the <a>nDep</a> dependents in array <a>pDep</a> in U arithmetic for the <a>nVar</a> variable in array <a>pVar</a> whose values are specified in <a>vVar</a> and write the result in <a>vDep</a> (or add the result to <a>vDep</a> if <a>add</a>==true) -- This function uses the subgraph for the dependent variables given in <a>opDep</a> as well as the preallocated array <a>opRes</a> of size <a>opDep.size()</a> to store intermediate results during the evaluation
  template <typename U> void eval
    ( std::list<const FFOp*>&opDep, U*opRes, const unsigned nDep, const FFVar*pDep,
      U*vDep, const unsigned nVar, const FFVar*pVar, U*vVar, const bool add=false );

  //! @brief Return current clock time
  static double cpuclock()
    { timeval time; gettimeofday(&time, 0) ;
      return time.tv_sec + time.tv_usec*1e-6; }

  /** @} */
   
protected:
  //! @brief Erase operation <a>op</a> in set <a>_Ops</a>
  bool _remove_operation
    ( FFOp*op );

  //! @brief Erase all operations in set <a>_Ops</a>
  void _clear_operations()
    { it_Ops ito = _Ops.begin();
      for( ; ito != _Ops.end(); ++ito ) delete *ito;
      _Ops.clear(); }

  //! @brief Reset all operations in set <a>_Ops</a>
  void _reset_operations() const
    { it_Ops ito = _Ops.begin();
      for( ; ito != _Ops.end(); ++ito ) (*ito)->flag( false ); }

  //! @brief Looks for the operation of type <a>top</a> with left and right operands <a>lop</a>, <a>rop</a> in set <a>_Ops</a> and adds it if not found
  FFOp* _insert_operation
    ( const typename FFOp::TYPE top, FFVar*lop, FFVar*rop=0 );

  //! @brief Looks for the binary operation of type <a>top</a> with left and right operands <a>Var1</a>, <a>Var2</a> in set <a>_Ops</a> and adds it if not found; also adds new auxiliary variable in set <a>_Vars</a> and update list of dependencies in both operands in <a>_Vars</a>
  static FFVar& _insert_binary_operation
    ( const typename FFOp::TYPE top, const FFDep&dep, const FFVar&Var1, const FFVar&Var2 );

  //! @brief Looks for the binary operation of type <a>top</a> with left and right operands <a>Cst1</a>, <a>Var2</a> in set <a>_Ops</a> and adds it if not found; also adds new auxiliary variable as well as constant <a>Cst1</a> in set <a>_Vars</a> and update list of dependencies in both operands in <a>_Vars</a>
  template <typename U> static FFVar&
  _insert_binary_operation
    ( const typename FFOp::TYPE top, const FFDep&dep, const U&Cst1, const FFVar&Var2 );

  //! @brief Looks for the binary operation of type <a>top</a> with left and right operands <a>Var1</a>, <a>Cst2</a> in set <a>_Ops</a> and adds it if not found; also adds new auxiliary variable as well as constant <a>Cst2</a> in set <a>_Vars</a> and update list of dependencies in both operands in <a>_Vars</a>
  template <typename U> static FFVar&
  _insert_binary_operation
    ( const typename FFOp::TYPE top, const FFDep&dep, const FFVar&Var1, const U&Cst2 );

  //! @brief Looks for the unary operation of type <a>top</a> with operand <a>Var1</a>, <a>Var2</a> in set <a>_Ops</a> and adds it if not found; also adds new auxiliary variable in set <a>_Vars</a> and update list of dependencies in both operands in <a>_Vars</a>
  static FFVar& _insert_unary_operation
    ( const typename FFOp::TYPE top, const FFDep&dep, const FFVar&Var );

  //! @brief Adds the auxiliary variable with dependency <a>dep</a> from operation <a>op</a>
  FFVar* _add_auxiliary
    ( const FFDep&dep, FFOp*pOp );

  //! @brief Looks for the real constant <a>x</a> and adds it if not found
  FFVar* _add_constant
    ( const double x );

  //! @brief Looks for the integer constant <a>n</a> and adds it if not found
  FFVar* _add_constant
    ( const int n );

  //! @brief Sets a variables to a constant and sets its numerical field
  FFVar* _set_constant
    ( const FFVar*pVar, const FFNum&num );

  //! @brief Unsets a variables to a constant
  FFVar* _unset_constant
    ( const FFVar*pVar );

  //! @brief Erase all variables in _Vars
  void _clear_variables()
    { it_Vars itv = _Vars.begin();
      for( ; itv != _Vars.end(); ++itv ) delete *itv;
      _Vars.clear(); }

  //! @brief Appends the auxiliary variable <a>pAux</a> and define it in _Ops with type <a>tOp</a>
  void _append_aux
    ( FFVar*pAux, typename FFOp::TYPE tOp );

  //! @brief Appends new auxiliary variable
  virtual void _append_aux
    ( FFVar*pAux );

  //! @brief Appends new original variable
  virtual void _append_var
    ( FFVar*pVar );

  //! @brief Search for the variable with identify <a>id</a> in <a>_Vars</a>
  FFVar* _find_var
    ( const typename FFVar::pt_idVar&id );

private:
  //! @brief Private methods to block default compiler methods
  FFGraph(const FFGraph&);
  FFGraph& operator=(const FFGraph&);

};

////////////////////////////////// FFNum ///////////////////////////////////////

inline std::ostream&
operator <<
( std::ostream&out, const FFNum&Num )
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

inline FFVar::FFVar
( FFGraph*dag )
: _dag( dag? dag: throw typename FFGraph::Exceptions( FFGraph::Exceptions::INIT )),
  _id( VAR, _dag->_nvar++ ), _num( 0./0. ), _dep(), _val( 0 ), _cst( false )
{ 
  // Initialize dependence
  _dep.indep(_id.second);

  // Insert new variable in set FFGraph::_Vars and corresponding operation in set FFGraph::_Ops
  FFVar* pVar = new FFVar( *this );
  FFOp* pOp = new FFOp( FFOp::VAR, 0, 0, pVar );
  _dag->_Ops.insert( pOp );
  pVar->_ops.first = _ops.first = pOp;
  _dag->_append_var( pVar );
}

inline FFVar::FFVar
( FFGraph*dag, const double d )
: _dag( dag? dag: throw typename FFGraph::Exceptions( FFGraph::Exceptions::INIT )),
  _id( VAR, _dag->_nvar++ ), _num( d ), _dep(), _val( 0 ), _cst( true )
{ 
  // Initialize dependence
  _dep.indep(_id.second);

  // Insert new variable in set FFGraph::_Vars and corresponding operation in set FFGraph::_Ops
  FFVar* pVar = new FFVar( *this );
  FFOp* pOp = new FFOp( FFOp::VAR, 0, 0, pVar );
  _dag->_Ops.insert( pOp );
  pVar->_ops.first = _ops.first = pOp;
  _dag->_append_var( pVar );
}

inline FFVar::FFVar
( FFGraph*dag, const int i )
: _dag( dag? dag: throw typename FFGraph::Exceptions( FFGraph::Exceptions::INIT )),
  _id( VAR, _dag->_nvar++ ), _num( i ), _dep(), _val( 0 ), _cst( true )
{ 
  // Initialize dependence
  _dep.indep(_id.second);

  // Insert new variable in set FFGraph::_Vars and corresponding operation in set FFGraph::_Ops
  FFVar* pVar = new FFVar( *this );
  FFOp* pOp = new FFOp( FFOp::VAR, 0, 0, pVar );
  _dag->_Ops.insert( pOp );
  pVar->_ops.first = _ops.first = pOp;
  _dag->_append_var( pVar );
}

inline FFVar::FFVar
( FFGraph*dag, const FFDep&dep, FFOp*op )
: _dag( dag ), _id( AUX, dag->_naux++ ), _num( 0./0. ), _dep(dep),
  _val ( 0 ), _cst( false )
{ _ops.first = op; }

inline FFVar::FFVar
( FFGraph*dag, const pt_idVar&id )
: _dag( dag ), _id( id ), _num( 0 ), _dep(), _val( 0 ), _cst( false )
{ _ops.first = 0; }

inline void FFVar::set
( const int i ) const
{ _num = i; _cst = true;
  if( !_dag ) return;
  _dag->_set_constant( this, _num );
  return;
}

inline void FFVar::set
( const double d ) const
{ _num = d; _cst = true;
  if( !_dag ) return;
  _dag->_set_constant( this, _num );
  return;
}

inline void FFVar::unset
() const
{ _cst = false;
  if( !_dag ) return;
  _dag->_unset_constant( this );
  return;
}

inline std::ostream&
operator <<
( std::ostream&out, const FFVar&Var)
{
  if( Var.id().second == FFVar::NOREF ) out << Var.num();
  else out << Var.name();
   // << " <= " << std::left << Var._num << "\t(" << Var._dag << ")";
  return out;
}

inline FFVar&
FFVar::operator=
( const FFVar&Var )
{
  if( this == &Var ) return *this;
  _id   = Var._id;
  _num  = Var._num;
  _dep  = Var._dep;
  _dag  = Var._dag;
  _val  = Var._val;
  _cst  = Var._cst;
  _ops  = Var._ops;
  return *this;
}

inline FFVar&
FFVar::operator=
( const int i )
{
  _id   = std::make_pair(CINT,NOREF);
  _num  = i;
  _dep  = i;
  _dag  = 0;
  _val  = 0;
  _cst  = true;
  _ops.first = 0;
  _ops.second.clear();
  return *this;
}

inline FFVar&
FFVar::operator=
( const double x )
{
  _id   = std::make_pair(CREAL,NOREF);
  _num  = x;
  _dep  = x;
  _dag  = 0;
  _val  = 0;
  _cst  = true;
  _ops.first = 0;
  _ops.second.clear();
  return *this;
}

inline FFVar
operator+
( const FFVar&Var )
{
  return Var;
}

template <typename U> inline FFVar&
FFVar::operator+=
( const U&Var )
{
  FFVar VarNew = *this + Var;
  *this = VarNew;
  return *this;
}

inline FFVar
operator+
( const FFVar&Var1, const FFVar&Var2 )
{ 
  if( &Var1 == &Var2 ) return( 2. * Var1 );

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
      return FFGraph::_insert_binary_operation( FFOp::PLUS, Var1._dep+Var2._dep, Var1, Var2 );
    }
  }
/*
  // Case either or both operands are (unreferenced) numeric constants
  if( Var1._id.second == FFVar::NOREF && Var2._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.n + Var2._num.n );
        case FFNum::REAL:  return( Var1._num.n + Var2._num.x );
      }
      case FFNum::REAL:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.x + Var2._num.n );
        case FFNum::REAL:  return( Var1._num.x + Var2._num.x );
      }
    }
  }
  if( Var1._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:   return( Var2 + Var1._num.n );
      case FFNum::REAL:  return( Var2 + Var1._num.x );
    }
  }
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Var1 + Var2._num.n );
      case FFNum::REAL:  return( Var1 + Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFGraph::_insert_binary_operation( FFOp::PLUS, Var1._dep+Var2._dep, Var1, Var2 );
*/
}

template <typename U> inline FFVar
operator+
( const U&Cst1, const FFVar&Var2 )
{
  // Case constant is zero
  if( Cst1 == U(0) ) return Var2;

  switch( Var2._id.first ){
  case FFVar::CREAL:
    return( Cst1 + Var2._num.x );
  case FFVar::CINT:
    return( Cst1 + Var2._num.n );
  default:
    return FFGraph::_insert_binary_operation( FFOp::PLUS, Cst1+Var2._dep, (double)Cst1, Var2 );
  }
/*
  // Case right operand is (unreferenced) numeric constrants
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Cst1 + Var2._num.n );
      case FFNum::REAL:  return( Cst1 + Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant Cst1 if not defined
  return FFGraph::_insert_binary_operation( FFOp::PLUS, Cst1+Var2._dep, (double)Cst1, Var2 );
*/
}

template <typename U> inline FFVar
operator+
( const FFVar&Var1, const U&Cst2 )
{
  return( Cst2 + Var1 );
}

inline FFVar
operator-
( const FFVar&Var )
{

  // Check if expression of type -(-X)
  if( Var._ops.first && Var._ops.first->type == FFOp::NEG )
    return *Var._ops.first->plop;

  switch( Var._id.first ){
  case FFVar::CREAL:
    return( -Var._num.x );
  case FFVar::CINT:
    return( -Var._num.n );
  default:
    return FFGraph::_insert_unary_operation( FFOp::NEG, -Var._dep, Var );
  }
/*
  // Case right operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( -Var._num.n );
      case FFNum::REAL:  return( -Var._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFGraph::_insert_unary_operation( FFOp::NEG, -Var._dep, Var );
*/
}

template <typename U> inline FFVar&
FFVar::operator-=
( const U&Var )
{
  FFVar VarNew = *this - Var;
  *this = VarNew;
  return *this;
}

inline FFVar
operator-
( const FFVar&Var1, const FFVar&Var2 )
{
  if( &Var1 == &Var2 ) return 0.;

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
      return FFGraph::_insert_binary_operation( FFOp::MINUS, Var1._dep-Var2._dep, Var1, Var2 );
    }
  }
/*
  // Case either or both operands are numeric constants
  if( Var1._id.second == FFVar::NOREF
   && Var2._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.n - Var2._num.n );
        case FFNum::REAL:  return( Var1._num.n - Var2._num.x );
      }
      case FFNum::REAL:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.x - Var2._num.n );
        case FFNum::REAL:  return( Var1._num.x - Var2._num.x );
      }
    }
  }

  if( Var1._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:   return( (double)Var1._num.n - Var2 );
      case FFNum::REAL:  return( Var1._num.x - Var2 );
    }
  }
  
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Var1 - (double)Var2._num.n );
      case FFNum::REAL:  return( Var1 - Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFGraph::_insert_binary_operation( FFOp::MINUS, Var1._dep-Var2._dep, Var1, Var2 );
*/
}

template <typename U> inline FFVar
operator-
( const FFVar&Var1, const U&Cst2 )
{
  // Case constant is zero
  if( Cst2 == U(0) ) return Var1;

  switch( Var1._id.first ){
  case FFVar::CREAL:
    return( Var1._num.x - Cst2 );
  case FFVar::CINT:
    return( Var1._num.n - Cst2 );
  default:
    return FFGraph::_insert_binary_operation( FFOp::MINUS, Var1._dep-Cst2, Var1, (double)Cst2 );
  }
/*
  // Case constant is zero
  if( Cst2 == 0. ) return Var1;

  // Case right operand is a numeric constant
  if( Var1._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:   return( Var1._num.n - Cst2 );
      case FFNum::REAL:  return( Var1._num.x - Cst2 );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant Cst2 if not defined
  return FFGraph::_insert_binary_operation( FFOp::MINUS, Var1._dep-Cst2, Var1, (double)Cst2 );
*/
}

template <typename U> inline FFVar
operator-
( const U&Cst1, const FFVar&Var2 )
{
  // Case constant is zero
  if( Cst1 == U(0) ) return -Var2;

  switch( Var2._id.first ){
  case FFVar::CREAL:
    return( Cst1 - Var2._num.x );
  case FFVar::CINT:
    return( Cst1 - Var2._num.n );
  default:
    return FFGraph::_insert_binary_operation( FFOp::MINUS, Cst1-Var2._dep, (double)Cst1, Var2 );
  }
/*
  // Case right operand is a numeric constant
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Cst1 - Var2._num.n );
      case FFNum::REAL:  return( Cst1 - Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant Cst1 if not defined
  return FFGraph::_insert_binary_operation( FFOp::MINUS, Cst1-Var2._dep, (double)Cst1, Var2 );
*/
}

template <typename U> inline FFVar&
FFVar::operator*=
( const U&Var )
{
  FFVar VarNew = *this * Var;
  *this = VarNew;
  return *this;
}

inline FFVar
operator*
( const FFVar&Var1, const FFVar&Var2 )
{
  if( &Var1 == &Var2 ) return sqr(Var1);

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
      return FFGraph::_insert_binary_operation( FFOp::TIMES, Var1._dep*Var2._dep, Var1, Var2 );
      //FFVar VarR = FFGraph::_insert_binary_operation( FFOp::TIMES, Var1._dep*Var2._dep, Var1, Var2 );
      //return !Var1._dag->options.DCDECOMPOSE? VarR: inter( VarR, (sqr(Var1+Var2)-sqr(Var1-Var2))/4. );
     }
    }
  }
/*
  // Case either or both operands are numeric constants
  if( Var1._id.second == FFVar::NOREF && Var2._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.n * Var2._num.n );
        case FFNum::REAL:  return( Var1._num.n * Var2._num.x );
      }
      case FFNum::REAL:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.x * Var2._num.n );
        case FFNum::REAL:  return( Var1._num.x * Var2._num.x );
      }
    }
  }

  if( Var1._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:   return( (double)Var1._num.n * Var2 );
      case FFNum::REAL:  return( Var1._num.x * Var2 );
    }
  }
  
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Var1 * (double)Var2._num.n );
      case FFNum::REAL:  return( Var1 * Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFGraph::_insert_binary_operation( FFOp::TIMES, Var1._dep*Var2._dep, Var1, Var2 );
*/
}

template <typename U> inline FFVar
operator*
( const FFVar&Var1, const U&Cst2 )
{
  return( Cst2 * Var1 );
}

template <typename U> inline FFVar
operator*
( const U&Cst1, const FFVar&Var2 )
{
  // Case constant is zero
  if( Cst1 == U(0) ) return 0.;
  // Case constant is one
  if( Cst1 == U(1) ) return Var2;
  // Case constant is one
  if( Cst1 == U(-1) ) return -Var2;

  switch( Var2._id.first ){
  case FFVar::CREAL:
    return( Cst1 * Var2._num.x );
  case FFVar::CINT:
    return( Cst1 * Var2._num.n );
  default:
    return FFGraph::_insert_binary_operation( FFOp::SCALE, Cst1*Var2._dep, (double)Cst1, Var2 );
  }
/*
  // Case right operand is a numeric constant
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Cst1 * Var2._num.n );
      case FFNum::REAL:  return( Cst1 * Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant Cst1 if not defined
  return FFGraph::_insert_binary_operation( FFOp::SCALE, Cst1*Var2._dep, (double)Cst1, Var2 );
*/
}

template <typename U> inline FFVar&
FFVar::operator/=
( const U&Var )
{
  FFVar VarNew = *this / Var;
  *this = VarNew;
  return *this;
}

inline FFVar
operator/
( const FFVar&Var1, const FFVar&Var2 )
{
  if( &Var1 == &Var2 ) return 1.;

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
      return FFGraph::_insert_binary_operation( FFOp::DIV, Var1._dep/Var2._dep, Var1, Var2 );
      //FFVar VarR = FFGraph::_insert_binary_operation( FFOp::DIV, Var1._dep/Var2._dep, Var1, Var2 );
      //return !Var1._dag->options.DCDECOMPOSE? VarR: inter( VarR, Var1*pow(Var2,-1) );
     }
    }
  }
/*
  // Case either or both operands are numeric constants
  if( Var1._id.second == FFVar::NOREF && Var2._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.n / Var2._num.n );
        case FFNum::REAL:  return( Var1._num.n / Var2._num.x );
      }
      case FFNum::REAL:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.x / Var2._num.n );
        case FFNum::REAL:  return( Var1._num.x / Var2._num.x );
      }
    }
  }

  if( Var1._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:   return( Var1._num.n / Var2 );
      case FFNum::REAL:  return( Var1._num.x / Var2 );
    }
  }
  
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Var1 / Var2._num.n );
      case FFNum::REAL:  return( Var1 / Var2._num.x );
    }
  }
  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFGraph::_insert_binary_operation( FFOp::DIV, Var1._dep/Var2._dep, Var1, Var2 );
*/
}

template <typename U> inline FFVar
operator/
( const FFVar&Var1, const U&Cst2 )
{
  return( ( 1. / Cst2 ) * Var1 );
}

template <typename U> inline FFVar
operator/
( const U&Cst1, const FFVar&Var2 )
{
  // Case constant is zero
  if( Cst1 == U(0) ) return 0.;

  switch( Var2._id.first ){
  case FFVar::CREAL:
    return( Cst1 / Var2._num.x );
  case FFVar::CINT:
    return( Cst1 / Var2._num.n );
  default:{
    //FFVar VarR = FFGraph::_insert_binary_operation( FFOp::DIV, Cst1/Var2._dep, (double)Cst1, Var2 );
    //return !Var2._dag->options.DCDECOMPOSE? VarR: inter( VarR, Cst1*pow(Var2,-1) );
    return FFGraph::_insert_binary_operation( FFOp::DIV, Cst1/Var2._dep, (double)Cst1, Var2 );
   }
  }
/*
  // Case right operand is a numeric constant
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Cst1 / Var2._num.n );
      case FFNum::REAL:  return( Cst1 / Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant Cst1 if not defined
  return FFGraph::_insert_binary_operation( FFOp::DIV, Cst1/Var2._dep, (double)Cst1, Var2 );
*/
}

inline FFVar
inv
( const FFVar&Var )
{
  return( pow( Var, -1 ) );
}

inline FFVar
max
( const FFVar&Var1, const FFVar&Var2 )
{
  if( &Var1 == &Var2 ) return Var1;

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
  return FFGraph::_insert_binary_operation( FFOp::MAXF, max(Var1._dep,Var2._dep), Var1, Var2 );
}

template <typename U> inline FFVar
max
( const FFVar&Var1, const U&Cst2 )
{
  return( max( Cst2, Var1 ) );
}

template <typename U> inline FFVar
max
( const U&Cst1, const FFVar&Var2 )
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
  return FFGraph::_insert_binary_operation( FFOp::MAXF, max((double)Cst1,Var2._dep), (double)Cst1, Var2 );
}

inline FFVar
max
( const unsigned nVar, const FFVar*pVar )
{
  if( nVar<=0 || !pVar ) return( 0 );
  
  FFVar VarR = pVar[0];
  for( unsigned int i=1; i<nVar; i++ ) VarR = max( VarR, pVar[i] );
  return( VarR );
}

inline FFVar
min
( const FFVar&Var1, const FFVar&Var2 )
{
  if( &Var1 == &Var2 ) return Var1;

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
  return FFGraph::_insert_binary_operation( FFOp::MINF, min(Var1._dep,Var2._dep), Var1, Var2 );
}

template <typename U> inline FFVar
min
( const FFVar&Var1, const U&Cst2 )
{
  return( min( Cst2, Var1 ) );
}

template <typename U> inline FFVar
min
( const U&Cst1, const FFVar&Var2 )
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
  return FFGraph::_insert_binary_operation( FFOp::MINF, min((double)Cst1,Var2._dep), (double)Cst1, Var2 );
}

inline FFVar
min
( const unsigned nVar, const FFVar*pVar )
{
  if( nVar<=0 || !pVar ) return( 0 );
  
  FFVar VarR = pVar[0];
  for( unsigned int i=1; i<nVar; i++ ) VarR = min( VarR, pVar[i] );
  return( VarR );
}

inline FFVar
inter
( const FFVar&Var1, const FFVar&Var2 )
{
  if( &Var1 == &Var2 ) return Var1;

  // Case either or both operands are numeric constants
  if( Var1._id.second == FFVar::NOREF && Var2._id.second == FFVar::NOREF ){
    if( Var1._num.val() != Var2._num.val() )
      throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTER ); 
    return Var1._num.val();
  }
  if( Var1._id.second == FFVar::NOREF )
    return inter( Var2, Var1._num.val() );
  if( Var2._id.second == FFVar::NOREF )
    return inter( Var1, Var2._num.val() );

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFGraph::_insert_binary_operation( FFOp::INTER, inter(Var1._dep,Var2._dep), Var1, Var2 );
}

template <typename U> inline FFVar
inter
( const U&Cst1, const FFVar&Var2 )
{
  // Case right operand is a numeric constant
  if( Var2._id.second == FFVar::NOREF ){
    if( Cst1 != Var2._num.val() )
      throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTER ); 
    return Cst1;
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFGraph::_insert_binary_operation( FFOp::INTER, inter((double)Cst1,Var2._dep), (double)Cst1, Var2 );
}

template <typename U> inline FFVar
inter
( const FFVar&Var1, const U&Cst2 )
{
  return( inter( Cst2, Var1 ) );
}

inline FFVar
exp
( const FFVar&Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::exp( Var._num.n ) );
      case FFNum::REAL:  return( std::exp( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFGraph::_insert_unary_operation( FFOp::EXP, exp(Var._dep), Var );
}

inline FFVar
log
( const FFVar&Var )
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
  return FFGraph::_insert_unary_operation( FFOp::LOG, log(Var._dep), Var );
}

inline FFVar
sqr
( const FFVar&Var )
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
  return FFGraph::_insert_unary_operation( FFOp::SQR, sqr(Var._dep), Var );
}

inline FFVar
sqrt
( const FFVar&Var )
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
  return FFGraph::_insert_unary_operation( FFOp::SQRT, sqrt(Var._dep), Var );
}

inline FFVar
fabs
( const FFVar&Var )
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
  return FFGraph::_insert_unary_operation( FFOp::FABS, fabs(Var._dep), Var );
}

inline FFVar
pow
( const FFVar&Var, const int iExp )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::pow( Var._num.n, iExp ) );
      case FFNum::REAL:  return( std::pow( Var._num.x, iExp ) );
    }
  }

  // Case integer exponent is 0,1, or 2
  if( iExp == 0 ) return( 1. );
  if( iExp == 1 ) return Var;
  if( iExp == 2 ) return sqr(Var);

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant iExp if not defined
  return FFGraph::_insert_binary_operation( FFOp::IPOW, pow(Var._dep,iExp), Var, iExp );
}

inline FFVar
cheb
( const FFVar&Var, const unsigned iOrd )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( mc::cheb( Var._num.n, iOrd ) );
      case FFNum::REAL:  return( mc::cheb( Var._num.x, iOrd ) );
    }
  }

  // Case integer exponent is 0 or 1
  if( iOrd == 0 ) return( 1. );
  if( iOrd == 1 ) return Var;

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant iOrd if not defined
  FFVar VarRecu = 2.*Var*cheb(Var,iOrd-1) - cheb(Var,iOrd-2);
  FFVar VarCheb = FFGraph::_insert_binary_operation( FFOp::CHEB, cheb(Var._dep,iOrd), Var, (int)iOrd );
  return inter( VarRecu, VarCheb );
}

inline FFVar
pow
( const FFVar&Var1, const double Cst2 )
{
  return exp( Cst2 * log( Var1 ) );
}

inline FFVar
pow
( const FFVar&Var1, const FFVar&Var2 )
{
  // Case exponent is integer constant
  if( Var2._id.second == FFVar::NOREF && Var2._num.t == FFNum::INT )
    return pow( Var1, Var2._num.n );

  return exp( Var2 * log( Var1 ) );
}

inline FFVar
pow
( const double Cst1, const FFVar&Var2 )
{
  // Case exponent is integer constant
  if( Var2._id.second == FFVar::NOREF && Var2._num.t == FFNum::INT )
    return std::pow( Cst1, Var2._num.n );

  return exp( Var2 * std::log( Cst1 ) );
}

inline FFVar
cos
( const FFVar&Var )
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
  return FFGraph::_insert_unary_operation( FFOp::COS, cos(Var._dep), Var );
}

inline FFVar
sin
( const FFVar&Var )
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
  return FFGraph::_insert_unary_operation( FFOp::SIN, sin(Var._dep), Var );
}

inline FFVar
tan
( const FFVar&Var )
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
  return FFGraph::_insert_unary_operation( FFOp::TAN, tan(Var._dep), Var );
}

inline FFVar
asin
( const FFVar&Var )
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
  return FFGraph::_insert_unary_operation( FFOp::ASIN, asin(Var._dep), Var );
}

inline FFVar
acos
( const FFVar&Var )
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
  return FFGraph::_insert_unary_operation( FFOp::ACOS, acos(Var._dep), Var );
}

inline FFVar
atan
( const FFVar&Var )
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
  return FFGraph::_insert_unary_operation( FFOp::ATAN, atan(Var._dep), Var );
}

inline FFVar
erf
( const FFVar&Var )
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
  return FFGraph::_insert_unary_operation( FFOp::ERF, erf(Var._dep), Var );
}

inline FFVar
erfc
( const FFVar&Var )
{
  return ( 1. - erf( Var ) );
}

inline FFVar
fstep
( const FFVar&Var )
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
  return FFGraph::_insert_unary_operation( FFOp::FSTEP, fstep(Var._dep), Var );
}

inline FFVar
bstep
( const FFVar&Var )
{
  return ( fstep( -Var ) );
}


////////////////////////////////// FFOp ////////////////////////////////////////

inline
FFOp::FFOp
( TYPE top, FFVar*lop, FFVar*rop, FFVar*res ):
  type( top ), pres( res ), _visited(false)
{
  // Order operands for commutative binary operators
  if( !lop
    || ( top != PLUS && top != TIMES && top != SCALE )
    || lt_FFVar()( lop, rop ) )
    { plop = lop; prop = rop; }
  else
    { plop = rop; prop = lop; }
}

inline std::ostream&
operator <<
( std::ostream&out, const FFOp&Op)
{
  switch( Op.type ){
    case FFOp::CNST:  out << Op.pres->num() << "\t"; break;
    case FFOp::VAR:   out << "VARIABLE"; break;
    case FFOp::PLUS:  out << FFVar::_name( Op.plop->id() ) << " + " << FFVar::_name( Op.prop->id() ) << "\t"; break;
    case FFOp::NEG:   out << "- " << FFVar::_name( Op.plop->id() ) << "\t" ; break;
    case FFOp::MINUS: out << FFVar::_name( Op.plop->id() ) << " - " << FFVar::_name( Op.prop->id() ) << "\t"; break;
    case FFOp::TIMES:
    case FFOp::SCALE: out << FFVar::_name( Op.plop->id() ) << " * " << FFVar::_name( Op.prop->id() ) << "\t"; break;
    case FFOp::DIV:   out << FFVar::_name( Op.plop->id() ) << " / " << FFVar::_name( Op.prop->id() ) << "\t"; break;
    case FFOp::MINF:  out << "MIN( " << FFVar::_name( Op.plop->id() ) << ", " << FFVar::_name( Op.prop->id() ) << " )"; break;
    case FFOp::MAXF:  out << "MAX( " << FFVar::_name( Op.plop->id() ) << ", " << FFVar::_name( Op.prop->id() ) << " )"; break;
    case FFOp::INTER: out << "INTER( " << FFVar::_name( Op.plop->id() ) << ", " << FFVar::_name( Op.prop->id() ) << " )"; break;
    case FFOp::EXP:   out << "EXP( " << FFVar::_name( Op.plop->id() ) << " )\t"; break;
    case FFOp::LOG:   out << "LOG( " << FFVar::_name( Op.plop->id() ) << " )\t"; break;
    case FFOp::SQR:   out << "SQR( " << FFVar::_name( Op.plop->id() ) << " )\t"; break;
    case FFOp::SQRT:  out << "SQRT( " << FFVar::_name( Op.plop->id() ) << " )\t"; break;
    case FFOp::FABS:  out << "FABS( " << FFVar::_name( Op.plop->id() ) << " )\t"; break;
    case FFOp::IPOW:  out << "POW( " << FFVar::_name( Op.plop->id() ) << ", " << FFVar::_name( Op.prop->id() ) << " )"; break;
    case FFOp::CHEB:  out << "CHEB( " << FFVar::_name( Op.plop->id() ) << ", " << FFVar::_name( Op.prop->id() ) << " )"; break;
    case FFOp::COS:   out << "COS( " << FFVar::_name( Op.plop->id() ) << " )\n"; break;
    case FFOp::SIN:   out << "SIN( " << FFVar::_name( Op.plop->id() ) << " )\n"; break;
    case FFOp::TAN:   out << "TAN( " << FFVar::_name( Op.plop->id() ) << " )\n"; break;
    case FFOp::ASIN:  out << "ASIN( " << FFVar::_name( Op.plop->id() ) << " )\n"; break;
    case FFOp::ACOS:  out << "ACOS( " << FFVar::_name( Op.plop->id() ) << " )\n"; break;
    case FFOp::ATAN:  out << "ATAN( " << FFVar::_name( Op.plop->id() ) << " )\n"; break;
    case FFOp::ERF:   out << "ERF( " << FFVar::_name( Op.plop->id() ) << " )\n"; break;
    case FFOp::FSTEP: out << "FSTEP( " << FFVar::_name( Op.plop->id() ) << " )\n"; break;
    default:;
  } 
  return out;
}

inline void
FFOp::propagate_subgraph
( std::list<const FFOp*>&Ops ) const
{
  if( _visited ) return;
  _visited = true;

  //Ops.push_front( this );
  if( plop && plop->ops().first ) plop->ops().first->propagate_subgraph( Ops );
  if( prop && prop->ops().first ) prop->ops().first->propagate_subgraph( Ops );
  Ops.push_back( this );
}

template <typename U> inline void
FFOp::reset_val_subgraph
( const U& U_dum ) const
{
  if( _visited ) return;
  _visited = true;

  if( plop && plop->ops().first ) plop->ops().first->reset_val_subgraph( U_dum );
  if( prop && prop->ops().first ) prop->ops().first->reset_val_subgraph( U_dum );
  if( pres && pres->val() ) pres->reset_val( U_dum );
}

template <typename U> inline void
FFOp::reset_val_subgraph
( const U& U_dum, const std::vector<const FFVar*>&vDep,
  const std::vector<const FFVar*>&vIndep ) const
{
  if( _visited ) return;
  _visited = true;

  if( plop && plop->ops().first ) plop->ops().first->reset_val_subgraph( U_dum, vDep, vIndep );
  if( prop && prop->ops().first ) prop->ops().first->reset_val_subgraph( U_dum, vDep, vIndep );
  if( pres && pres->val() ){
    // Do not reset _val field of independent variables
    typename std::vector<const FFVar*>::const_iterator iti = vIndep.begin();
    for( ; iti!=vIndep.end(); ++iti ) if( pres->id() == (*iti)->id() ) return;
    // Do not reset _val field of dependent variables
    typename std::vector<const FFVar*>::const_iterator itd = vDep.begin();
    for( ; itd!=vDep.end(); ++itd ) if( pres->id() == (*itd)->id() ) return;
    pres->reset_val( U_dum );
  }
}

template <typename U> inline void
FFOp::evaluate
( const U& U_dum ) const
{
  switch( type ){
   case FFOp::VAR:
    if( !pres->cst() ){ // do not override constant value if set
    //if( !pres->cst() || pres->val() ){ // override constant value if set
      if( !pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::MISSVAR );
      return;
    }
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::CONST );

   case FFOp::CNST:
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    switch( pres->num().t ){
      case FFNum::INT:  pres->val() = new U( pres->num().n ); return;
      case FFNum::REAL: pres->val() = new U( pres->num().x ); return;
    }

   case FFOp::PLUS:
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( *static_cast<U*>( plop->val() )
                       + *static_cast<U*>( prop->val() ) );
#ifdef MC__FFUNC_DEBUG_EVAL
    std::cout << "FFOp::PLUS:" << std::endl;
    std::cout << "plop: " << plop << "   plop->val(): " << plop->val() << std::endl;
    std::cout << "prop: " << prop << "   prop->val(): " << prop->val() << std::endl;
    std::cout << "pres: " << pres << "   pres->val(): " << pres->val() << std::endl;
#endif
    return;

   case FFOp::NEG:
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( - *static_cast<U*>( plop->val() ) );
    return;

   case FFOp::MINUS:
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( *static_cast<U*>( plop->val() )
                       - *static_cast<U*>( prop->val() ) );
    return;

   case FFOp::SCALE:
   case FFOp::TIMES:
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( *static_cast<U*>( plop->val() )
                       * *static_cast<U*>( prop->val() ) );
#ifdef MC__FFUNC_DEBUG_EVAL
    std::cout << "FFOp::TIMES:" << std::endl;
    std::cout << "plop: " << plop << "   plop->val(): " << plop->val() << std::endl;
    std::cout << "prop: " << prop << "   prop->val(): " << prop->val() << std::endl;
    std::cout << "pres: " << pres << "   pres->val(): " << pres->val() << std::endl;
#endif
    return;

   case FFOp::DIV:  
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( *static_cast<U*>( plop->val() )
                       / *static_cast<U*>( prop->val() ) );
    return;

   case FFOp::IPOW:
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::pow( *static_cast<U*>( plop->val() ),
                                     prop->num().n ) );
    return;

   case FFOp::CHEB:
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::cheb( *static_cast<U*>( plop->val() ),
                                      prop->num().n ) );
    return;

   case FFOp::EXP:  
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::exp( *static_cast<U*>( plop->val() ) ) );
    return;

   case FFOp::LOG:  
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::log( *static_cast<U*>( plop->val() ) ) );
    return;

   case FFOp::SQR:  
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::sqr( *static_cast<U*>( plop->val() ) ) );
#ifdef MC__FFUNC_DEBUG_EVAL
    std::cout << "FFOp::SQR:" << std::endl;
    std::cout << "plop: " << plop << "   plop->val(): " << plop->val() << std::endl;
    std::cout << "pres: " << pres << "   pres->val(): " << pres->val() << std::endl;
#endif
    return;

   case FFOp::SQRT: 
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::sqrt( *static_cast<U*>( plop->val() ) ) );
    return;

   case FFOp::COS:  
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::cos( *static_cast<U*>( plop->val() ) ) );
    return;

   case FFOp::SIN:  
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::sin( *static_cast<U*>( plop->val() ) ) );
    return;

   case FFOp::TAN:  
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::tan( *static_cast<U*>( plop->val() ) ) );
    return;

   case FFOp::ACOS: 
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::acos( *static_cast<U*>( plop->val() ) ) );
    return;

   case FFOp::ASIN: 
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::asin( *static_cast<U*>( plop->val() ) ) );
    return;

   case FFOp::ATAN: 
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::atan( *static_cast<U*>( plop->val() ) ) );
    return;

   case FFOp::ERF:  
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::erf( *static_cast<U*>( plop->val() ) ) );
    return;

   case FFOp::FABS: 
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::fabs( *static_cast<U*>( plop->val() ) ) );
    return;

   case FFOp::FSTEP:
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::fstep( *static_cast<U*>( plop->val() ) ) );
    return;

   case FFOp::MINF: 
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::min( *static_cast<U*>( plop->val() ),
                                     *static_cast<U*>( prop->val() ) ) );
    return;

   case FFOp::MAXF: 
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U( Op<U>::max( *static_cast<U*>( plop->val() ),
                                     *static_cast<U*>( prop->val() ) ) );
    return;

   case FFOp::INTER: 
    //if( pres->val() ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
    pres->val() = new U;
    if( !Op<U>::inter( *static_cast<U*>( pres->val() ), *static_cast<U*>( plop->val() ),
                       *static_cast<U*>( prop->val() ) ) )
      throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTER );
    return;

   default:
     throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
  }
}

template <typename U> inline void
FFOp::evaluate
( U* pUres ) const
{
  switch( type ){
   case FFOp::VAR:
    if( !pres->cst() ) break; // do not override constant value if set
    //if( !pres->cst() || pUres ) break; // override constant value if set
    //if( pUres ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::CONST );

   case FFOp::CNST:
    switch( pres->num().t ){
      case FFNum::INT:  *pUres = pres->num().n; break;
      case FFNum::REAL: *pUres = pres->num().x; break;
    }
    break;

   case FFOp::PLUS:
    *pUres = *static_cast<U*>( plop->val() ) + *static_cast<U*>( prop->val() );
    break;

   case FFOp::NEG:
    *pUres = - *static_cast<U*>( plop->val() );
    break;

   case FFOp::MINUS:
    *pUres = *static_cast<U*>( plop->val() ) - *static_cast<U*>( prop->val() );
    break;

   case FFOp::SCALE:
   case FFOp::TIMES:
    *pUres = *static_cast<U*>( plop->val() ) * *static_cast<U*>( prop->val() );
    break;

   case FFOp::DIV:  
    *pUres = *static_cast<U*>( plop->val() ) / *static_cast<U*>( prop->val() );
    break;

   case FFOp::IPOW:
    *pUres = Op<U>::pow( *static_cast<U*>( plop->val() ), prop->num().n );
    break;

   case FFOp::CHEB:
    *pUres = Op<U>::cheb( *static_cast<U*>( plop->val() ), prop->num().n );
    break;

   case FFOp::EXP:  
    *pUres = Op<U>::exp( *static_cast<U*>( plop->val() ) );
    break;

   case FFOp::LOG:  
    *pUres = Op<U>::log( *static_cast<U*>( plop->val() ) );
    break;

   case FFOp::SQR:  
    *pUres = Op<U>::sqr( *static_cast<U*>( plop->val() ) );
    break;

   case FFOp::SQRT: 
    *pUres = Op<U>::sqrt( *static_cast<U*>( plop->val() ) );
    break;

   case FFOp::COS:  
    *pUres = Op<U>::cos( *static_cast<U*>( plop->val() ) );
    break;

   case FFOp::SIN:  
    *pUres = Op<U>::sin( *static_cast<U*>( plop->val() ) );
    break;

   case FFOp::TAN:  
    *pUres = Op<U>::tan( *static_cast<U*>( plop->val() ) );
    break;

   case FFOp::ACOS: 
    *pUres = Op<U>::acos( *static_cast<U*>( plop->val() ) );
    break;

   case FFOp::ASIN: 
    *pUres = Op<U>::asin( *static_cast<U*>( plop->val() ) );
    break;

   case FFOp::ATAN: 
    *pUres = Op<U>::atan( *static_cast<U*>( plop->val() ) );
    break;

   case FFOp::ERF:  
    *pUres = Op<U>::erf( *static_cast<U*>( plop->val() ) );
    break;

   case FFOp::FABS: 
    *pUres = Op<U>::fabs( *static_cast<U*>( plop->val() ) );
    break;

   case FFOp::FSTEP:
    *pUres = Op<U>::fstep( *static_cast<U*>( plop->val() ) );
    break;

   case FFOp::MINF: 
    *pUres = Op<U>::min( *static_cast<U*>( plop->val() ), *static_cast<U*>( prop->val() ) );
    break;

   case FFOp::MAXF: 
    *pUres = Op<U>::max( *static_cast<U*>( plop->val() ), *static_cast<U*>( prop->val() ) );
    break;

   case FFOp::INTER: 
    if( !Op<U>::inter( *pUres, *static_cast<U*>( plop->val() ), *static_cast<U*>( prop->val() ) ) )
      throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTER );
    break;

   default:
     throw typename FFGraph::Exceptions( FFGraph::Exceptions::INTERN );
  }

  pres->val() = pUres;
  //std::cout << "evaluation of " << *pres << ":\n" << *pUres;
  return;
}

inline void
FFOp::generate_dot_script
( std::ostream&os ) const
{
  if( _visited ) return;
  _visited = true;

  if( plop && plop->ops().first ) plop->ops().first->generate_dot_script( os );
  if( prop && prop->ops().first ) prop->ops().first->generate_dot_script( os );
  append_dot_script( os );
}

inline void
FFOp::append_dot_script
( std::ostream&os ) const
{
  std::ostringstream var_color; var_color << "red";
  std::ostringstream aux_color; aux_color << "blue";
  std::ostringstream op_color;  op_color  << "black";
  switch( type ){
  case FFOp::VAR:   return append_dot_script_variable( os, pres->cst()?true:false, 14 );
  case FFOp::CNST:  return append_dot_script_variable( os, true,  14 );
  case FFOp::PLUS:  return append_dot_script_factor( os, " + ",   false, 18 );
  case FFOp::NEG:   return append_dot_script_factor( os, " - ",   true,  18 );
  case FFOp::MINUS: return append_dot_script_factor( os, " - ",   false, 18 );
  case FFOp::SCALE: return append_dot_script_factor( os, " x ",   false, 18 );
  case FFOp::TIMES: return append_dot_script_factor( os, " x ",   false, 18 );
  case FFOp::DIV:   return append_dot_script_factor( os, " / ",   false, 18 );
  case FFOp::MINF:  return append_dot_script_factor( os, "min",   false, 14 );
  case FFOp::MAXF:  return append_dot_script_factor( os, "max",   false, 14 );
  case FFOp::INTER: return append_dot_script_factor( os, "inter", false, 14 );
  case FFOp::IPOW:  return append_dot_script_factor( os, "pow",   false, 14, true );
  case FFOp::CHEB:  return append_dot_script_factor( os, "cheb",  false, 14, true );
  case FFOp::EXP:   return append_dot_script_factor( os, "exp",   true,  14 );
  case FFOp::LOG:   return append_dot_script_factor( os, "log",   true,  14 );
  case FFOp::FABS:  return append_dot_script_factor( os, "fabs",  true,  14 );
  case FFOp::SQR:   return append_dot_script_factor( os, "sqr",   true,  14 );
  case FFOp::SQRT:  return append_dot_script_factor( os, "sqrt",  true,  14 );
  case FFOp::COS:   return append_dot_script_factor( os, "cos",   true,  14 );
  case FFOp::SIN:   return append_dot_script_factor( os, "sin",   true,  14 );
  case FFOp::TAN:   return append_dot_script_factor( os, "tan",   true,  14 );
  case FFOp::ACOS:  return append_dot_script_factor( os, "acos",  true,  14 );
  case FFOp::ASIN:  return append_dot_script_factor( os, "asin",  true,  14 );
  case FFOp::ATAN:  return append_dot_script_factor( os, "atan",  true,  14 );
  case FFOp::ERF:   return append_dot_script_factor( os, "erf",   true,  14 );
  case FFOp::FSTEP: return append_dot_script_factor( os, "fstep", true,  14 );
  default: os << "/* a factor was not displayed */\n";
  }
  return;
}

inline void
FFOp::append_dot_script_factor
( std::ostream&os, const std::string&fname, const bool unary,
  const unsigned int fontsize, const bool dotted ) const
{
  std::ostringstream op_color; op_color  << "black";

  os << "  " << pres->name() << " [shape=record,fontname=\"Arial\",color="
     << op_color.str().c_str() << ",label=\"<f0> " << fname.c_str() << "|<f1> "
     << pres->name() << "\"];\n";
  //os << "  \"" << plop->name() << "\":f0 -> \"" << pres->name() << "\":f0;\n";
  os << "  " << plop->name() << " -> " << pres->name() << " [arrowsize=0.7];\n";
  if( unary ) return;
  os << "  " << prop->name() << " -> " << pres->name() << " [arrowsize=0.7"
     << (dotted? ",style=dotted];\n": "];\n");
  //os << "  \"" << prop->name() << "\":f0 -> \"" << pres->name() << "\":f0 "
  //   << (dotted? "[style=dotted];\n": ";\n");
}

inline void
FFOp::append_dot_script_variable
( std::ostream&os, const bool constant, const unsigned int fontsize ) const
{
  std::ostringstream var_color; var_color << "red";
  std::ostringstream cst_color; cst_color << "blue";

  if( constant )
    os << "  " << pres->name() << " [shape=record,fontname=\"Arial\",color="
     << cst_color.str().c_str() << ",label=\"<f0> " << pres->num() << "|<f1> "
     << pres->name() << "\"];\n"; 
  else
    os << "  " << pres->name() << " [shape=ellipse,fontname=\"Arial\",color="
       << (var_color.str().c_str()) << "];\n";
    //   << (var_color.str().c_str()) << ",label=\"<f0> " << pres->name() << "\"];\n";
}

inline bool
FFOp::is_univariate
() const
{
  switch( type ){
  case FFOp::PLUS: case FFOp::MINUS: case FFOp::TIMES: case FFOp::DIV:
  case FFOp::MINF: case FFOp::MAXF: case FFOp::INTER:
    return false;
  case FFOp::NEG:  case FFOp::SCALE: case FFOp::IPOW: case FFOp::CHEB:
  case FFOp::EXP:  case FFOp::LOG:   case FFOp::FABS: case FFOp::SQR:
  case FFOp::SQRT: case FFOp::COS:   case FFOp::SIN:  case FFOp::TAN:
  case FFOp::ACOS: case FFOp::ASIN:  case FFOp::ATAN: case FFOp::ERF:
  case FFOp::FSTEP:
  default:
    return true;
  }
}

///////////////////////////////// FFGraph //////////////////////////////////////

inline std::ostream&
operator <<
( std::ostream&out, const FFGraph&dag)
{
  typename FFGraph::t_Vars Vars = dag._Vars;
  typename FFGraph::it_Vars itv = Vars.begin();

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

inline FFVar&
FFGraph::_insert_binary_operation
( const typename FFOp::TYPE top, const FFDep&dep, const FFVar&Var1, const FFVar&Var2 )
{
  if( Var1._dag != Var2._dag ) throw Exceptions( Exceptions::DAG );
  FFGraph *pdag = Var1._dag;
  FFVar *pVar1 = Var1._ops.first->pres, *pVar2 = Var2._ops.first->pres;

  FFOp *pOp = pdag->_insert_operation( top, pVar1, pVar2 );
  if( pOp->pres ) return *pOp->pres;
  pVar1->_ops.second.push_back( pOp );
  pVar2->_ops.second.push_back( pOp );
  pOp->pres = pdag->_add_auxiliary( dep, pOp );
  return *pOp->pres;
}

template <typename U> inline FFVar&
FFGraph::_insert_binary_operation
( const typename FFOp::TYPE top, const FFDep&dep, const U&Cst1, const FFVar&Var2 )
{
  FFGraph *pdag = Var2._dag;
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
FFGraph::_insert_binary_operation
( const typename FFOp::TYPE top, const FFDep&dep, const FFVar&Var1, const U&Cst2 )
{
  FFGraph *pdag = Var1._dag;
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
FFGraph::_insert_unary_operation
( const typename FFOp::TYPE top, const FFDep&dep, const FFVar&Var )
{
  FFGraph* pdag = Var._dag;
  FFVar *pVar = Var._ops.first->pres;

  FFOp* pOp = pdag->_insert_operation( top, Var._ops.first->pres );
  if( pOp->pres ) return *pOp->pres;
  pVar->_ops.second.push_back( pOp );
  pOp->pres = pdag->_add_auxiliary( dep, pOp );
  return *pOp->pres;
}

inline FFOp*
FFGraph::_insert_operation
( const typename FFOp::TYPE top, FFVar*lop, FFVar*rop )
{
  FFOp* op = new FFOp( top, lop, rop );
  typename FFGraph::it_Ops itop = _Ops.find( op );
  if( itop!=_Ops.end() ){ delete op; return *itop; }
  _Ops.insert( op );
  return op;
}

inline bool
FFGraph::_remove_operation
( FFOp* op )
{
  typename FFGraph::it_Ops itop = _Ops.find( op );
  if( itop==_Ops.end() ) return false;
  delete op;
  _Ops.erase( itop );
  return true;
}
 
inline FFVar*
FFGraph::_add_auxiliary
( const FFDep&dep, FFOp*pOp )
{
  pOp->pres = new FFVar( this, dep, pOp );
  _append_aux( pOp->pres );
  return pOp->pres;
}

inline FFVar*
FFGraph::_set_constant
( const FFVar*pVar, const FFNum&num )
{
  it_Vars itVar = _Vars.find( const_cast<FFVar*>(pVar) );
  if( itVar!=_Vars.end() ){
    (*itVar)->_num = num;
    (*itVar)->_cst=true;
  }
  return *itVar;
}

inline FFVar*
FFGraph::_unset_constant
( const FFVar*pVar )
{
  it_Vars itVar = _Vars.find( const_cast<FFVar*>(pVar) );
  if( itVar!=_Vars.end() ) (*itVar)->_cst=false;
  return *itVar;
}

inline FFVar*
FFGraph::_add_constant
( const double x )
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
FFGraph::_add_constant
( const int n )
{
  // Check if real constant x already defined in _Vars
  FFVar* pAux = new FFVar( n );
  it_Vars iAux = _Vars.find( pAux );
  if( iAux!=_Vars.end() ){ delete pAux; return *iAux; }

  // Otherwise, append constant n
  _append_aux( pAux, FFOp::CNST );
  return pAux;
}

inline void
FFGraph::_append_aux
( FFVar*pAux, typename FFOp::TYPE tOp )
{
  FFOp*pOp = new FFOp( tOp, 0, 0, pAux );
  _Ops.insert( pOp );
  pAux->dag() = this;
  pAux->ops().first = pOp;
  pAux->id().second = _naux++;
  _append_aux( pAux );
}

inline void
FFGraph::_append_aux
( FFVar*pAux )
{
  _Vars.insert( pAux );
}

inline void
FFGraph::_append_var
( FFVar*pVar )
{
  _Vars.insert( pVar );
}   

inline FFVar*
FFGraph::_find_var
( const typename FFVar::pt_idVar&id )
{
  FFVar* pVar = new FFVar( this, id );
  it_Vars iVar = _Vars.find( pVar );
  delete pVar;
  return( iVar==_Vars.end()? 0: *iVar );
}

inline std::list<const FFOp*>
FFGraph::subgraph
( const unsigned int nDep, const FFVar*pDep ) const
{
  std::list<const FFOp*> Ops;
  _reset_operations();
  for( unsigned int i=0; i<nDep; i++ ){
    if( !pDep[i].ops().first ) continue;
    pDep[i].ops().first->propagate_subgraph( Ops );
  }
  return Ops;
}

inline std::list<const FFOp*>
FFGraph::subgraph
( const std::vector<const FFVar*>&vDep ) const
{
  //return subgraph( vDep.size(), vDep.data() );

  std::list<const FFOp*> Ops;
  _reset_operations();
  typename std::vector<const FFVar*>::const_iterator it = vDep.begin();
  for( ; it!=vDep.end(); ++it ){
    if( !(*it)->ops().first ) continue;
    (*it)->ops().first->propagate_subgraph( Ops );
  }
  return Ops;
}

inline void
FFGraph::output
( const std::list<const FFOp*>&Ops, std::ostream&os )
{
  os << ( !Ops.empty()? "\nFACTORS IN SUBGRAPH:\n": "\nNO FACTORS IN SUBGRAPH\n" );
  typename std::list<const FFOp*>::const_iterator ito = Ops.begin();
  for( ; ito!=Ops.end(); ++ito )
    os << "  " << *(*ito)->pres << "\t" << "<=  " << **ito << std::endl;
}

inline void
FFGraph::dot_script
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
FFGraph::dot_script
( const std::vector<const FFVar*>&vDep, std::ostream&os ) const
{
  //dot_script( vDep.size(), *vDep.data(), os )

  _reset_operations();
  os << "\ndigraph G {\nnode [shape=record];\n";
  typename std::vector<const FFVar*>::const_iterator it = vDep.begin();
  for( ; it!=vDep.end(); ++it ){
    if( !(*it)->ops().first ) continue;
    (*it)->ops().first->generate_dot_script( os );
  }
  os << "}\n";
}

inline std::tuple< unsigned, const unsigned*, const unsigned*, const FFVar* >
FFGraph::SFAD
( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
  const FFVar* const pIndep, const FFVar* const pDir )
{
  if( !nDep || !nIndep ) return std::make_tuple(0,(unsigned*)0,(unsigned*)0, (FFVar*)0);  // Nothing to do!
  assert( pDep && pIndep );

  std::vector<const FFVar*> vDep, vIndep, vDir;
  for( unsigned i=0; i<nDep; i++ )   vDep.push_back( pDep+i );
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );
  for( unsigned i=0; pDir && i<nIndep; i++ ) vDir.push_back( pDir+i );

  auto vDep_F = SFAD( vDep, vIndep, vDir );
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

inline std::tuple< unsigned, const unsigned*, const unsigned*, const FFVar* >
FFGraph::SFAD
( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
  const FFVar* const pIndep, const bool LUopt )
{
  if( !nDep || !nIndep ) return std::make_tuple(0,(unsigned*)0,(unsigned*)0, (FFVar*)0);  // Nothing to do!
  assert( pDep && pIndep );

  std::vector<const FFVar*> vDep, vIndep;
  for( unsigned i=0; i<nDep; i++ )   vDep.push_back( pDep+i );
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );

  auto vDep_F = SFAD( vDep, vIndep );
  unsigned nDep_F = 0;
  for( unsigned iv=0; iv<std::get<0>(vDep_F).size(); iv++ ){
    if( (  LUopt && std::get<0>(vDep_F)[iv] < std::get<1>(vDep_F)[iv] )
     || ( !LUopt && std::get<0>(vDep_F)[iv] > std::get<1>(vDep_F)[iv] ) )
      continue;
    nDep_F++;
  }
  unsigned* iDep_F = new unsigned[ nDep_F ];
  unsigned* jDep_F = new unsigned[ nDep_F ];
  FFVar* pDep_F = new FFVar[ nDep_F ];
  for( unsigned iv=0, ie=0; iv<std::get<0>(vDep_F).size(); iv++ ){
    if( (  LUopt && std::get<0>(vDep_F)[iv] < std::get<1>(vDep_F)[iv] )
     || ( !LUopt && std::get<0>(vDep_F)[iv] > std::get<1>(vDep_F)[iv] ) )
      continue;
    iDep_F[ie] = std::get<0>(vDep_F)[iv];
    jDep_F[ie] = std::get<1>(vDep_F)[iv];
    pDep_F[ie++] = *std::get<2>(vDep_F)[iv];
  }
  return std::make_tuple( nDep_F, iDep_F, jDep_F, pDep_F );
}

inline const FFVar*
FFGraph::FAD
( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
  const FFVar* const pIndep, const bool transp )
{
  if( !nDep || !nIndep ) return 0;  // Nothing to do!
  assert( pDep && pIndep );

  std::vector<const FFVar*> vDep, vIndep;
  for( unsigned i=0; i<nDep; i++ )   vDep.push_back( pDep+i );
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );

  std::vector<const FFVar*> vDep_F = FAD( vDep, vIndep );
  FFVar* pDep_F = new FFVar[ vDep_F.size() ];
  typename std::vector<const FFVar*>::const_iterator it = vDep_F.begin();
  if( !transp )
    for( unsigned k=0; it!=vDep_F.end(); ++it, k++ ) pDep_F[k] = **it;
  else
    for( unsigned i=0; it!=vDep_F.end(); i++ )
      for( unsigned j=0; j<nIndep; ++it, j++ ) pDep_F[i+j*nDep] = **it;

  return pDep_F;
}

inline const FFVar*
FFGraph::FAD
( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
  const FFVar* const pIndep, const FFVar* const pDir )
{
  if( !nDep || !nIndep ) return 0;  // Nothing to do!
  assert( pDep && pIndep );

  std::vector<const FFVar*> vDep, vIndep, vDir;
  for( unsigned i=0; i<nDep; i++ )   vDep.push_back( pDep+i );
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );
  for( unsigned i=0; pDir && i<nIndep; i++ ) vDir.push_back( pDir+i );

  std::vector<const FFVar*> vDep_F = FAD( vDep, vIndep, vDir );
  FFVar* pDep_F = new FFVar[ vDep_F.size() ];
  typename std::vector<const FFVar*>::const_iterator it = vDep_F.begin();
  for( unsigned k=0; it!=vDep_F.end(); ++it, k++ ) pDep_F[k] = **it;

  return pDep_F;
}

inline std::vector<const FFVar*>
FFGraph::FAD
( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
  const std::vector<const FFVar*>&vDir )
{
  auto sDep_F = SFAD( vDep, vIndep, vDir );
  const FFVar*pZero = _add_constant( 0. );
  std::vector<const FFVar*> vDep_F( vIndep.size()*vDep.size(), pZero );
  for( unsigned ie(0); ie<std::get<2>(sDep_F).size(); ie++ ){
    unsigned pDep_F = std::get<0>(sDep_F)[ie]*vIndep.size()+std::get<1>(sDep_F)[ie];
    vDep_F[pDep_F] = std::get<2>(sDep_F)[ie];
  }
  return vDep_F;
}

inline std::vector<const FFVar*>
FFGraph::FAD
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

inline std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> >
FFGraph::SFAD
( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep,
  const std::vector<const FFVar*>&vDir )
{
  // Nothing to do!
  if( !vIndep.size() || !vDep.size() ) return std::make_tuple( std::vector<unsigned>(),
    std::vector<unsigned>(), std::vector<const FFVar*>() );
  assert( !vDir.size() || vIndep.size() == vDir.size() );
  //fadbad::F<mc::FFVar> FFVar_dum();

  // Initialize of all independent variables participating in the dependent ones
  it_Vars itv = _Vars.begin();
  for( ; itv!=_Vars.end() && (*itv)->_id.first<=FFVar::VAR; ++itv ){
    fadbad::F<mc::FFVar>* pX_F = new fadbad::F<mc::FFVar>( **itv );
    auto iti = vIndep.begin();
    auto itd = vDir.begin();
    for( unsigned int i=0; iti!=vIndep.end(); ++iti, ++itd, i++ )
      if( (*itv)->id().second == (*iti)->id().second ){
        if( vDir.size() ) pX_F->diff( 0, 1 ) = **itd;
        else pX_F->diff( i, vIndep.size() ); 
      }
    // Attach fadbad::F<mc::FFVar>* variable to corresponding variable in _Vars
    (*itv)->val() = pX_F;
  }
  // THIS IS DOING A BIT TOO MUCH WORK AS ONLY THE INDEPENDENT VARIABLES
  // PARTICIPATING IN THE DEPENDENTS SHOULD BE TAKEN INTO ACCOUNT REALLY
  // (E.G., THIS COULD BE DONE USING THE .dep() FIELD IN THE DEPENDENTS)

  // Evaluate dependents given by vDep in fadbad::F type
  std::list<const FFOp*> opDep = subgraph( vDep );
  auto ito = opDep.begin();
  for( ; ito!=opDep.end(); ++ito ){
    _curOp = *ito;
    _curOp->evaluate( fadbad::F<mc::FFVar>() );
  }
  // Retreive dependents variables in fadbad::F as given by vDep
  std::tuple< std::vector<unsigned>, std::vector<unsigned>,
              std::vector<const FFVar*> > vDep_F; // <- vector holding the results in sparse format
  auto itd = vDep.begin();
  for( unsigned i=0; itd!=vDep.end(); ++itd, i++ ){
    // Obtain pointer to dependent variable in FFGraph
    FFVar* pF = _find_var( (*itd)->id() );
    auto iti = vIndep.begin();
    // Push corresponding evaluation in fadbad::F into result vector
    for( unsigned j=0; iti!=vIndep.end(); ++iti, j++ ){
      if( !pF ) continue;
      // THE FOLLOWING MATCHING IS NECESSARY BECAUSE THE VARIABLES CREATED
      // BY FFOp::evaluate ARE NOT THE SAME AS THOSE STORED IN FFGraph
      fadbad::F<mc::FFVar>* pF_F = static_cast<fadbad::F<mc::FFVar>*>( pF->val() );
      const FFVar* pdFdX = _find_var( pF_F->deriv(j).id() );
      if( !pdFdX ){
        const FFNum& num = pF_F->deriv(j).num();
        switch( num.t ){
          case FFNum::INT:  if( num.n )       pdFdX = _add_constant( num.n ); break;
          case FFNum::REAL: if( num.x != 0. ) pdFdX = _add_constant( num.x ); break;
        }
        if( !pdFdX ) continue;
      }
      std::get<0>(vDep_F).push_back( i ); // add row index
      std::get<1>(vDep_F).push_back( j ); // add column index
      std::get<2>(vDep_F).push_back( pdFdX ); // add Jacobian element
      if( vDir.size() ) break; // interrupt if directional derivatives requested
    }
  }

  // Reset mc::FFVAR_val field to NULL
  itv = _Vars.begin();
  for( ; itv!=_Vars.end() && (*itv)->_id.first<=FFVar::VAR; ++itv )
    (*itv)->reset_val( fadbad::F<mc::FFVar>() );

  _reset_operations();
  itd = vDep.begin();
  for( ; itd!=vDep.end(); ++itd ){
    if( !(*itd)->ops().first ) continue;
    (*itd)->ops().first->reset_val_subgraph( fadbad::F<mc::FFVar>() );
  }

  return vDep_F;
}

inline std::tuple< unsigned, const unsigned*, const unsigned*, const FFVar* >
FFGraph::SBAD
( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
  const FFVar* const pIndep )
{
  if( !nDep || !nIndep ) return std::make_tuple(0,(unsigned*)0,(unsigned*)0, (FFVar*)0);  // Nothing to do!
  assert( pDep && pIndep );

  std::vector<const FFVar*> vDep, vIndep;
  for( unsigned i=0; i<nDep; i++ )   vDep.push_back( pDep+i );
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );

  auto vDep_B = SBAD( vDep, vIndep );
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

inline std::tuple< unsigned, const unsigned*, const unsigned*, const FFVar* >
FFGraph::SBAD
( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
  const FFVar* const pIndep, const bool LUopt )
{
  if( !nDep || !nIndep ) return std::make_tuple(0,(unsigned*)0,(unsigned*)0, (FFVar*)0);  // Nothing to do!
  assert( pDep && pIndep );

  std::vector<const FFVar*> vDep, vIndep;
  for( unsigned i=0; i<nDep; i++ )   vDep.push_back( pDep+i );
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );

  auto vDep_B = SBAD( vDep, vIndep );
  unsigned nDep_B = 0;
  for( unsigned iv=0; iv<std::get<0>(vDep_B).size(); iv++ ){
    if( (  LUopt && std::get<0>(vDep_B)[iv] < std::get<1>(vDep_B)[iv] )
     || ( !LUopt && std::get<0>(vDep_B)[iv] > std::get<1>(vDep_B)[iv] ) )
      continue;
    nDep_B++;
  }
  unsigned* iDep_B = new unsigned[ nDep_B ];
  unsigned* jDep_B = new unsigned[ nDep_B ];
  FFVar* pDep_B = new FFVar[ nDep_B ];
  for( unsigned iv=0, ie=0; iv<std::get<0>(vDep_B).size(); iv++ ){
    if( (  LUopt && std::get<0>(vDep_B)[iv] < std::get<1>(vDep_B)[iv] )
     || ( !LUopt && std::get<0>(vDep_B)[iv] > std::get<1>(vDep_B)[iv] ) )
      continue;
    iDep_B[ie] = std::get<0>(vDep_B)[iv];
    jDep_B[ie] = std::get<1>(vDep_B)[iv];
    pDep_B[ie++] = *std::get<2>(vDep_B)[iv];
  }
  return std::make_tuple( nDep_B, iDep_B, jDep_B, pDep_B );
}

inline const FFVar*
FFGraph::BAD
( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
  const FFVar* const pIndep, const bool transp )
{
  if( !nDep || !nIndep ) return 0;  // Nothing to do!
  assert( pDep && pIndep );

  std::vector<const FFVar*> vDep, vIndep;
  for( unsigned i=0; i<nDep; i++ )   vDep.push_back( pDep+i );
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );

  std::vector<const FFVar*> vDep_B = BAD( vDep, vIndep );
  FFVar* pDep_B = new FFVar[ vDep_B.size() ];
  typename std::vector<const FFVar*>::const_iterator it = vDep_B.begin();
  if( !transp )
    for( unsigned k=0; it!=vDep_B.end(); ++it, k++ ) pDep_B[k] = **it;
  else
    for( unsigned i=0; it!=vDep_B.end(); i++ )
      for( unsigned j=0; j<nIndep; ++it, j++ ) pDep_B[i+j*nDep] = **it;

  return pDep_B;
}

inline std::vector<const FFVar*>
FFGraph::BAD
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
/*
inline std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> >
FFGraph::SBAD
( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep )
{
  // Nothing to do!
  if( !vIndep.size() || !vDep.size() ) return std::make_tuple( std::vector<unsigned>(),
    std::vector<unsigned>(), std::vector<const FFVar*>() );

  // Initialize of all independent variables participating in the dependent ones
  it_Vars itv = _Vars.begin();
  for( ; itv!=_Vars.end() && (*itv)->_id.first<=FFVar::VAR; ++itv ){
    fadbad::B<mc::FFVar>* pX_B = new fadbad::B<mc::FFVar>( **itv );
    // Attach fadbad::B<mc::FFVar>* variable to corresponding variable in _Vars
    (*itv)->val() = pX_B;
    //std::cout << "*itv: " << *itv << "   pX_B: " << pX_B << std::endl;
  }
  // THIS IS DOING A BIT TOO MUCH WORK AS ONLY THE INDEPENDENT VARIABLES
  // PARTICIPATING IN THE DEPENDENTS SHOULD BE TAKEN INTO ACCOUNT REALLY
  // (E.G., THIS COULD BE DONE USING THE .dep() FIELD IN THE DEPENDENTS)

  // Propagate dependents given by vDep in fadbad::B type
  std::list<const FFOp*> opDep = subgraph( vDep );
  typename std::list<const FFOp*>::iterator ito = opDep.begin();
  for( ; ito!=opDep.end(); ++ito ) (*ito)->evaluate( fadbad::B<mc::FFVar>() );

  // Propagate derivatives of dependents using fadbad::B
  typename std::vector<const FFVar*>::const_iterator itd = vDep.begin();
  //itd = vDep.begin();
  for( unsigned j=0; itd!=vDep.end(); ++itd, j++ ){
    // Obtain pointer to dependent variable in FFGraph
    FFVar* pF = _find_var( (*itd)->id() );
    //std::cout << "pF: " << pF;
    if( !pF ) continue;
    fadbad::B<mc::FFVar>* pF_B = static_cast<fadbad::B<mc::FFVar>*>( pF->val() );
    pF_B->diff( j, vDep.size() );
    *pF_B = 0.;
    //fadbad::B<mc::FFVar> F_B = *static_cast<fadbad::B<mc::FFVar>*>( pF->val() );
    //F_B.diff( j, vDep.size() );
    //std::cout << "   pF_B: " << pF_B << std::endl;
  }

  // Erase pointers to intermediate so they go out of scope
  // OTHERWISE BADIFF COMPLAINS
  _reset_operations();
  itd = vDep.begin();
  //typename std::vector<const FFVar*>::const_iterator itd = vDep.begin();
  for( ; itd!=vDep.end(); ++itd ){
    if( !(*itd)->ops().first ) continue;
    (*itd)->ops().first->reset_val_subgraph( fadbad::B<mc::FFVar>(), vDep, vIndep );
  }

  // Retreive dependents variables in fadbad::B as given by vDep
  std::tuple< std::vector<unsigned>, std::vector<unsigned>,
              std::vector<const FFVar*> > vDep_B; // <- vector holding the results in sparse format
  itd = vDep.begin();
  for( unsigned i=0; itd!=vDep.end(); ++itd, i++ ){
    // Obtain pointer to dependent variable in FFGraph
    FFVar* pF = _find_var( (*itd)->id() );
    typename std::vector<const FFVar*>::const_iterator iti = vIndep.begin();
    // Push corresponding evaluation in fadbad::B into result vector
    for( unsigned j=0; iti!=vIndep.end(); ++iti, j++ ){
      if( !pF ) continue;
      // Obtain pointer to independent variable in FFGraph
      FFVar* pX = _find_var( (*iti)->id() );
      // THE FOLLOWING MATCHING IS NECESSARY BECAUSE THE VARIABLES CREATED
      // BY FFOp::evaluate ARE NOT THE SAME AS THOSE STORED IN FFGraph
      fadbad::B<mc::FFVar>* pX_B = static_cast<fadbad::B<mc::FFVar>*>( pX->val() );
      std::cout << "before\n";
      FFVar* dXj = &pX_B->d(i);
      std::cout << "after\n";
      const FFVar* pdFdX = 0;
      // Check if dependent matches the independent
      //if( pX == pF )
      //  pdFdX = _add_constant( 1. );
      //else{
        pdFdX = _find_var( dXj->id() );
        if( !pdFdX ){
          const FFNum& num = dXj->num();
          switch( num.t ){
            case FFNum::INT:  if( num.n )       pdFdX = _add_constant( num.n ); break;
            case FFNum::REAL: if( num.x != 0. ) pdFdX = _add_constant( num.x ); break;
          }
          if( !pdFdX ) continue;
        }
      //}
      std::get<0>(vDep_B).push_back( i ); // add row index
      std::get<1>(vDep_B).push_back( j ); // add column index
      std::get<2>(vDep_B).push_back( pdFdX ); // add Jacobian element
    }
  }

  // Reset mc::FFVAR::_val field to NULL in dependent and independent variables
  itv = _Vars.begin();
  for( ; itv!=_Vars.end() && (*itv)->_id.first<=FFVar::VAR; ++itv )
    (*itv)->reset_val( fadbad::B<mc::FFVar>() );

  itd = vDep.begin();
  for( ; itd!=vDep.end(); ++itd ){
    FFVar*pF = _find_var( (*itd)->id() );
    if( pF ) pF->reset_val( fadbad::B<mc::FFVar>() );
  }

  return vDep_B;
}
*/

inline std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<const FFVar*> >
FFGraph::SBAD
( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep )
{
  // Nothing to do!
  if( !vIndep.size() || !vDep.size() ) return std::make_tuple( std::vector<unsigned>(),
    std::vector<unsigned>(), std::vector<const FFVar*>() );

  // Initialize of all independent variables participating in the dependent ones
  std::vector< mc::FFVar > vVars, vDeps;
  std::vector< fadbad::B<mc::FFVar> > vVars_B, vDeps_B( vDep.size() );
  for( auto itv = _Vars.begin(); itv!=_Vars.end() && (*itv)->_id.first<=FFVar::VAR; ++itv ){
    vVars.push_back( **itv );
    vVars_B.push_back( **itv );
  }
  for( auto itd = vDep.begin(); itd!=vDep.end(); ++itd )
    vDeps.push_back( **itd );
  // THIS IS DOING A BIT TOO MUCH WORK AS ONLY THE INDEPENDENT VARIABLES
  // PARTICIPATING IN THE DEPENDENTS SHOULD BE TAKEN INTO ACCOUNT REALLY
  // (E.G., THIS COULD BE DONE USING THE .dep() FIELD IN THE DEPENDENTS)

  // Propagate dependents given by vDep in fadbad::B type
  eval( vDeps.size(), vDeps.data(), vDeps_B.data(), vVars.size(), vVars.data(), vVars_B.data() );
  auto itd_B = vDeps_B.begin();
  for( unsigned i=0; itd_B!=vDeps_B.end(); ++itd_B, i++ )
    (*itd_B).diff( i, vDeps_B.size() );

  // Retreive dependents variables in fadbad::B as given by vDep
  std::tuple< std::vector<unsigned>, std::vector<unsigned>,
              std::vector<const FFVar*> > vDep_B; // <- vector holding the results in sparse format

  auto itd = vDeps.begin();
  for( unsigned i=0; itd!=vDeps.end(); ++itd, i++ ){
    // Obtain pointer to dependent variable in FFGraph
    FFVar* pF = _find_var( (*itd).id() );
    if( !pF ) continue;

    auto itv = vVars.begin();
    for( unsigned k=0; itv!=vVars.end(); ++itv, k++ ){
      const FFVar* pdFdX = 0;
      unsigned j=0;
      for( auto iti=vIndep.begin(); iti!=vIndep.end(); ++iti, j++ ){
        if( (*iti)->id() == (*itv).id() ){
          FFVar dXj = vVars_B[k].d(i);
          pdFdX = _find_var( dXj.id() );
          if( !pdFdX ){
            const FFNum& num = dXj.num();
            switch( num.t ){
              case FFNum::INT:  if( num.n )       pdFdX = _add_constant( num.n ); break;
              case FFNum::REAL: if( num.x != 0. ) pdFdX = _add_constant( num.x ); break;
            }
          }
          break;
        }
      }
      if( !pdFdX ) continue;
      std::get<0>(vDep_B).push_back( i ); // add row index
      std::get<1>(vDep_B).push_back( j ); // add column index
      std::get<2>(vDep_B).push_back( pdFdX ); // add Jacobian element
    }
  }

  return vDep_B;
}

inline const FFVar*
FFGraph::TAD
( const unsigned int ordermax, const unsigned nDep, const FFVar* const pDep,
  const unsigned nVar, const FFVar* const pVar, const FFVar* const pIndep )
{
  if( !nDep || !nVar ) return 0;  // Nothing to do!
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

inline std::vector<const FFVar*>
FFGraph::TAD
( const unsigned int ordermax, const std::vector<const FFVar*>&vDep,
  const std::vector<const FFVar*>&vVar, const FFVar* const pIndep )
{
  // Check dependent and independent vector sizes
  if( !vVar.size() || !vDep.size() || vVar.size() != vDep.size() )
    return std::vector<const FFVar*>();
  //fadbad::T<mc::FFVar> FFVar_dum();
  std::vector<const FFVar*> vDep_T; // <- vector holding the results

  // Initialize of all independent variables participating in the dependent ones
  fadbad::T<mc::FFVar>** pX_T = new fadbad::T<mc::FFVar>*[ vVar.size() ];
  it_Vars itv = _Vars.begin();
  for( ; itv!=_Vars.end() && (*itv)->_id.first<=FFVar::VAR; ++itv ){
    fadbad::T<mc::FFVar>* pXi_T = new fadbad::T<mc::FFVar>( **itv );
    typename std::vector<const FFVar*>::const_iterator iti = vVar.begin();
    // Time
    if( pIndep && (*itv)->id().second == pIndep->id().second )
      (*pXi_T)[1] = 1.;
    // Time-dependents
    for( unsigned int i=0; iti!=vVar.end(); ++iti, i++ ){
      if( (*itv)->id().second == (*iti)->id().second ){
        pX_T[i] = pXi_T;
        vDep_T.push_back( *itv ); // <- Append 0th-order Taylor coefficient of ith-dependent to result vector
#ifdef MC__FFUNC_DEBUG_TAD
        std::cout << "FFGraph::TAD *** f(" << i << ")[0] = "
                  << **itv << "  (" << *itv << ")\n";
#endif
      }
    }
    // Attach fadbad::F<mc::FFVar>* variable to corresponding variable in _Vars
    (*itv)->val() = pXi_T;
  }
  // THIS IS DOING A BIT TOO MUCH WORK AS ONLY THE INDEPENDENT VARIABLES
  // PARTICIPATING IN THE DEPENDENTS SHOULD BE TAKEN INTO ACCOUNT REALLY
  // (E.G., THIS COULD BE DONE USING THE .dep() FIELD IN THE DEPENDENTS)

  // Evaluate dependents given by vDep in fadbad::T type
  std::list<const FFOp*> opDep = subgraph( vDep );
  typename std::list<const FFOp*>::iterator ito = opDep.begin();
  for( ; ito!=opDep.end(); ++ito ){
    _curOp = *ito;
    _curOp->evaluate( fadbad::T<mc::FFVar>() );
  }

  // Set pointers to the dependents
  fadbad::T<mc::FFVar>** pF_T = new fadbad::T<mc::FFVar>*[ vDep.size() ];
  typename std::vector<const FFVar*>::const_iterator itd = vDep.begin();
  for( unsigned j=0; itd!=vDep.end(); ++itd, j++ ){
    FFVar*pF = _find_var( (*itd)->id() );
    pF_T[j] = ( pF? static_cast<fadbad::T<mc::FFVar>*>( pF->val() ): 0 );
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
      mc::FFVar Xjq = (*pF_T[j])[q]/double(q+1);
      //mc::FFVar& Xjq = (*pF_T[j])[q];
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
        vDep_T.push_back( pXjq ); // <- Append (q+1)th-order Taylor coefficient of jth-dependent to result vector
      }
#ifdef MC__FFUNC_DEBUG_TAD
      std::cout << "FFGraph::TAD *** f(" << j << ")[" << q+1 << "] = "
                << *pXjq << "  (" << pXjq << ")\n";
#endif
    }
  }

  // Reset mc::FFVAR_val field to NULL
  itv = _Vars.begin();
  for( ; itv!=_Vars.end() && (*itv)->_id.first<=FFVar::VAR; ++itv )
    (*itv)->reset_val( fadbad::T<mc::FFVar>() );

  _reset_operations();
  itd = vDep.begin();
  for( ; itd!=vDep.end(); ++itd ){
    if( !(*itd)->ops().first ) continue;
    (*itd)->ops().first->reset_val_subgraph( fadbad::T<mc::FFVar>() );
  }

  delete[] pX_T;
  delete[] pF_T;

  return vDep_T;
}

inline const FFVar*
FFGraph::compose
( const unsigned nDepOut, const FFVar*pDepOut, const unsigned nDepIn,
  const FFVar*pVarOut,  const FFVar*pDepIn )
{
  if( !nDepOut || !nDepIn ) return pDepOut;  // Nothing to do!
  assert( pDepOut && pVarOut && pDepIn );

  std::vector<const FFVar*> vDepOut;
  std::vector< std::pair<const FFVar*,const FFVar*> > vDepIn;
  for( unsigned i=0; i<nDepOut; i++ ) vDepOut.push_back( pDepOut+i );
  for( unsigned i=0; i<nDepIn; i++ ) vDepIn.push_back( std::make_pair(pVarOut+i,pDepIn+i) );

  std::vector<const FFVar*> vDepComp = compose( vDepOut, vDepIn );

  FFVar* pDepComp = new FFVar[ vDepComp.size() ];
  typename std::vector<const FFVar*>::const_iterator it = vDepComp.begin();
  for( unsigned k=0; it!=vDepComp.end(); ++it, k++ ) pDepComp[k] = **it;

  return pDepComp;
}

inline std::vector<const FFVar*>
FFGraph::compose
( const std::vector<const FFVar*>&vDepOut,
  const std::vector< std::pair<const FFVar*,const FFVar*> >&vDepIn )
{
  // Nothing to do!
  if( !vDepIn.size() || !vDepOut.size() ) return vDepOut;

  // Initialize of all independent variables participating in the dependent ones
  std::vector< std::pair<const FFVar*,FFVar> > vIndep;
  it_Vars itv = _Vars.begin();
  for( ; itv!=_Vars.end() && (*itv)->_id.first<=FFVar::VAR; ++itv ){
    typename std::vector< std::pair<const FFVar*,const FFVar*> >::const_iterator iti = vDepIn.begin();
    bool match = false;
    // Pair with corresponding dependent if matching
    for( ; !match && iti!=vDepIn.end(); ++iti )
      if( (*itv)->id().second == (*iti).first->id().second ){
        FFVar DepIn = *(*iti).second; // to avoid constness issue
        vIndep.push_back( std::make_pair( (*iti).first, DepIn ) ); //_find_var( (*iti).second->id() );
#ifdef MC__FFUNC_DEBUG_COMPOSE
        std::cout << *(*iti).first << " <-> " << DepIn << std::endl;
#endif
        match = true;
      }
    // Pair with itself otherwise
    if( !match ){
      vIndep.push_back( std::make_pair( *itv, **itv ) );
#ifdef MC__FFUNC_DEBUG_COMPOSE
      std::cout << **itv << " <-> " << **itv << std::endl;
#endif
    }
  }

  std::vector<FFVar> vDepComp = eval( vDepOut, vIndep );
  std::vector<const FFVar*> vDepNew;
  typename std::vector<FFVar>::iterator itd = vDepComp.begin();
  for( ; itd!=vDepComp.end(); ++itd ){
    FFVar* pF = _find_var( (*itd).id() );
    if( !pF ){
      const FFNum& num = (*itd).num();
      switch( num.t ){
        case FFNum::INT:  vDepNew.push_back( _add_constant( num.n ) ); continue;
        case FFNum::REAL: vDepNew.push_back( _add_constant( num.x ) ); continue;
      }
    }
    else vDepNew.push_back( pF );
  }
  return vDepNew;
}

inline void
FFGraph::decompose
( const FFVar&pDep, unsigned&nLin, double*cLin, FFVar*pLin,
  unsigned&nNLin, double*cNLin, FFVar*pNLin )
{
}

template <typename U> inline std::vector<U>
FFGraph::eval
( const std::vector<const FFVar*>&vDep,
  const std::vector< std::pair<const FFVar*,U> >&vVar )
{
  // Nothing to do!
  if( !vDep.size() ) return std::vector<U>();

  // Generate subgraph -- This can be the most time consuming step!!!
  std::list<const FFOp*> opDep = subgraph( vDep );

  return eval( opDep, vDep, vVar );
}

template <typename U> inline std::vector<U>
FFGraph::eval
( std::list<const FFOp*>&opDep, const std::vector<const FFVar*>&vDep,
  const std::vector< std::pair<const FFVar*,U> >&vVar )
{
  // Nothing to do!
  if( !vDep.size() ) return std::vector<U>();
#ifdef MC__FFUNC_CPU_EVAL
  double cputime;
#endif

  // Initialize all independent variables participating in the dependent ones
#ifdef MC__FFUNC_CPU_EVAL
  cputime = -cpuclock();
#endif
  auto iti = vVar.begin();
  for( ; iti!=vVar.end(); ++iti ){
    FFVar* pF = _find_var( (*iti).first->id() );
    if( pF ){
      pF->val() = new U( (*iti).second ); //const_cast<U*>( &(*iti).second );
#ifdef MC__FFUNC_DEBUG_EVAL
      std::cout << (*iti).first << "  " << (*iti).second << std::endl;
#endif
    }
  }
#ifdef MC__FFUNC_CPU_EVAL
  cputime += cpuclock();
  std::cout << "\nIndep. init. time: " << std::fixed << cputime << std::endl;
#endif
  // THIS IS DOING A BIT TOO MUCH WORK AS ONLY THE INDEPENDENT VARIABLES
  // PARTICIPATING IN THE DEPENDENTS SHOULD BE TAKEN INTO ACCOUNT STRICTLY
  // (E.G., THIS COULD BE DONE USING THE .dep() FIELD IN THE DEPENDENTS)

  std::vector<U> vDep_U; // <- vector holding the results
  bool ffexcp = false, galexcp = false;
  FFGraph::Exceptions ffexcpobj;
  try{
#ifdef MC__FFUNC_CPU_EVAL
    cputime = -cpuclock();
#endif
    // Evaluate dependents given by vDep in U type
    auto ito = opDep.begin();
    for( ; ito!=opDep.end(); ++ito ){
      _curOp = *ito;
      _curOp->evaluate( U() );
    }
#ifdef MC__FFUNC_CPU_EVAL
    cputime += cpuclock();
    std::cout << "Evaluation time: " << std::fixed << cputime << std::endl;
#endif

    // Retreive dependents variables in U type as given by vDep
#ifdef MC__FFUNC_CPU_EVAL
    cputime = -cpuclock();
#endif
    auto itd = vDep.begin();
    for( ; itd!=vDep.end(); ++itd ){
      // Obtain pointer to dependent variable in FFGraph
      FFVar* pF = _find_var( (*itd)->id() );
      // Push corresponding evaluation in U type into result vector
      if( pF ) vDep_U.push_back( U( *static_cast<U*>( pF->val() ) ) );
      else switch( (*itd)->num().t ){
        case FFNum::INT:
          vDep_U.push_back( (*itd)->num().n );
          break;
        case FFNum::REAL:
          vDep_U.push_back( (*itd)->num().x );
          break;
      }
    }
#ifdef MC__FFUNC_CPU_EVAL
    cputime += cpuclock();
    std::cout << "Dep. collect. time: " << std::fixed << cputime << std::endl;
#endif
  }
  catch( FFGraph::Exceptions &eObj ){
    ffexcp = true; ffexcpobj = eObj;
  }
  catch(...){
    galexcp = true;
  }

  // Reset mc::FFVAR_val field to NULL
#ifdef MC__FFUNC_CPU_EVAL
  cputime = -cpuclock();
#endif
  iti = vVar.begin();
  for( ; iti!=vVar.end(); ++iti ){
    FFVar* pF = _find_var( (*iti).first->id() );
    if( pF ) pF->reset_val( U() );
  }
  //it_Vars itv = _Vars.begin();
  //for( ; itv!=_Vars.end() && (*itv)->_id.first<=FFVar::VAR; ++itv )
  //  (*itv)->reset_val( U() );
  
  auto ito = opDep.begin();
  for( ; ito!=opDep.end(); ++ito ) (*ito)->flag( false );
  auto itd = vDep.begin();
  for( ; itd!=vDep.end(); ++itd ){
    if( !(*itd)->ops().first ) continue;
    (*itd)->ops().first->reset_val_subgraph( U() );
  }
#ifdef MC__FFUNC_CPU_EVAL
  cputime += cpuclock();
  std::cout << "Clean-up time: " << std::fixed << cputime << std::endl;
#endif

  if( ffexcp )  throw ffexcpobj;
  if( galexcp ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::EVAL );
  return vDep_U;
}

template <typename U> inline std::vector<U>
FFGraph::eval
( std::list<const FFOp*>&opDep, U*opRes, const std::vector<const FFVar*>&vDep,
  const std::vector< std::pair<const FFVar*,U> >&vVar )
{
  // Nothing to do!
  if( !vDep.size() ) return std::vector<U>();
  assert( opDep.size() && opRes );

  // Propagate values in U arithmetic through subgraph
#ifdef MC__FFUNC_CPU_EVAL
  double cputime = -cpuclock();
#endif
  auto ito = opDep.begin();
  for( U*pUres=opRes; ito!=opDep.end(); ++ito, pUres++ ){
    // Initialize variable using values in vVar
    if( (*ito)->type == FFOp::VAR ){
      FFVar* pF = 0;
      auto iti = vVar.begin();
      for( ; iti!=vVar.end(); ++iti ){
        if( (*ito)->pres->id() == (*iti).first->id() ){
          pF = (*ito)->pres; //_find_var( (*iti).first->id() );
          *pUres = (*iti).second;
          break;
        }
      }
      if( !pF ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::MISSVAR );
    }
    _curOp = *ito;
    _curOp->evaluate( pUres );
  }
#ifdef MC__FFUNC_CPU_EVAL
  cputime += cpuclock();
  std::cout << "\nEvaluation time: " << std::fixed << cputime << std::endl;
#endif

  // Retreive dependents variables in U type as given by vDep
#ifdef MC__FFUNC_CPU_EVAL
  cputime = -cpuclock();
#endif
  std::vector<U> vDep_U; // <- vector holding the results
  auto itd = vDep.begin();
  for( ; itd!=vDep.end(); ++itd ){
    // Obtain pointer to dependent variable in FFGraph
    FFVar* pF = _find_var( (*itd)->id() );
    // Push corresponding evaluation in U type into result vector
    if( pF ) vDep_U.push_back( U( *static_cast<U*>( pF->val() ) ) );
    else switch( (*itd)->num().t ){
      case FFNum::INT:
        vDep_U.push_back( (*itd)->num().n );
        break;
      case FFNum::REAL:
        vDep_U.push_back( (*itd)->num().x );
        break;
    }
  }
#ifdef MC__FFUNC_CPU_EVAL
  cputime += cpuclock();
  std::cout << "Dep. collect. time: " << std::fixed << cputime << std::endl;
#endif

  return vDep_U;
}

template <typename U> inline void
FFGraph::eval
( const unsigned nDep, const FFVar*pDep, U*vDep, const unsigned nVar,
  const FFVar*pVar, const U*vVar, const bool add )
{
  // Nothing to do!
  if( !nDep ) return;

  // Generate subgraph -- This can be the most time consuming step!!!
  std::list<const FFOp*> opDep = subgraph( nDep, pDep );

  return eval( opDep, nDep, pDep, vDep, nVar, pVar, vVar, add );
}

template <typename U> inline void
FFGraph::eval
( std::list<const FFOp*>&opDep, const unsigned nDep, const FFVar*pDep,
  U*vDep, const unsigned nVar, const FFVar*pVar, const U*vVar, const bool add )
{
  return eval( opDep, nDep, pDep, vDep, std::list<unsigned>(1,nVar),
               std::list<const FFVar*>(1,pVar), std::list<const U*>(1,vVar), add );
}

template <typename U> inline void
FFGraph::eval
( const unsigned nDep, const FFVar*pDep, U*vDep, const std::list<unsigned>&nVar,
  const std::list<const FFVar*>&pVar, const std::list<const U*>&vVar, const bool add )
{
  // Nothing to do!
  if( !nDep ) return;

  // Generate subgraph -- This can be the most time consuming step!!!
  std::list<const FFOp*> opDep = subgraph( nDep, pDep );

  return eval( opDep, nDep, pDep, vDep, nVar, pVar, vVar, add );
}

template <typename U> inline void
FFGraph::eval
( std::list<const FFOp*>&opDep, const unsigned nDep, const FFVar*pDep,
  U*vDep, const std::list<unsigned>&nVar, const std::list<const FFVar*>&pVar,
  const std::list<const U*>&vVar, const bool add )
{
  // Nothing to do!
  if( !nDep ) return;
  assert( pDep && vDep );
  const unsigned nIndep = nVar.size();
  assert( pVar.size() == nIndep && vVar.size() == nIndep );
  //assert( !nVar || ( pVar && vVar ) );
#ifdef MC__FFUNC_CPU_EVAL
  double cputime;
#endif

  // Initialize all independent variables participating in the dependent ones
#ifdef MC__FFUNC_CPU_EVAL
  cputime = -cpuclock();
#endif
  auto itnVar = nVar.begin(); auto itpVar = pVar.begin(); auto itvVar = vVar.begin();
  for( ; itnVar != nVar.end(); ++itnVar, ++itpVar, ++itvVar ){
    for( unsigned i=0; i<(*itnVar); i++ ){
      FFVar* pF = _find_var( (*itpVar)[i].id() );
      if( pF ){
        pF->val() = new U( (*itvVar)[i] );
        //std::cerr << "creating pF->val():" << pF->val() << std::endl;
#ifdef MC__FFUNC_DEBUG_EVAL
        std::cout << (*itpVar)[i] << "  " << (*itvVar)[i] << std::endl;
#endif
      }
    }
  }
#ifdef MC__FFUNC_CPU_EVAL
  cputime += cpuclock();
  std::cout << "\nIndep. init. time: " << std::fixed << cputime << std::endl;
#endif
  // THIS IS DOING A BIT TOO MUCH WORK AS ONLY THE INDEPENDENT VARIABLES
  // PARTICIPATING IN THE DEPENDENTS SHOULD BE TAKEN INTO ACCOUNT STRICTLY
  // (E.G., THIS COULD BE DONE USING THE .dep() FIELD IN THE DEPENDENTS)

  bool ffexcp = false, galexcp = false;
  FFGraph::Exceptions ffexcpobj;
  try{
#ifdef MC__FFUNC_CPU_EVAL
    cputime = -cpuclock();
#endif
    // Evaluate dependents given by vDep in U type
    auto ito = opDep.begin();
    for( ; ito!=opDep.end(); ++ito ){
      _curOp = *ito;
      _curOp->evaluate( U() );
    }
#ifdef MC__FFUNC_CPU_EVAL
    cputime += cpuclock();
    std::cout << "Evaluation time: " << std::fixed << cputime << std::endl;
#endif

    // Retreive dependents variables in U type as given by vDep
#ifdef MC__FFUNC_CPU_EVAL
    cputime = -cpuclock();
#endif
    for( unsigned i=0; i<nDep; i++ ){
      // Obtain pointer to dependent variable in FFGraph
      FFVar* pF = _find_var( pDep[i].id() );
      // Write/add corresponding evaluation in U type into/to result vector
      if( !add && pF ) vDep[i] = *static_cast<U*>( pF->val() );
      else if( pF )   vDep[i] += *static_cast<U*>( pF->val() );
      else if( !add ) switch( pDep[i].num().t ){
        case FFNum::INT:
          vDep[i] = pDep[i].num().n;
          break;
        case FFNum::REAL:
          vDep[i] = pDep[i].num().x;
          break;
      }
      else switch( pDep[i].num().t ){
        case FFNum::INT:
          vDep[i] += pDep[i].num().n;
          break;
        case FFNum::REAL:
          vDep[i] += pDep[i].num().x;
          break;
      }
    }
#ifdef MC__FFUNC_CPU_EVAL
    cputime += cpuclock();
    std::cout << "Dep. collect. time: " << std::fixed << cputime << std::endl;
#endif
  }
  catch( FFGraph::Exceptions &eObj ){
    ffexcp = true; ffexcpobj = eObj;
  }
  catch(...){
    galexcp = true;
  }

  // Reset mc::FFVAR_val field to NULL
#ifdef MC__FFUNC_CPU_EVAL
  cputime = -cpuclock();
#endif
  itnVar = nVar.begin(); itpVar = pVar.begin(); itvVar = vVar.begin();
  for( ; itnVar != nVar.end(); ++itnVar, ++itpVar, ++itvVar ){
    for( unsigned i=0; i<(*itnVar); i++ ){
      FFVar* pF = _find_var( (*itpVar)[i].id() );
      if( pF ){
        //std::cerr << "erasing pF->val():" << pF->val() << std::endl;
        pF->reset_val( U() );
      }
    }
  }

  auto ito = opDep.begin();
  for( ; ito!=opDep.end(); ++ito ) (*ito)->flag( false );
  for( unsigned i=0; i<nDep; i++ ){
    if( !pDep[i].ops().first ) continue;
    pDep[i].ops().first->reset_val_subgraph( U() );
  }
#ifdef MC__FFUNC_CPU_EVAL
  cputime += cpuclock();
  std::cout << "Clean-up time: " << std::fixed << cputime << std::endl;
#endif

  if( ffexcp )  throw ffexcpobj;
  if( galexcp ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::EVAL );
  return;
}
/*
template <typename U> inline void
FFGraph::eval
( std::list<const FFOp*>&opDep, const unsigned nDep, const FFVar*pDep,
  U*vDep, const unsigned nVar, const FFVar*pVar, const U*vVar, const bool add )
{
  // Nothing to do!
  if( !nDep ) return;
  assert( pDep && vDep );
  assert( !nVar || ( pVar && vVar ) );
#ifdef MC__FFUNC_CPU_EVAL
  double cputime;
#endif

  // Initialize all independent variables participating in the dependent ones
#ifdef MC__FFUNC_CPU_EVAL
  cputime = -cpuclock();
#endif
  for( unsigned int i=0; i<nVar; i++ ){
    FFVar* pF = _find_var( pVar[i].id() );
    if( pF ){
      pF->val() = new U( vVar[i] );
#ifdef MC__FFUNC_DEBUG_EVAL
      std::cout << pVar[i] << "  " << vVar[i] << std::endl;
#endif
    }
  }
#ifdef MC__FFUNC_CPU_EVAL
  cputime += cpuclock();
  std::cout << "\nIndep. init. time: " << std::fixed << cputime << std::endl;
#endif
  // THIS IS DOING A BIT TOO MUCH WORK AS ONLY THE INDEPENDENT VARIABLES
  // PARTICIPATING IN THE DEPENDENTS SHOULD BE TAKEN INTO ACCOUNT STRICTLY
  // (E.G., THIS COULD BE DONE USING THE .dep() FIELD IN THE DEPENDENTS)

  bool ffexcp = false, galexcp = false;
  FFGraph::Exceptions ffexcpobj;
  try{
#ifdef MC__FFUNC_CPU_EVAL
    cputime = -cpuclock();
#endif
    // Evaluate dependents given by vDep in U type
    typename std::list<const FFOp*>::iterator ito = opDep.begin();
    for( ; ito!=opDep.end(); ++ito ) (*ito)->evaluate( U() );
#ifdef MC__FFUNC_CPU_EVAL
    cputime += cpuclock();
    std::cout << "Evaluation time: " << std::fixed << cputime << std::endl;
#endif

    // Retreive dependents variables in U type as given by vDep
#ifdef MC__FFUNC_CPU_EVAL
    cputime = -cpuclock();
#endif
    for( unsigned i=0; i<nDep; i++ ){
      // Obtain pointer to dependent variable in FFGraph
      FFVar* pF = _find_var( pDep[i].id() );
      // Write/add corresponding evaluation in U type into/to result vector
      if( !add && pF ) vDep[i] = *static_cast<U*>( pF->val() );
      else if( pF )   vDep[i] += *static_cast<U*>( pF->val() );
      else if( !add ) switch( pDep[i].num().t ){
        case FFNum::INT:
          vDep[i] = pDep[i].num().n;
          break;
        case FFNum::REAL:
          vDep[i] = pDep[i].num().x;
          break;
      }
      else switch( pDep[i].num().t ){
        case FFNum::INT:
          vDep[i] += pDep[i].num().n;
          break;
        case FFNum::REAL:
          vDep[i] += pDep[i].num().x;
          break;
      }
    }
#ifdef MC__FFUNC_CPU_EVAL
    cputime += cpuclock();
    std::cout << "Dep. collect. time: " << std::fixed << cputime << std::endl;
#endif
  }
  catch( FFGraph::Exceptions &eObj ){
    ffexcp = true; ffexcpobj = eObj;
  }
  catch(...){
    galexcp = true;
  }

  // Reset mc::FFVAR_val field to NULL
#ifdef MC__FFUNC_CPU_EVAL
  cputime = -cpuclock();
#endif
  for( unsigned int i=0; i<nVar; i++ ){
    FFVar* pF = _find_var( pVar[i].id() );
    if( pF ) pF->reset_val( U() );
  }
  //it_Vars itv = _Vars.begin();
  //for( ; itv!=_Vars.end() && (*itv)->_id.first<=FFVar::VAR; ++itv )
  //  (*itv)->reset_val( U() );
  
  typename std::list<const FFOp*>::iterator ito = opDep.begin();
  for( ; ito!=opDep.end(); ++ito ) (*ito)->flag( false );
  for( unsigned i=0; i<nDep; i++ ){
    if( !pDep[i].ops().first ) continue;
    pDep[i].ops().first->reset_val_subgraph( U() );
  }
#ifdef MC__FFUNC_CPU_EVAL
  cputime += cpuclock();
  std::cout << "Clean-up time: " << std::fixed << cputime << std::endl;
#endif

  if( ffexcp )  throw ffexcpobj;
  if( galexcp ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::EVAL );
  return;
}
*/
template <typename U> inline void
FFGraph::eval
( std::list<const FFOp*>&opDep, U*opRes, const unsigned nDep, const FFVar*pDep,
  U*vDep, const unsigned nVar, const FFVar*pVar, U*vVar, const bool add )
{
  // Nothing to do!
  if( !nDep ) return;
  //assert( opDep.size() && opRes ); <-- Don't test otherwise fails with constant dependents
  assert( pDep && vDep );
  assert( !nVar || ( pVar && vVar ) );
#ifdef MC__FFUNC_CPU_EVAL
  double cputime;
#endif

  // Propagate values in U arithmetic through subgraph
#ifdef MC__FFUNC_CPU_EVAL
  double cputime = -cpuclock();
#endif
  typename std::list<const FFOp*>::iterator ito = opDep.begin();
  for( U*pUres=opRes; ito!=opDep.end(); ++ito, pUres++ ){
    // Initialize variable using values in vVar
    if( (*ito)->type == FFOp::VAR ){
      FFVar* pF = 0;
      for( unsigned i=0; i<nVar; i++ ){
        if( (*ito)->pres->id() == pVar[i].id() ){
          pF = (*ito)->pres; //_find_var( (*iti).first->id() );
          *pUres = vVar[i];
          break;
        }
      }
      if( !pF ) throw typename FFGraph::Exceptions( FFGraph::Exceptions::MISSVAR );
    }
    _curOp = *ito;
    _curOp->evaluate( pUres );
  }
#ifdef MC__FFUNC_CPU_EVAL
  cputime += cpuclock();
  std::cout << "\nEvaluation time: " << std::fixed << cputime << std::endl;
#endif

  // Retreive dependents variables in U type as given by vDep
#ifdef MC__FFUNC_CPU_EVAL
  cputime = -cpuclock();
#endif
  for( unsigned i=0; i<nDep; i++ ){
    // Obtain pointer to dependent variable in FFGraph
    FFVar* pF = _find_var( pDep[i].id() );
    // Write/add corresponding evaluation in U type into/to result vector
    if( !add && pF ) vDep[i] = *static_cast<U*>( pF->val() );
    else if( pF )   vDep[i] += *static_cast<U*>( pF->val() );
    else if( !add ) switch( pDep[i].num().t ){
      case FFNum::INT:
        vDep[i] = pDep[i].num().n;
        break;
      case FFNum::REAL:
        vDep[i] = pDep[i].num().x;
        break;
    }
    else switch( pDep[i].num().t ){
      case FFNum::INT:
        vDep[i] += pDep[i].num().n;
        break;
      case FFNum::REAL:
        vDep[i] += pDep[i].num().x;
        break;
    }
  }
#ifdef MC__FFUNC_CPU_EVAL
  cputime += cpuclock();
  std::cout << "Dep. collect. time: " << std::fixed << cputime << std::endl;
#endif

  return;
}

inline CPPL::dssmatrix
FFGraph::depmap
( const unsigned nDep, const FFVar* const pDep, const unsigned nIndep,
  const FFVar* const pIndep )
{
  // generate subgraph
  std::list<const FFOp*> opDep = subgraph( nDep, pDep );
#ifdef MC__FFUNC_DEBUG_DEPMAP
  output( opDep );
#endif

  return depmap( opDep, nDep, pDep, nIndep, pIndep );
}

inline CPPL::dssmatrix
FFGraph::depmap
( std::list<const FFOp*>&opDep, const unsigned nDep, const FFVar* const pDep,
  const unsigned nIndep, const FFVar* const pIndep )
{
  std::vector<const FFVar*> vDep, vIndep;
  for( unsigned i=0; i<nDep; i++ )   vDep.push_back( pDep+i );
  for( unsigned i=0; i<nIndep; i++ ) vIndep.push_back( pIndep+i );

  return depmap( opDep, vDep, vIndep );
}

inline CPPL::dssmatrix
FFGraph::depmap
( const std::vector<const FFVar*>&vDep, const std::vector<const FFVar*>&vIndep ) 
{
  // generate subgraph
  std::list<const FFOp*> opDep = subgraph( vDep );
#ifdef MC__FFUNC_DEBUG_DEPMAP
  output( opDep );
#endif

  return depmap( opDep, vDep, vIndep );
}

inline CPPL::dssmatrix
FFGraph::depmap
( std::list<const FFOp*>&opDep, const std::vector<const FFVar*>&vDep,
  const std::vector<const FFVar*>&vIndep ) 
{

  bool contains_product = false;
  long ndepmap;
  std::map< const FFVar* , long > aux_map;
  long vindex = 0 , prod_index = 0;

  for (std::vector<const FFVar*>::const_iterator i=vIndep.begin(); i!=vIndep.end(); ++i ){
    aux_map[ _find_var((*i)->_id) ] = vindex; // <-- !!!Consider adding safeguards here!!!
    vindex += 1;
  }
   	
  for(std::list<const FFOp*>::const_iterator i=opDep.begin(); i!=opDep.end(); ++i)
    switch((*i) -> pres ->id().first){
      case FFVar::AUX:
        if( (*i)->type == mc::FFOp::TIMES && !contains_product ){
          contains_product = true;
          prod_index = vindex;
          vindex  += 6;  // Add block for product
        }
        aux_map[ (*i)->pres ] = vindex;
        vindex  += 1; 
        break;

      case FFVar::VAR: case FFVar::CINT: case FFVar::CREAL: default:
        break;
    }

#ifdef MC__FFUNC_DEBUG_DEPMAP
  std::cout << "AUX MAP: " << std::endl;
  for( std::map<const FFVar*,long>::const_iterator i=aux_map.begin(); i!=aux_map.end(); ++i )
    std::cout << *(i->first) << " , " << i->second << std::endl;
#endif

  ndepmap = aux_map.size();
  if( contains_product ) ndepmap += 6; 
  CPPL::dssmatrix depmap(ndepmap);

  //initialize independent variables
  for( unsigned long i=0 ; i<vIndep.size(); ++i ){
    for( unsigned long j=i ; j<vIndep.size(); ++j ) depmap.put( i, j, 1 );
    depmap.put( i, i, 1 );
  }

  //initialize dependent variables
  for( std::vector<const FFVar*>::const_iterator i_it=vDep.begin(); i_it!=vDep.end(); ++i_it ){
    long i = aux_map.find( _find_var((*i_it)->_id) )->second;
    for( std::vector<const FFVar*>::const_iterator j_it = i_it ; j_it != vDep.end(); ++j_it ){
      long j = aux_map.find( _find_var((*j_it)->_id) )->second;
      depmap.put( i, j, 1 );
    }
    depmap.put( i, i, 1);
  }

  long idep = ndepmap - 1;
  contains_product = false;
  for(std::list<const FFOp*>::const_reverse_iterator i_it=opDep.rbegin(); i_it!=opDep.rend(); ++i_it){
    if( (*i_it)->pres->id().first == FFVar::AUX ){   
      if( !(*i_it)->is_univariate() ){
        long k, l;
        if( (*i_it)->plop->id().first == FFVar::AUX || (*i_it)->plop->id().first == FFVar::VAR )
          k = aux_map.find((*i_it)->plop)->second;	
        if( (*i_it)->prop->id().first == FFVar::AUX || (*i_it)->prop->id().first == FFVar::VAR )
          l = aux_map.find((*i_it)->prop)->second;
        else
          l = -1;
        if( k != -1 && l != -1 ) depmap.put(k,l,1);
        for( long j=0; j<idep ; ++j ){
          if( depmap.isListed( idep, j ) ){
            if( k != -1 ) depmap.put( j, k, 1 );
            if( l != -1 ) depmap.put( j, l, 1 ) ;
          }
        }
        if( k != -1 ) depmap.put( k, k ,1 );
        if( l != -1 ) depmap.put( l, l ,1 );
#ifdef MC__FFUNC_DEBUG_DEPMAP
        std::cout << " Row Info: \n";
        std::cout << " Operation : " << *((*i_it)->pres) << " = " <<*(*i_it) << std::endl; 
        std::cout << "i = " << idep << " k = " << k << " l = " << l << std::endl;
#endif
        if( (*i_it)->type == FFOp::TIMES && ( idep == prod_index+6 ) ) 
          idep -= 6; // jump product block if it is the first occurrence of a product
      }
      else{
        long k;
        if( (*i_it)->plop->id().first == FFVar::AUX || (*i_it)->plop->id().first == FFVar::VAR )
          k = aux_map.find((*i_it)->plop)->second;	
        else
          k = -1;
        for( long j=0; j<idep; ++j )
          if( depmap.isListed( idep, j ) && k != -1 ) depmap.put( j, k, 1 );
        if( k != -1 ) depmap.put( k, k, 1 );
#ifdef MC__FFUNC_DEBUG_DEPMAP
        std::cout << " Row Info: \n";
        std::cout << " Operation : " << *((*i_it)->pres) << " = " <<*(*i_it) << std::endl; 
        std::cout << "i = " << idep << " k = " << k << " l = " << -1 << std::endl;
#endif
      }			
      idep -= 1;
    }
  }
  
  return depmap;
}

} // namespace mc

namespace mc
{

//! @brief Specialization of the structure mc::Op to allow usage of the type mc::FFVar as a template parameter in other MC++ types
template <> struct Op< mc::FFVar >
{
  typedef mc::FFVar FV;
  static FV point( const double c ) { return FV(c); }
  static FV zeroone() { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF ); }
  static void I(FV& x, const FV&y) { x = y; }
  static double l(const FV& x) { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF ); }
  static double u(const FV& x) { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF ); }
  static double abs (const FV& x) { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF );  }
  static double mid (const FV& x) { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF );  }
  static double diam(const FV& x) { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF ); }
  static FV inv (const FV& x) { return mc::inv(x);  }
  static FV sqr (const FV& x) { return mc::sqr(x);  }
  static FV sqrt(const FV& x) { return mc::sqrt(x); }
  static FV log (const FV& x) { return mc::log(x);  }
  static FV xlog(const FV& x) { return x*mc::log(x); }
  static FV fabs(const FV& x) { return mc::fabs(x); }
  static FV exp (const FV& x) { return mc::exp(x);  }
  static FV sin (const FV& x) { return mc::sin(x);  }
  static FV cos (const FV& x) { return mc::cos(x);  }
  static FV tan (const FV& x) { return mc::tan(x);  }
  static FV asin(const FV& x) { return mc::asin(x); }
  static FV acos(const FV& x) { return mc::acos(x); }
  static FV atan(const FV& x) { return mc::atan(x); }
  static FV erf (const FV& x) { return mc::erf(x);  }
  static FV erfc(const FV& x) { return mc::erfc(x); }
  static FV fstep(const FV& x) { return mc::fstep(x); }
  static FV bstep(const FV& x) { return mc::bstep(x); }
  static FV hull(const FV& x, const FV& y) { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF ); }
  static FV min (const FV& x, const FV& y) { return mc::min(x,y);  }
  static FV max (const FV& x, const FV& y) { return mc::max(x,y);  }
  static FV arh (const FV& x, const double k) { return mc::exp(-k/x); }
  static FV cheb(const FV& x, const unsigned n) { return mc::cheb(x,n); }
  template <typename X, typename Y> static FV pow(const X& x, const Y& y) { return mc::pow(x,y); }
  static FV monomial (const unsigned int n, const FV* x, const int* k) { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF ); }
  static bool inter(FV& xIy, const FV& x, const FV& y) { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF ); }
  static bool eq(const FV& x, const FV& y) { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF ); }
  static bool ne(const FV& x, const FV& y) { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF ); }
  static bool lt(const FV& x, const FV& y) { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF ); }
  static bool le(const FV& x, const FV& y) { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF ); }
  static bool gt(const FV& x, const FV& y) { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF ); }
  static bool ge(const FV& x, const FV& y) { throw typename FFGraph::Exceptions( FFGraph::Exceptions::UNDEF ); }
};

} // namespace mc

namespace fadbad
{
  //! @brief Specialization of the structure fadbad::Op to allow usage of the type mc::FFVar as a template parameter of the classes fadbad::F, fadbad::B and fadbad::T in FADBAD++
  template <> struct Op<mc::FFVar>{
    typedef double Base;
    typedef mc::FFVar FV;
    static Base myInteger( const int i ) { return Base(i); }
    static Base myZero() { return myInteger(0); }
    static Base myOne() { return myInteger(1);}
    static Base myTwo() { return myInteger(2); }
    static double myPI() { return mc::PI; }
    static FV myPos( const FV& x ) { return  x; }
    static FV myNeg( const FV& x ) { return -x; }
    template <typename U> static FV& myCadd( FV& x, const U& y ) { return x+=y; }
    template <typename U> static FV& myCsub( FV& x, const U& y ) { return x-=y; }
    template <typename U> static FV& myCmul( FV& x, const U& y ) { return x*=y; }
    template <typename U> static FV& myCdiv( FV& x, const U& y ) { return x/=y; }
    static FV myInv( const FV& x ) { return mc::inv( x ); }
    static FV mySqr( const FV& x ) { return mc::pow( x, 2 ); }
    template <typename X, typename Y> static FV myPow( const X& x, const Y& y ) { return mc::pow( x, y ); }
    static FV myCheb( const FV& x, const unsigned n ) { return mc::cheb( x, n ); }
    static FV mySqrt( const FV& x ) { return mc::sqrt( x ); }
    static FV myLog( const FV& x ) { return mc::log( x ); }
    static FV myExp( const FV& x ) { return mc::exp( x ); }
    static FV mySin( const FV& x ) { return mc::sin( x ); }
    static FV myCos( const FV& x ) { return mc::cos( x ); }
    static FV myTan( const FV& x ) { return mc::tan( x ); }
    static FV myAsin( const FV& x ) { return mc::asin( x ); }
    static FV myAcos( const FV& x ) { return mc::acos( x ); }
    static FV myAtan( const FV& x ) { return mc::atan( x ); }
    static bool myEq( const FV& x, const FV& y ) { throw std::runtime_error("fadbad::Op<FFVar>::myEq -- operation not permitted"); }
    static bool myNe( const FV& x, const FV& y ) { throw std::runtime_error("fadbad::Op<FFVar>::myNe -- operation not permitted"); }
    static bool myLt( const FV& x, const FV& y ) { throw std::runtime_error("fadbad::Op<FFVar>::myLt -- operation not permitted"); }
    static bool myLe( const FV& x, const FV& y ) { throw std::runtime_error("fadbad::Op<FFVar>::myLe -- operation not permitted"); }
    static bool myGt( const FV& x, const FV& y ) { throw std::runtime_error("fadbad::Op<FFVar>::myGt -- operation not permitted"); }
    static bool myGe( const FV& x, const FV& y ) { throw std::runtime_error("fadbad::Op<FFVar>::myGe -- operation not permitted"); }
  };
} // end namespace fadbad

#endif
