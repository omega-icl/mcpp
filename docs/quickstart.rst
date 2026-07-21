Quickstart
==========

Building an expression DAG
--------------------------

Functions are defined on a directed acyclic graph (DAG) environment
``FFGraph``, with participating variables of type ``FFVar``:

.. code-block:: python

   import pymcpp

   DAG = pymcpp.FFGraph()
   X = pymcpp.FFVar(DAG, "X")
   Y = pymcpp.FFVar(DAG, "Y")
   F = pymcpp.exp(X * Y) - 2 * X**2 + 3
   F.set("F")

   SGF = DAG.subgraph([F])
   DAG.output(SGF)          # print the operations in the subgraph

Automatic differentiation
-------------------------

Derivatives are generated symbolically on the DAG, by backward
(``bdiff``) or forward (``fdiff``) automatic differentiation:

.. code-block:: python

   DFDXY = DAG.bdiff([F], [X, Y])     # gradient of F w.r.t. X and Y
   SGD = DAG.subgraph(DFDXY[2])

Evaluation in different arithmetics
-----------------------------------

The same subgraph can be evaluated in real arithmetic or in any of the
bounding arithmetics, simply by passing values of the corresponding type:

.. code-block:: python

   # Real arithmetic
   print(DAG.eval(SGD, DFDXY[2], [X, Y], [1.0, 1.0]))

   # Interval arithmetic: rigorous bounds on the gradient over a box
   IX = pymcpp.Interval(0.0, 1.0)
   IY = pymcpp.Interval(1.0, 2.0)
   print(DAG.eval(SGD, DFDXY[2], [X, Y], [IX, IY]))

Vectorized evaluation over several points is available via ``veval``:

.. code-block:: python

   print(DAG.veval(SGF, [F], [X, Y], [[1.0, 1.0], [2.0, 2.0]]))

Where to go next
----------------

Each bounding arithmetic is illustrated in a dedicated notebook in
:doc:`tutorials`. The full list of classes and functions is in the
:doc:`api`.
