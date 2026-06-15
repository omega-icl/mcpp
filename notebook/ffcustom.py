import pymc

G = pymc.FFGraph()
X = G.add_vars( 2 )

def prod( x ):
  return [ x[0]*x[1] ]

def dprod( x ):
  return [ x[1], x[0] ]

OpDY = pymc.FFCustom()
OpDY.set_D_eval( dprod )

OpY = pymc.FFCustom()
OpY.set_D_eval( prod )
OpY.set_I_eval( prod )
OpY.set_deriv( OpDY, 0 )

Y = OpY( X, 1 )
print( Y, " = ", Y.str() )

print( G )

[DY] = G.eval( [Y], X, [2,3] )
print( "DY = ", DY )

[IY] = G.eval( [Y], X, [pymc.Interval(1,3), pymc.Interval(2,4)] )
print( "IY = ", IY )

dYdX = G.fdiff( [Y], X )
[ print( dYdXi, " = ", dYdXi.str() ) for dYdXi in dYdX[2] ]

dDYdX = G.eval( dYdX[2], X, [2,3] )
print( dDYdX )

print( G )
