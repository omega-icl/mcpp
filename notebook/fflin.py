import pymc
G = pymc.FFGraph()
X = G.add_vars( 3 )
Sum = pymc.FFLin()

Sum1 = Sum( X, 1. )
print( Sum1.str() )

Sum1 = Sum( X, 1., 2. )
print( Sum1.str() )

Sum1 = Sum( X, 2. )
print( Sum1.str() )

Sum1 = Sum( X, 2., 1. )
print( Sum1.str() )

print( G )

Sum1 = Sum( X, [-1,-1,-1] )
print( Sum1.str() )

print( G )

Sum1 = Sum( X, -1. );
print( Sum1.str() )

print( G )

