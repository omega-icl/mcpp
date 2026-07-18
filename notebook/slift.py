import pymcpp

DAG = pymcpp.FFGraph()
X = DAG.add_vars( 2, "X" )
F = X[0]**3 + pymcpp.sqrt( X[0]**2 + X[1]**2 )
print( F.str() )

SE = pymcpp.SLift( DAG )
SE.process( [F], True ) #False )
print( SE )
#print( SE.aux )
#print( SE.aux[0][0].str() )
#print( SE.aux[1][0].str() )
print( SE.OpLift )
print( SE.AuxLift )

mon1 = pymcpp.FFMon( { X[0]: 2, X[1]: 1 } )
mon2 = pymcpp.FFMon( { X[0]: 1, X[1]: 1 } )
print( mon1.display(), mon2.display() )

pymcpp.FFPoly.options.BASIS = pymcpp.FFPoly.options.MONOM;
pol1 = pymcpp.FFPoly( { mon1: 2, mon2: -1 } )
print( pol1, " at (-0.5,1.5) = ", pol1.eval( { X[0]: -0.5, X[1]: 1.5 } ) )
print( pol1.factor( X[0] ) )
print( pol1.factor( X[1] ) )

pymcpp.FFPoly.options.BASIS = pymcpp.FFPoly.options.CHEB;
PX = [pymcpp.FFPoly( X[0] ), pymcpp.FFPoly( X[1] )]
PF = PX[0]**2 * PX[1]**3
PF += ( pymcpp.FFMon(X[0],2), 2. )
print( PF )
print( PF.coefmon )

PF.convert( pymcpp.FFPoly.options.MONOM )
pymcpp.FFPoly.options.BASIS = pymcpp.FFPoly.options.MONOM;
print( PF )
