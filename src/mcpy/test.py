import mcpy

def dag_test():

  # Define DAG environment
  DAG = mcpy.FFGraph()

  # Define variables and dependents
  X = mcpy.FFVar(DAG,"X")
  Y = mcpy.FFVar(DAG,"Y")
  C = mcpy.FFVar(3)
  F = mcpy.exp(X*Y)-2*X**2+C
  F.set( "F" )
  
  # Subgraph and dot script
  SGF = DAG.subgraph( [F] )
  DAG.output( SGF )
  DAG.dot_script( [F], "F.dot" )

  # Backward differentiation
  DFDXY = DAG.bdiff( [F], [X,Y] )
  print( DFDXY )
  SGDFDXY = DAG.subgraph( DFDXY[2] )
  DAG.output( SGDFDXY )

  # Forward second-order differentiation
  D2FDXY2 = DAG.fdiff( DFDXY[2], [X,Y] )
  print( D2FDXY2 )
  SGD2FDXY2 = DAG.subgraph( D2FDXY2[2] )
  DAG.output( SGD2FDXY2 )

  # Dependent evaluation in various arithmetics
  print( "grad F @(1,1): ", DAG.eval( SGDFDXY, DFDXY[2], [X,Y], [1.,1.] ) )
  print( "grad F @([0,1],[1,2]): ", DAG.eval( SGDFDXY, DFDXY[2], [X,Y], [mcpy.Interval(0.,1.),mcpy.Interval(1.,2.)] ) )

  # Setting certain variables
  Y.set(2.);
  print( "grad F @(1,2): ", DAG.eval( DFDXY[2], [X], [1.] ) )
  print( "grad F @([0,1],2): ", DAG.eval( DFDXY[2], [X], [mcpy.Interval(0.,1.)] ) )
  Y.unset()
  
dag_test()

