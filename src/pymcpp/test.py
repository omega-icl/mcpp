import pymcpp

def dag_test1():

  # Define DAG environment
  DAG = pymcpp.FFGraph()

  # Define variables and dependents
  X = pymcpp.FFVar(DAG,"X")
  Y = pymcpp.FFVar(DAG,"Y")
  C = pymcpp.FFVar(3)
  F = pymcpp.exp(X*Y)-2*X**2+C
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
  print( "grad F @([0,1],[1,2]): ", DAG.eval( SGDFDXY, DFDXY[2], [X,Y], [pymcpp.Interval(0.,1.),pymcpp.Interval(1.,2.)] ) )

  # Setting certain variables
  Y.set(2.);
  print( "grad F @(1,2): ", DAG.eval( DFDXY[2], [X], [1.] ) )
  print( "grad F @([0,1],2): ", DAG.eval( DFDXY[2], [X], [pymcpp.Interval(0.,1.)] ) )
  Y.unset()


def dag_test2():

  # Define DAG environment
  DAG = pymcpp.FFGraph()

  # Define variables and dependents
  X = pymcpp.FFVar(DAG,"X")
  Y = pymcpp.FFVar(DAG,"Y")
  Z = pymcpp.FFVar(DAG,"Z")
  C = pymcpp.FFVar(3)
  F = Z-pymcpp.exp(X*Y)-2*X**2+C
  G = pymcpp.exp(Y*Z)
  F.set( "F" )
  G.set( "G" )

  # Subgraph and dot script
  SG = DAG.subgraph( [F,G] )
  DAG.output( SG )
  DAG.dot_script( [F,G], "FG.dot" )

  # Backward differentiation
  DFGDXYZ = DAG.fdiff( [F,G], [X,Y,Z] )
  print( DFGDXYZ )

  # Forward second-order differentiation
  D2FGDXYZ2 = DAG.fdiff( DFGDXYZ[2], [X,Y,Z] )
  print( D2FGDXYZ2 )


def dag_test3():

  # Define DAG environment
  DAG = pymcpp.FFGraph()

  # Define variables and dependents
  X = pymcpp.FFVar(DAG,"X")
  Y = pymcpp.FFVar(DAG,"Y")
  F = pymcpp.pow(X,2)+X*Y+4

  # Evaluation in interval arithmetic
  IX, IY = pymcpp.Interval(-0.8,-0.3), pymcpp.Interval(6.,9.)
  [IF] = DAG.eval( [F], [X,Y], [IX,IY] )
  print( "IX: ", IX, "IY: ", IY, "IF: ", IF )

  # Constraint propagation in interval arithmetic
  IF = pymcpp.Interval(0.)
  [IX,IY], [IF] = DAG.reval( [F], [IF], [X,Y], [IX,IY], pymcpp.Interval(-1,1)*1e20 )
  print( "IX: ", IX, "IY: ", IY, "IF: ", IF )


def dag_test4():

  # Define DAG environment
  DAG = pymcpp.FFGraph()

  # Define variables and dependents
  X = pymcpp.FFVar(DAG,"X")
  Y = pymcpp.FFVar(DAG,"Y")
  C = pymcpp.FFVar(3)
  F = pymcpp.exp(X*Y)-2*X**2+C
  F.set( "F" )

  # Subgraph and dot script
  SGF = DAG.subgraph( [F] )
  DAG.output( SGF )

  # Dependent evaluation in various arithmetics
  print( DAG.veval( SGF, [F], [X,Y], [[1.,1.],[2.,2.]] ) )

#dag_test1()
#dag_test2()
#dag_test3()
dag_test4()
