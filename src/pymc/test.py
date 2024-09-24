import pymc

def dag_test1():

  # Define DAG environment
  DAG = pymc.FFGraph()

  # Define variables and dependents
  X = pymc.FFVar(DAG,"X")
  Y = pymc.FFVar(DAG,"Y")
  C = pymc.FFVar(3)
  F = pymc.exp(X*Y)-2*X**2+C
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
  print( "grad F @([0,1],[1,2]): ", DAG.eval( SGDFDXY, DFDXY[2], [X,Y], [pymc.Interval(0.,1.),pymc.Interval(1.,2.)] ) )

  # Setting certain variables
  Y.set(2.);
  print( "grad F @(1,2): ", DAG.eval( DFDXY[2], [X], [1.] ) )
  print( "grad F @([0,1],2): ", DAG.eval( DFDXY[2], [X], [pymc.Interval(0.,1.)] ) )
  Y.unset()


def dag_test2():

  # Define DAG environment
  DAG = pymc.FFGraph()

  # Define variables and dependents
  X = pymc.FFVar(DAG,"X")
  Y = pymc.FFVar(DAG,"Y")
  Z = pymc.FFVar(DAG,"Z")
  C = pymc.FFVar(3)
  F = Z-pymc.exp(X*Y)-2*X**2+C
  G = pymc.exp(Y*Z)
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
  DAG = pymc.FFGraph()

  # Define variables and dependents
  X = pymc.FFVar(DAG,"X")
  Y = pymc.FFVar(DAG,"Y")
  F = pymc.pow(X,2)+X*Y+4

  # Evaluation in interval arithmetic
  IX, IY = pymc.Interval(-0.8,-0.3), pymc.Interval(6.,9.)
  [IF] = DAG.eval( [F], [X,Y], [IX,IY] )
  print( "IX: ", IX, "IY: ", IY, "IF: ", IF )
  
  # Constraint propagation in interval arithmetic
  IF = pymc.Interval(0.)
  [IX,IY], [IF] = DAG.reval( [F], [IF], [X,Y], [IX,IY], pymc.Interval(-1,1)*1e20 )
  print( "IX: ", IX, "IY: ", IY, "IF: ", IF )


def dag_test4():

  # Define DAG environment
  DAG = pymc.FFGraph()

  # Define variables and dependents
  X = pymc.FFVar(DAG,"X")
  Y = pymc.FFVar(DAG,"Y")
  C = pymc.FFVar(3)
  F = pymc.exp(X*Y)-2*X**2+C
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
