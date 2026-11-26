// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

// Simulation of a surrogate of the Ackley function of two variables on [0,1]^2, saved as a dictionary by
// torch.save: a piecewise-linear Kolmogorov-Arnold network (ackley2d_KAN-*.pt) or a dense ReLU network
// (ackley2d_ANN-*.pt).
//  1. The file is read once by MC++, without Torch (mc::NETIO::read_dict), and the dictionary is given to an
//     mc::KAN or an mc::MLP (set_data), which maps the raw inputs to the raw outputs: the scalings of the
//     file are folded in.
//  2. The network is inserted in a DAG as an FFKAN or FFMLP operation.
//  3. The DAG is evaluated at N random points of [0,1]^2, written to <file>.out with the values of the
//     Ackley function: x1, x2, network, Ackley. For gnuplot:
//       gnuplot -p -e 'file="ackley2d_KAN-TITO-h5-g8x8.out"' sim_NN.gnu
// notebook/sim_NN.ipynb does the same in Python.
//
// Usage: ./sim_NN <file.pt> [N=1000]

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "interval.hpp"
#include "ffkan.hpp"
#include "ffmlp.hpp"

typedef mc::Interval I;

// Ackley function of two variables
double ackley( double const x1, double const x2 )
{
  return -20. * std::exp( -0.2 * std::sqrt( 0.5 * ( x1*x1 + x2*x2 ) ) )
         - std::exp( 0.5 * ( std::cos( 2.*mc::PI*x1 ) + std::cos( 2.*mc::PI*x2 ) ) ) + std::exp( 1. ) + 20.;
}

int main( int argc, char** argv )
{
  if( argc < 2 ){ std::cerr << "Usage: " << argv[0] << " <file.pt> [N=1000]" << std::endl; return 1; }
  std::string const file = argv[1];
  size_t const N = ( argc > 2? std::stoul( argv[2] ): 1000 );

  // Network: the file is read once, then given to a KAN or to a dense network
  mc::NETIO::Dict D;
  std::string msg;
  if( !mc::NETIO::read_dict( file, D, msg ) ){ std::cerr << file << ": " << msg << std::endl; return 1; }
  bool const isKAN = ( D.type == "KANmodel" );
  mc::KAN<I> KAN;
  mc::MLP<I> ANN;
  if( !( isKAN? KAN.set_data( D, &msg ): ANN.set_data( D, &msg ) ) ){ std::cerr << file << ": " << msg << std::endl; return 1; }
  std::cout << file << ": " << D.type << ", widths";
  for( size_t n : D.widths ) std::cout << " " << n;
  if( isKAN ) std::cout << ", " << KAN.nwei() << " segments";
  std::cout << std::endl;

  // Network in a DAG
  mc::FFGraph DAG;
  std::vector<mc::FFVar> X = DAG.add_vars( 2, "X" );
  mc::FFKAN<I> OpKAN;
  mc::FFMLP<I> OpANN;
  std::vector<mc::FFVar> Z = ( isKAN? OpKAN( X, &KAN ): OpANN( X, &ANN ) );

  // Simulation at N random points of [0,1]^2
  std::mt19937 gen( 0 );
  std::uniform_real_distribution<double> U( 0., 1. );
  std::vector<std::vector<double>> xs( N, std::vector<double>( 2 ) ), zs;
  for( auto& x : xs ) x = { U( gen ), U( gen ) };
  DAG.veval( Z, zs, X, xs );

  // Network and Ackley values
  std::string const fout = file.substr( 0, file.rfind( '.' ) ) + ".out";
  std::ofstream os( fout );
  os.precision( 10 );
  double sqr = 0., emax = 0.;
  for( size_t p=0; p<N; p++ ){
    double const f = ackley( xs[p][0], xs[p][1] ), e = zs[p][0] - f;
    os << xs[p][0] << " " << xs[p][1] << " " << zs[p][0] << " " << f << std::endl;
    sqr += e*e;
    emax = std::max( emax, std::fabs( e ) );
  }
  std::cout << "  " << N << " points written to " << fout << ": x1, x2, network, Ackley" << std::endl
            << "  error of the network: RMS " << std::sqrt( sqr / N ) << ", largest " << emax << std::endl;
  return 0;
}
