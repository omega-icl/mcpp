#!/usr/bin/env python3
# Fail if a repaired wheel bundles a library that THIRD_PARTY_NOTICES.txt does not cover (2026-10-07).
#   python tools/ci/check_notices.py THIRD_PARTY_NOTICES.txt wheelhouse/*.whl
# The libraries a wheel bundles depend on the image it was built on (auditwheel/delocate copy whatever the module
# links), so the notices are checked against the wheels themselves rather than a list kept by hand.
import re, sys, zipfile
# bundled library (file name prefix) -> the "Name:" of its notice block
COVER = {
  "libarmadillo": "Armadillo", "libarpack": "ARPACK-NG", "libopenblas": "OpenBLAS", "libblas": "LAPACK",
  "liblapack": "LAPACK", "libgfortran": "GCC runtime", "libquadmath": "GCC runtime", "libgomp": "GCC runtime", "libsatlas": "ATLAS", "libtatlas": "ATLAS",
  "libsuperlu": "SuperLU", "libsundials": "SUNDIALS", "libsuitesparseconfig": "SuiteSparse SuiteSparse_config",
  "libamd": "SuiteSparse AMD", "libcamd": "SuiteSparse CAMD", "libcolamd": "SuiteSparse COLAMD",
  "libccolamd": "SuiteSparse CCOLAMD", "libbtf": "SuiteSparse BTF", "libklu": "SuiteSparse KLU",
  "libcholmod": "SuiteSparse CHOLMOD", "libumfpack": "SuiteSparse UMFPACK", "libspqr": "SuiteSparse SPQR",
}
def main():
    notices, wheels = sys.argv[ 1 ], sys.argv[ 2: ]
    names = re.findall( r"^Name: (.*)$", open( notices ).read(), flags = re.M )
    bad = []
    for w in wheels:
        libs = sorted( { re.split( r"[-.]", n.rsplit( "/", 1 )[ -1 ], 1 )[ 0 ] for n in zipfile.ZipFile( w ).namelist()
                         if re.search( r"\.(libs|dylibs)/", n ) } )
        for lib in libs:
            if not lib: continue                       # the .libs/ directory entry itself
            key = next( ( k for k in COVER if lib.startswith( k ) ), None )
            if key is None:
                bad.append( "%s: bundles %s, which has no entry in this script's table" % ( w.rsplit( "/", 1 )[ -1 ], lib ) )
            elif not any( n.startswith( COVER[ key ] ) for n in names ):
                bad.append( "%s: bundles %s, but %s has no notice in %s" % ( w.rsplit( "/", 1 )[ -1 ], lib, COVER[ key ], notices ) )
        print( "%s: %d bundled libraries checked" % ( w.rsplit( "/", 1 )[ -1 ], len( libs ) ) )
    for b in bad: print( "  MISSING " + b )
    return 1 if bad else 0
sys.exit( main() )
