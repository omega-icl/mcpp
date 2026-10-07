#!/usr/bin/env python3
# Fail if a repaired wheel bundles a library that THIRD_PARTY_NOTICES.txt does not cover (2026-10-07).
#   python tools/ci/check_notices.py THIRD_PARTY_NOTICES.txt WHEEL_OR_DIRECTORY [...]
# The libraries a wheel bundles depend on the platform and the build image (auditwheel, delocate and delvewheel copy
# whatever the module links), so the notices are checked against the wheels themselves.  A directory argument is
# searched for *.whl here, not by the shell, so the same command works under bash and PowerShell.  Bundled libraries
# are found wherever the repair tools put them -- <pkg>.libs/ (auditwheel), <pkg>.dylibs/ (delocate),
# <pkg>-<ver>.data/platlib/ or <pkg>.libs/ (delvewheel) -- as any shared library below the wheel's root.
import glob, os, re, sys, zipfile

# bundled library (normalised name prefix) -> the "Name:" of its notice block
COVER = {
  "libarmadillo": "Armadillo", "libarpack": "ARPACK-NG", "libopenblas": "OpenBLAS", "libblas": "LAPACK",
  "liblapack": "LAPACK", "libsatlas": "ATLAS", "libtatlas": "ATLAS",
  "libgfortran": "GCC runtime", "libquadmath": "GCC runtime", "libgomp": "GCC runtime", "libgcc_s": "GCC runtime",
  "libomp": "LLVM OpenMP runtime",
  "msvcp140": "Microsoft Visual C++ runtime", "vcruntime140": "Microsoft Visual C++ runtime",
  "vcomp140": "Microsoft Visual C++ runtime", "concrt140": "Microsoft Visual C++ runtime",
  "libsuperlu": "SuperLU", "libsundials": "SUNDIALS", "libsuitesparseconfig": "SuiteSparse SuiteSparse_config",
  "libamd": "SuiteSparse AMD", "libcamd": "SuiteSparse CAMD", "libcolamd": "SuiteSparse COLAMD",
  "libccolamd": "SuiteSparse CCOLAMD", "libbtf": "SuiteSparse BTF", "libklu": "SuiteSparse KLU",
  "libcholmod": "SuiteSparse CHOLMOD", "libumfpack": "SuiteSparse UMFPACK", "libspqr": "SuiteSparse SPQR",
}
LIB = re.compile( r"\.(dll|dylib)$|\.so(\.\d+)*$", re.I )

def bundled( wheel ):
    out = set()
    for n in zipfile.ZipFile( wheel ).namelist():
        if "/" not in n or n.endswith( "/" ) or ".dist-info/" in n: continue      # the module itself sits at the root
        base = n.rsplit( "/", 1 )[ -1 ]
        if not LIB.search( base ): continue
        name = re.sub( r"-([0-9a-f]{32}|[0-9a-f]{8})(?=\.)", "", base )             # delvewheel's or auditwheel's hash
        name = re.split( r"[.]", name, 1 )[ 0 ]                                     # extension and version
        name = re.sub( r"-[\d.]+$", "", name ).lower()                              # e.g. libgcc_s_seh-1
        out.add( name )
    return sorted( out )

def main():
    notices, args = sys.argv[ 1 ], sys.argv[ 2: ]
    wheels = []
    for a in args:
        wheels += sorted( glob.glob( os.path.join( a, "*.whl" ) ) ) if os.path.isdir( a ) else [ a ]
    if not wheels:
        print( "check_notices.py: no wheel found in %s" % args ); return 1
    names = re.findall( r"^Name: (.*)$", open( notices, encoding = "utf-8" ).read(), flags = re.M )
    bad = []
    for w in wheels:
        libs = bundled( w ); short = os.path.basename( w )
        for lib in libs:
            key = next( ( k for k in COVER if lib.startswith( k ) ), None )
            if key is None:
                bad.append( "%s: bundles %s, which has no entry in this script's table" % ( short, lib ) )
            elif not any( n.startswith( COVER[ key ] ) for n in names ):
                bad.append( "%s: bundles %s, but '%s' has no notice in %s" % ( short, lib, COVER[ key ], notices ) )
        print( "%s: %d bundled libraries checked (%s)" % ( short, len( libs ), " ".join( libs ) ) )
    for b in bad: print( "  MISSING " + b )
    return 1 if bad else 0

sys.exit( main() )
