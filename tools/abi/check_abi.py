#!/usr/bin/env python3
# check_abi.py -- enforce MC++'s patch-release rule (INSTALL.md, "Releasing") on the ABI fingerprint.
#
#   check_abi.py FINGERPRINT_BINARY [--version X.Y.Z] [--record]
#
# Runs the abi_fingerprint program (cmake --build build --target abi_fingerprint) and compares its output with
# tools/abi/fingerprint-X.Y-<platform>.txt, the fingerprint recorded for this MINOR version on this platform:
#   * identical                         -> OK (exit 0)
#   * different                         -> FAIL (exit 1): a patch release must not change the layout of the classes that
#                                          modules built on pymcpp share with it, nor the pybind11 version -- release
#                                          a new minor version and record its fingerprint
#   * none recorded, X.Y.0              -> FAIL, asking to record it (--record) with the new minor version
#   * none recorded for this platform   -> skipped with a note (only some platforms are recorded; a source change
#                                          shows on any of them)
# --record writes the file: allowed for X.Y.0, or when none is recorded yet for X.Y on this platform.
# The version defaults to pyproject.toml's.
import os, platform, re, subprocess, sys, difflib

HERE = os.path.dirname( os.path.abspath( __file__ ) )
ROOT = os.path.dirname( os.path.dirname( HERE ) )

def main( argv ):
  args = [ a for a in argv if not a.startswith( "--" ) ]
  record = "--record" in argv
  version = None
  if "--version" in argv:
    version = argv[ argv.index( "--version" ) + 1 ]
    args = [ a for a in args if a != version ]
  if len( args ) != 1:
    print( __doc__ or "usage: check_abi.py FINGERPRINT_BINARY [--version X.Y.Z] [--record]" ); return 2
  if version is None:
    m = re.search( r'^version\s*=\s*"([0-9]+\.[0-9]+\.[0-9]+)"', open( os.path.join( ROOT, "pyproject.toml" ) ).read(), re.M )
    if not m: print( "check_abi: no version in pyproject.toml" ); return 2
    version = m.group( 1 )
  major, minor, patch = ( int( x ) for x in version.split( "." ) )
  plat = "%s-%s" % ( sys.platform.replace( "darwin", "macos" ).replace( "win32", "windows" ), platform.machine().lower() )
  path = os.path.join( HERE, "fingerprint-%d.%d-%s.txt" % ( major, minor, plat ) )
  got = subprocess.run( [ args[0] ], capture_output = True, text = True, check = True ).stdout

  if record:
    if os.path.exists( path ) and patch > 0 and open( path ).read() != got:
      print( "check_abi: refusing to re-record %s at patch release %s: the fingerprint changed -- that needs a new "
             "minor version" % ( os.path.basename( path ), version ) ); return 1
    open( path, "w" ).write( got )
    print( "check_abi: recorded %s (MC++ %s)" % ( os.path.basename( path ), version ) ); return 0

  if not os.path.exists( path ):
    recorded = [ f for f in os.listdir( HERE ) if f.startswith( "fingerprint-%d.%d-" % ( major, minor ) ) ]
    if recorded:
      print( "check_abi: no fingerprint recorded for %s on %s (recorded: %s) -- skipped" % ( version, plat, ", ".join( recorded ) ) )
      return 0
    print( "check_abi: no fingerprint recorded for MC++ %d.%d.  %s" % ( major, minor,
           "A new minor version: record it with --record." if patch == 0 else
           "A patch release needs its minor version's fingerprint: record it at %d.%d.0." % ( major, minor ) ) )
    return 1

  want = open( path ).read()
  if got == want:
    print( "check_abi: OK -- MC++ %s keeps the %d.%d fingerprint (%s)" % ( version, major, minor, os.path.basename( path ) ) )
    return 0
  print( "check_abi: FAIL -- MC++ %s changes the fingerprint recorded for %d.%d (%s):" % ( version, major, minor, os.path.basename( path ) ) )
  for line in difflib.unified_diff( want.splitlines(), got.splitlines(), "recorded", "built", lineterm = "" ):
    print( "  " + line )
  print( "A patch release must keep the layout of the classes modules built on pymcpp share with it, and the pybind11\n"
         "version (INSTALL.md, \"Releasing\").  Release this as %d.%d.0 and record its fingerprint (--record).  If the\n"
         "sources did not change these classes, the CI toolchain did (its compiler, standard library or Armadillo):\n"
         "re-record with the same minor version only after checking that." % ( major, minor + 1 ) )
  return 1

if __name__ == "__main__":
  sys.exit( main( sys.argv[1:] ) )
