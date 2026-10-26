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
    # The record must come from the wheels' toolchain, and the one thing a local build gets wrong without
    # noticing is the pybind11 it was linked against: the 5.1 record was first taken with 3.0.1 from the
    # machine while the wheels build with pyproject's pinned 3.0.4, and CI failed on that line alone
    # (2026-10-26).  Refuse to record a pybind11 other than the pin.
    pin = re.search( r'"pybind11==([0-9]+\.[0-9]+\.[0-9]+)"', open( os.path.join( ROOT, "pyproject.toml" ) ).read() )
    built = re.search( r'^pybind11 ([0-9]+\.[0-9]+\.[0-9]+)', got, re.M )
    if pin and built and pin.group( 1 ) != built.group( 1 ):
      print( "check_abi: refusing to record: the fingerprint was built against pybind11 %s, the wheels use the pinned "
             "%s (pyproject.toml).  The usual cause is PyTorch, whose include directory carries its own pybind11 and "
             "shadows extern/pybind11 when ENABLE_TORCH is on.  Record in the wheels' configuration, in a build "
             "directory of its own:\n"
             "    cmake -S . -B build-abi -DENABLE_TORCH=OFF -DENABLE_HSL=OFF -DPYMCPP_VENDORED_PYBIND11=ON\n"
             "    cmake --build build-abi --target abi_fingerprint\n"
             "    python3 tools/abi/check_abi.py build-abi/abi_fingerprint --record"
             % ( built.group( 1 ), pin.group( 1 ) ) ); return 1
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
  # The pybind11 line holds a version and an INTERNALS number.  What decides whether two modules can exchange
  # objects is the internals number -- pybind11 keeps it fixed across patch releases and bumps it when the
  # internals ABI changes -- so that is what binds.  A different version with the same internals is noted,
  # not failed: a Torch-enabled build compiles the fingerprint against PyTorch's bundled pybind11 (3.0.1 in
  # torch 2.x) while the wheels use the pinned 3.0.4, and every class layout agrees (2026-10-26).
  pb = re.compile( r'^pybind11 (\S+) internals (\d+)$', re.M )
  mw, mg = pb.search( want ), pb.search( got )
  note = ""
  if mw and mg and mw.group( 2 ) == mg.group( 2 ) and mw.group( 1 ) != mg.group( 1 ):
    note = ( "  note: built against pybind11 %s, recorded with %s -- the same internals (%s), so the ABI is the same;\n"
             "  the record is taken in the wheels' configuration (INSTALL.md, \"Releasing\")" % ( mg.group( 1 ), mw.group( 1 ), mg.group( 2 ) ) )
    got = got.replace( mg.group( 0 ), mw.group( 0 ) )
  if got == want:
    print( "check_abi: OK -- MC++ %s keeps the %d.%d fingerprint (%s)" % ( version, major, minor, os.path.basename( path ) ) )
    if note: print( note )
    return 0
  print( "check_abi: FAIL -- MC++ %s changes the fingerprint recorded for %d.%d (%s):" % ( version, major, minor, os.path.basename( path ) ) )
  for line in difflib.unified_diff( want.splitlines(), got.splitlines(), "recorded", "built", lineterm = "" ):
    print( "  " + line )
  print( "A patch release must keep the layout of the classes modules built on pymcpp share with it, and the pybind11\n"
         "internals version (INSTALL.md, \"Releasing\").  Release this as %d.%d.0 and record its fingerprint (--record).  If the\n"
         "sources did not change these classes, the CI toolchain did (its compiler, standard library or Armadillo):\n"
         "re-record with the same minor version only after checking that." % ( major, minor + 1 ) )
  return 1

if __name__ == "__main__":
  sys.exit( main( sys.argv[1:] ) )
