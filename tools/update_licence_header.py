#!/usr/bin/env python3
# Copyright (C) Benoit Chachuat, Imperial College London.
# All Rights Reserved.
# This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
# SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
#
# update_licence_header.py -- state the licence (EPL-2.0 with GPL-2.0-or-later as a Secondary License) in the header of
# every .hpp and .cpp file of the project (2026-10-07).
#
#   python tools/update_licence_header.py [ROOT] [--dry-run | --check]
#
# ROOT defaults to the repository holding this script.  Only files whose header carries the project's copyright line
# are changed, so vendored third-party code keeps its own licence; the directories .git, extern, src/3rdparty and
# build* are not searched at all.  In such a file, the licence line
#     // This code is published under the Eclipse Public License.        (or "... License 2.0.")
# becomes
#     // This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
#     // SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
# Nothing else is touched (line endings preserved).  Running it again changes nothing.  A file with the copyright line
# but no recognisable licence line is REPORTED and left alone, rather than guessed at.
#   --dry-run   list what would change, write nothing
#   --check     exit 1 if any file still needs the update (for CI)
import os, re, sys

COPYRIGHT = "// Copyright (C) Benoit Chachuat, Imperial College London."
OLD = re.compile( r"^// This code is published under the Eclipse Public License( 2\.0)?\.[ \t]*$" )
NEW = [ "// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.",
        "// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later" ]
SKIP_DIRS = { ".git", "extern" }
HEADER_LINES = 6                     # the copyright and licence lines sit at the top of the file

def candidates( root ):
    for d, dirs, files in os.walk( root ):
        rel = os.path.relpath( d, root )
        dirs[:] = sorted( x for x in dirs
                          if x not in SKIP_DIRS and not x.startswith( "build" )
                          and os.path.normpath( os.path.join( rel, x ) ) != os.path.normpath( "src/3rdparty" ) )
        for f in sorted( files ):
            if f.endswith( ( ".hpp", ".cpp" ) ):
                yield os.path.join( d, f )

def main():
    args = [ a for a in sys.argv[ 1: ] if not a.startswith( "--" ) ]
    flags = { a for a in sys.argv[ 1: ] if a.startswith( "--" ) }
    if flags - { "--dry-run", "--check" } or len( args ) > 1:
        sys.exit( "usage: update_licence_header.py [ROOT] [--dry-run | --check]" )
    root = args[ 0 ] if args else os.path.dirname( os.path.dirname( os.path.abspath( __file__ ) ) )
    write = not ( flags & { "--dry-run", "--check" } )
    updated, current, unrecognised, foreign = [], [], [], 0
    for path in candidates( root ):
        raw = open( path, "rb" ).read()
        try:
            text = raw.decode( "utf-8" )
        except UnicodeDecodeError:
            continue
        nl = "\r\n" if "\r\n" in text else "\n"
        lines = text.split( nl )
        head = lines[ :HEADER_LINES ]
        rel = os.path.relpath( path, root )
        if COPYRIGHT not in head:
            foreign += 1                 # not the project's file (or no header): left alone
            continue
        if NEW[ 0 ] in head:
            current.append( rel ); continue
        hits = [ i for i, l in enumerate( head ) if OLD.match( l ) ]
        if len( hits ) != 1:
            unrecognised.append( rel ); continue
        i = hits[ 0 ]
        lines[ i:i + 1 ] = NEW
        updated.append( rel )
        if write:
            open( path, "wb" ).write( nl.join( lines ).encode( "utf-8" ) )
    verb = "updated" if write else "to update"
    print( "update_licence_header.py: %d %s, %d already current, %d without the project's copyright line (left alone)"
           % ( len( updated ), verb, len( current ), foreign ) )
    if flags & { "--dry-run" }:
        for r in updated: print( "  " + r )
    for r in unrecognised:
        print( "  UNRECOGNISED licence line (left alone): " + r )
    if "--check" in flags:
        return 1 if ( updated or unrecognised ) else 0
    return 1 if unrecognised else 0

sys.exit( main() )
