CONTENTS
     inlude\: CPPLapack header files
        lib\: BLAS and LAPACK libraries for 32bit Windows
      lib64\: BLAS and LAPACK libraries for 64bit Windows
  vc_sample\: A sample project for Visual C++ 2008 Express (vcproj)
              and Visual C++ 2010 Express (vcxproj)
  README.txt: This file

UNINSTALL
  Simplly remove the "c:\cpplapack" directory.

LIBRARIES
  The BLAS and LAPACK libraries are NOT products of CPPLapack.
  They are originally distributed at http://www.netlib.org/clapack/ .
  Authors of CPPLapack deeply appliciate the developers of BLAS and LAPACK.

USAGE
  Add "c:\cpplapack\include" to the include path of your project.
  You may also need to add either "c:\cpplapack\lib" or "c:\cpplapack\lib64"
  to the library path of your project and then link with the "libf2c.lib", 
  "BLAS.lib", and "clapack.lib".

DOCUMENTATION
  An online documentation is available at http://cpplapack.sourceforge.net/ .

LICENCE
  CPPLapack is a GPL software.
