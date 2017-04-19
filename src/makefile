# This is the main makefile for MC++
# - `make install' creates the subsdirectory `include' (if not already
#   existing, and creates links toheader files therein
# - `make test' runs unit testing for various components of MC++ using
#   CppUnit
# - `make cleaninstall' removes the links and cleans src/test

export version = 2.0

export rootpath = $(shell cd .. ; pwd)
export srcpath  = $(rootpath)/src
export incpath  = $(rootpath)/include
export mcpath   = $(srcpath)/mc
export testpath = $(srcpath)/test

#####

install: dirs mc_install

dirs :
	@if test ! -e $(incpath)/; then \
		echo creating header files directory; \
		mkdir $(incpath); \
	fi
	@if test ! -e $(libpath)/; then \
		echo creating header files directory; \
		mkdir $(libpath); \
	fi

mc_install:
	-(cd $(mcpath); make -f mc.mk install)

#####

test: mc_test

mc_test:
	-(cd $(testpath); make -f test.mk)

#####

cleaninstall: mc_cleaninstall

mc_cleaninstall:
	-(cd $(mcpath);   make -f mc.mk cleaninstall)
	-(cd $(testpath); make -f test.mk clean)
