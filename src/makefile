# This is the main makefile for MC++
# - `make install' creates the subsdirectory `include' (if not already
#   existing, and creates links toheader files therein
# - `make cleaninstall' removes the links and cleans src/test

export version = 2.2

export rootpath = $(shell cd .. ; pwd)
export srcpath  = $(rootpath)/src
export incpath  = $(rootpath)/include
export mcpath   = $(srcpath)/mc

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

cleaninstall: mc_cleaninstall

mc_cleaninstall:
	-(cd $(mcpath);   make -f mc.mk cleaninstall)
