# This is the main makefile for MC++
# - `make all' compiles the reates the subsdirectory `include' (if not already
#   existing, and creates links toheader files therein
# - `make install' creates the subsdirectories `include' and 'lib' (if not already
#   existing), and creates smybolic links therein
# - `make cleaninstall' removes the smybolic links and cleans src/test

export version = 3.0

export rootpath = $(shell cd .. ; pwd)
export srcpath  = $(rootpath)/src
export incpath  = $(rootpath)/include
export libpath  = $(rootpath)/lib
export mcpath   = $(srcpath)/mc
export mcpypath = $(srcpath)/mcpy

#####

all: mc_all mcpy_all

mc_all:
	-(cd $(mcpath); make -f mc.mk all)

mcpy_all:
	-(cd $(mcpypath); make -f mcpy.mk all)

#####

install: dirs mc_install mcpy_install

dirs :
	@if test ! -e $(incpath)/; then \
		echo creating include directory; \
		mkdir $(incpath); \
	fi
	@if test ! -e $(libpath)/; then \
		echo creating lib directory; \
		mkdir $(libpath); \
	fi

mc_install:
	-(cd $(mcpath); make -f mc.mk install)

mcpy_install:
	-(cd $(mcpypath); make -f mcpy.mk install)

#####

clean: mc_clean mcpy_clean

mc_clean:
	-(cd $(mcpath);   make -f mc.mk clean)

mcpy_clean:
	-(cd $(mcpypath); make -f mcpy.mk clean)

#####

cleaninstall: clean mc_cleaninstall mcpy_cleaninstall

mc_cleaninstall:
	-(cd $(mcpath);   make -f mc.mk cleaninstall)
	
mcpy_cleaninstall:
	-(cd $(mcpypath); make -f mcpy.mk cleaninstall)

