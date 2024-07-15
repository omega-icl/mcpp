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
export pymcpath = $(srcpath)/pymc

#####

all: mc_all pymc_all

mc_all:
	-(cd $(mcpath); make -f mc.mk all)

pymc_all:
	-(cd $(pymcpath); make -f pymc.mk all)

#####

install: dirs mc_install pymc_install

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

pymc_install:
	-(cd $(pymcpath); make -f pymc.mk install)

#####

clean: mc_clean pymc_clean

mc_clean:
	-(cd $(mcpath);   make -f mc.mk clean)

pymc_clean:
	-(cd $(pymcpath); make -f pymc.mk clean)

#####

cleaninstall: clean mc_cleaninstall pymc_cleaninstall

mc_cleaninstall:
	-(cd $(mcpath);   make -f mc.mk cleaninstall)
	
pymc_cleaninstall:
	-(cd $(pymcpath); make -f pymc.mk cleaninstall)

