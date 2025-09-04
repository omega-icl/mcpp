# This is the main makefile for MC++
# Run `make help' for a list of available targets

export version = 4.0

export rootpath = $(shell cd .. ; pwd)
export srcpath  = $(rootpath)/src
export incpath  = $(rootpath)/include
export libpath  = $(rootpath)/lib
export mcpath   = $(srcpath)/mc
export pymcpath = $(srcpath)/pymc

#####

install: dirs_create mc_install pymc_install

dirs_create:
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

clean: pymc_clean

pymc_clean:
	-(cd $(pymcpath); make -f pymc.mk clean)

#####

uninstall: mc_uninstall pymc_uninstall

mc_uninstall:
	-(cd $(mcpath); make -f mc.mk uninstall)

pymc_uninstall:
	-(cd $(pymcpath); make -f pymc.mk uninstall)

#####

help:
	@echo
	@echo "MC++ Makefile (ver. $(version))"
	@echo "Available targets:"
	@echo "  install      - Create 'include' and 'lib' directories, compile a shared library,"
	@echo "                 create symbolic links to header files and compiled libraries"
	@echo "  clean        - Remove object files"
	@echo "  uninstall    - Remove libraries, object files, and symbolic links"
	@echo "  help         - Show this help message"
	@echo
