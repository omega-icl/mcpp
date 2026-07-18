# This is the main makefile for MC++
# Run `make help' for a list of available targets

export version = 5.0

export rootpath = $(shell cd .. ; pwd)
export srcpath  = $(rootpath)/src
export incpath  = $(rootpath)/include
export libpath  = $(rootpath)/lib
export mcpath   = $(srcpath)/mc
export pymcpppath = $(srcpath)/pymcpp

#####

install: dirs_create mc_install pymcpp_install

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

pymcpp_install:
	-(cd $(pymcpppath); make -f pymcpp.mk install)

#####

clean: pymcpp_clean

pymcpp_clean:
	-(cd $(pymcpppath); make -f pymcpp.mk clean)

#####

uninstall: mc_uninstall pymcpp_uninstall

mc_uninstall:
	-(cd $(mcpath); make -f mc.mk uninstall)

pymcpp_uninstall:
	-(cd $(pymcpppath); make -f pymcpp.mk uninstall)

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
