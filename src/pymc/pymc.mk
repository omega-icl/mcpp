# This makefile compiles the python interface and creates a symbolic link
# to the library in $(libpath)

include $(srcpath)/makeoptions.mk

#####

libobjs = mcfunc.o interval.o mccormick.o supmodel.o \
          ffunc.o fflin.o ffvect.o ffcustom.o \
          pymc.o
libname = pymc.so

#####

install: dispBuild $(libname) dispInstall
	@if test ! -e $(libpath)/$(libname); then \
		echo creating symolic link to shared library $(libname); \
		cd $(libpath); ln -s $(pymcpath)/$(libname) $(libname); \
	fi
	@echo

$(libname): $(libobjs)
	$(CPP) -shared -Wl,--export-dynamic $(libobjs) -o $(libname)

%.o : %.cpp
	$(CPP) $(INC_MC) $(INC_PYBIND11) $(FLAG_CPP) $(FLAG_MC) -fPIC -c $< -o $@

dispBuild:
	@echo
	@(echo '***Compiling PYMC library (ver.' $(version)')***')
	@echo

dispInstall:
	@echo
	@(echo '***Installing PYMC library (ver.' $(version)')***')
	@echo

#####

clean: dispClean
	rm -f $(libobjs) $(libname)

dispClean:
	@echo
	@(echo '***Cleaning PYMC directory (ver.' $(version)')***')
	@echo

#####

uninstall: dispUninstall
	cd $(libpath); rm -f $(libname)

dispUninstall:
	@echo
	@(echo '***Uninstalling PYMC library (ver.' $(version)')***')
	@echo
