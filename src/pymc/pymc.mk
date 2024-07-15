# This makefile compiles the python interface and creates a symbolic link
# to the libray in $(libpath)

include $(srcpath)/makeoptions.mk

#####

libobjs = interval.o ffunc.o main.o
libname = pymc.so

#####

install: dispBuild $(libname) dispInstall
	@if test ! -e $(libpath)/$(libname); then \
		echo creating symolic link to shared library $(libname); \
		cd $(libpath) ; ln -s $(pymcpath)/$(libname) $(libname); \
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
	rm -fi $(libobjs) $(libname)

dispClean:
	@echo
	@(echo '***Cleaning PYMC directory (ver.' $(version)')***')
	@echo

#####

cleandist: dispCleanInstall
	rm -f $(libname)
	-(cd $(libpath) ; rm -f $(libname))

dispCleanInstall:
	@echo
	@(echo '***Uninstalling PYMC library (ver.' $(version)')***')
	@echo

