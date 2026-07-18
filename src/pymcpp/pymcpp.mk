# This makefile compiles the python interface and creates a symbolic link
# to the library in $(libpath)

include $(srcpath)/makeoptions.mk

#####

libobjs = mcfunc.o interval.o mccormick.o specbnd.o \
          tmodel.o cmodel.o scmodel.o sicmodel.o \
          supmodel.o polimage.o ellimage.o \
          smon.o spoly.o \
          ffdep.o ffinv.o ffunc.o ffmon.o ffpoly.o slift.o \
          fflin.o ffmlp.o ffdagext.o ffvect.o ffcustom.o \
          pymcpp.o
libname = pymcpp.so

#####

install: dispBuild $(libname) dispInstall
	@if test ! -e $(libpath)/$(libname); then \
		echo creating symolic link to shared library $(libname); \
		cd $(libpath); ln -s $(pymcpppath)/$(libname) $(libname); \
	fi
	@echo

$(libname): $(libobjs)
	$(CPP) -shared -Wl,--export-dynamic $(libobjs) $(LIB_MC) -o $(libname)

%.o : %.cpp
	$(CPP) $(INC_PYBIND11) $(INC_MC) $(FLAG_CPP) $(FLAG_MC) -fPIC -c $< -o $@

dispBuild:
	@echo
	@(echo '***Compiling PYMCPP library (ver.' $(version)')***')
	@echo

dispInstall:
	@echo
	@(echo '***Installing PYMCPP library (ver.' $(version)')***')
	@echo

#####

clean: dispClean
	rm -fi $(libobjs) $(libname)

dispClean:
	@echo
	@(echo '***Cleaning PYMCPP directory (ver.' $(version)')***')
	@echo

#####

uninstall: dispUninstall
	rm -f $(libobjs) $(libname)
	cd $(libpath); rm -f $(libname)

dispUninstall:
	@echo
	@(echo '***Uninstalling PYMCPP library (ver.' $(version)')***')
	@echo
