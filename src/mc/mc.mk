# This makefile creates symbolic links to the header files in $(incpath)

include $(srcpath)/makeoptions.mk

#####

incobjs = mcfunc.hpp mctime.hpp mclapack.hpp \
          mcop.hpp mcboost.hpp mcprofil.hpp mcfilib.hpp mcfadbad.hpp \
	  interval.hpp ismodel.hpp mccormick.hpp specbnd.hpp \
	  ellipsoid.hpp ellimage.hpp polimage.hpp \
          polymodel.hpp tmodel.hpp cmodel.hpp smon.hpp scmodel.hpp sicmodel.hpp \
          spoly.hpp squad.hpp \
          univarmodels.hpp amsmodel.hpp \
          ffdep.hpp ffinv.hpp ffunc.hpp sred.hpp slift.hpp selim.hpp

#####

install: dispInstall
	@for INC in $(incobjs); do \
		if test ! -e $(incpath)/$$INC; then \
			echo creating symbolic link to header file $$INC; \
			cd $(incpath); ln -s $(mcpath)/$$INC $$INC; \
		fi; \
	done
	@echo

dispInstall:
	@echo
	@(echo '***Installing MC++ library (ver.' $(version)')***')
	@echo

#####

cleaninstall: dispCleanInstall
	cd $(incpath) ; rm $(incobjs)

dispCleanInstall:
	@echo
	@(echo '***Uninstalling MC++ library (ver.' $(version)')***')
	@echo
