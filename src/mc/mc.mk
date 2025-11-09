# This makefile creates symbolic links to the header files in $(incpath)

include $(srcpath)/makeoptions.mk

#####

incobjs = mcfunc.hpp mctime.hpp mclapack.hpp \
          mcop.hpp mcboost.hpp mcprofil.hpp mcfilib.hpp mcfadbad.hpp \
	  interval.hpp mccormick.hpp specbnd.hpp \
	  supmodel.hpp pwcu.hpp pwlu.hpp \
	  ellipsoid.hpp ellimage.hpp polimage.hpp \
	  polymodel.hpp tmodel.hpp cmodel.hpp smon.hpp scmodel.hpp sicmodel.hpp \
	  spoly.hpp squad.hpp \
	  ffdep.hpp ffinv.hpp ffunc.hpp ffexpr.hpp mchsl.hpp sred.hpp slift.hpp selim.hpp \
	  fflin.hpp ffspol.hpp ffvect.hpp ffextern.hpp ffdagext.hpp ffmlp.hpp ffmlpreg.hpp ffcustom.hpp

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

uninstall: dispUninstall
	cd $(incpath); rm -f $(incobjs)

dispUninstall:
	@echo
	@(echo '***Uninstalling MC++ library (ver.' $(version)')***')
	@echo
