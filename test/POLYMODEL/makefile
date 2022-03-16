PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

all: TM-1D TM-2D CM-1D CM-2D CM-SCM-2D SCM-SICM-2D SCM-QUAD-2D

TM-1D : TM-1D.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

TM-2D : TM-2D.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

CM-1D : CM-1D.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

CM-2D : CM-2D.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

CM-SCM-2D : CM-SCM-2D.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

SCM-SICM-2D : SCM-SICM-2D.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

SCM-QUAD-2D : SCM-QUAD-2D.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

%.o: %.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $< -o $@

#####

clean :
	rm -f QUAD-TEST SPARSEEXPR-TEST *.o *.out *.eps *.dot *.svg *~

