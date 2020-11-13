PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

LIB_MC = -L$(PATH_MC)/lib $(LIB_PROFIL) $(LIB_FILIB) $(LIB_LAPACK) $(LIB_HSL)
#INC_MC = -I$(PATH_MC)/include $(INC_FADBAD) $(INC_LAPACK) $(INC_PROFIL) $(INC_FILIB)
INC_MC = -I$(PATH_MC)/src/mc $(INC_FADBAD) $(INC_LAPACK) $(INC_PROFIL) $(INC_FILIB) $(INC_HSL)

#####

all: TM-1D TM-2D CM-1D CM-2D CM-SCM-2D SCM-QUAD-2D

TM-1D : TM-1D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_MC)

TM-2D : TM-2D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_MC)

CM-1D : CM-1D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_MC)

CM-2D : CM-2D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_MC)

CM-SCM-2D : CM-SCM-2D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_MC)

CM-SCM-3D : CM-SCM-3D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_MC)

SCM-QUAD-2D : SCM-QUAD-2D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_MC)

%.o: %.cpp 
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_MC) $< -o $@

#####

graph : gendot gensvg

gendot :
	dot -Tsvg -O $(shell ls *.dot)

gensvg :
	@for FILE in $(shell ls *.svg); do \
	  echo $$FILE; \
	  (display $$FILE &); \
	done

#####

clean :
	rm -f TM-1D TM-2D CM-1D CM-2D CM-SCM-2D SCM-QUAD-2D *.o *.out *.eps *.dot *.svg *~

