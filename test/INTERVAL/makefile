PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

LIB_MC = -L$(PATH_MC)/lib $(LIB_PROFIL) $(LIB_BOOST) $(LIB_LAPACK) $(LIB_HSL)
INC_MC = -I$(PATH_MC)/src/mc $(INC_FADBAD) $(INC_LAPACK) $(INC_PROFIL) $(INC_BOOST) $(INC_HSL)

#####

all: INT-1D INT-2D

INT-?D : INT-?D.o
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
	rm -f INT-?D *.o *.out *.eps *.dot *.svg *~

