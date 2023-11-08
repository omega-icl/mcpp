PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

PATH_GUROBI = $(GUROBI_HOME)
LIB_MIP     = -L$(PATH_GUROBI)/lib -lgurobi_g++5.2 -lgurobi100 -pthread
INC_MIP     = -I$(PATH_GUROBI)/include
FLAG_MIP    = -DMC__USE_GUROBI

#####

all: LINREL_build LINREL_build2 LINREL_build3 LINREL_1d LINREL_2d

LINREL_% : LINREL_%.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC) $(LIB_MIP)

PWLREL_% : PWLREL_%.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC) $(LIB_MIP)

%.o: %.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $(FLAG_MIP) $(INC_MIP) $< -o $@

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
	rm -f LINREL_build LINREL_build2 LINREL_build3 LINREL_?d PWLREL_?d
	rm -f *.o *~ *.eps *.lp *.dot *.svg *.out

