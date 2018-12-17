PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#FLAGS_OPT = -Ofast

LIB_MC = -L$(PATH_MC)/lib $(LIB_FADBAD) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_LAPACK) $(LIB_HSL)
INC_MC = -I$(PATH_MC)/src/mc $(INC_FADBAD) $(INC_LAPACK) $(INC_PROFIL) $(INC_FILIB) $(INC_HSL)
FLAGS_MC = $(FLAGS_FILIB) $(FLAGS_HSL)

PATH_GUROBI = $(GUROBI_HOME)
LIB_GUROBI  = -L$(PATH_GUROBI)/lib -lgurobi_g++5.2 -lgurobi81 -pthread
INC_GUROBI  = -I$(PATH_GUROBI)/include

FLAGS_ALL = $(FLAGS_CPP) $(FLAGS_OPT) $(FLAGS_MC) -DMC__USE_GUROBI
LIB_ALL   = $(LIB_MC) $(LIB_GUROBI)
INC_ALL   = $(INC_MC) $(INC_GUROBI)

#####

all: LINREL_build LINREL_1d LINREL_2d

LINREL_% : LINREL_%.o
	g++ -o $@ $< $(LIB_ALL)

LINREL_%.o: LINREL_%.cpp
	g++ -c $(FLAGS_ALL) -I../../src $(INC_ALL) -o $@ $<

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
	rm -f LINREL_build LINREL_?d *.o *~ *.eps *.lp *.dot *.svg *.out

