# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
PATH_MC = $(shell cd ../../ ; pwd)
LIB_MC = -llapack -lblas
INC_MC = -I$(PATH_MC)/include
include $(PATH_MC)/src/makeoptions.mk

PATH_GUROBI = $(GUROBI_HOME)
LIB_GUROBI = -L$(PATH_GUROBI)/lib -lgurobi_c++ -lgurobi56 -pthread
INC_GUROBI = -I$(PATH_GUROBI)/include

FLAG_ALL = $(FLAGS_CPP)
LIB_ALL  = $(LIB_MC) $(LIB_LAPACK) $(LIB_FADBAD) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_GUROBI)
INC_ALL  = $(INC_MC) $(INC_LAPACK) $(INC_FADBAD) $(INC_PROFIL) $(INC_FILIB) $(INC_GUROBI)

#####

all: test0 test1D test2D

test% : test%.o
	g++ -o $@ $< $(LIB_ALL)

test%.o: test%.cpp
	g++ -c $(FLAG_ALL) -I../../src $(INC_ALL) -o $@ $<

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
	rm -f test? test?D *.o *~ *.eps *.lp *.dot *.svg *.out

