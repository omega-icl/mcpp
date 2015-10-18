# Executable name
EXE = POLTEST

PATH_MC = $(shell cd ../../ ; pwd)
LIB_MC = -L$(PATH_MC)/lib
INC_MC = -I$(PATH_MC)/include

PATH_GUROBI = $(GUROBI_HOME)
LIB_GUROBI = -L$(PATH_GUROBI)/lib -lgurobi_c++ -lgurobi56 -pthread
INC_GUROBI = -I$(PATH_GUROBI)/include

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

OBJS = test.o

#####

$(EXE) : $(OBJS)
	$(LINK) $(PROF) $(FLAGS_LINK) -o $(EXE) $(OBJS) $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_LAPACK) $(LIB_GUROBI)

test.o: test.cpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_MC) $(INC_FADBAD) $(INC_LAPACK) $(INC_PROFIL) $(INC_FILIB) $(INC_GUROBI) -o test.o test.cpp

test2D.o: test2D.cpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_MC) $(INC_FADBAD) $(INC_LAPACK) $(INC_PROFIL) $(INC_FILIB) $(INC_GUROBI) -o test2D.o test2D.cpp

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
	rm -f $(EXE) $(OBJS) *.dot *.svg *.out *.eps *.lp *~

