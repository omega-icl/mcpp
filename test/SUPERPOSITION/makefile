PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

PATH_GUROBI = $(GUROBI_HOME)
LIB_GUROBI  = -L$(PATH_GUROBI)/lib -lgurobi_g++5.2 -lgurobi91 -pthread
INC_GUROBI  = -I$(PATH_GUROBI)/include

#####

all: test1 test2 test3

test1 : test1.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test2 : test2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

test3 : test3.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC) $(LIB_GUROBI)

testAM : testAM.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

%.o: %.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $(INC_GUROBI) $< -o $@

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
	rm -f test? *.o *.out *.eps *.dot *.svg *~

