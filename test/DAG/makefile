PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

PATH_GUROBI = $(GUROBI_HOME)
LIB_MIP     = -L$(PATH_GUROBI)/lib -lgurobi_g++5.2 -lgurobi100 -pthread
INC_MIP     = -I$(PATH_GUROBI)/include
FLAG_MIP    = -DMC__USE_GUROBI

#####

all: DAG_doxygen DAG_compose DAG_differentiate DAG_evaluate DAG_external DAG_store DAG_analyze

DAG_doxygen : DAG_doxygen.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

DAG_compose : DAG_compose.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

DAG_differentiate : DAG_differentiate.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

DAG_evaluate : DAG_evaluate.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

DAG_external : DAG_external.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

DAG_store : DAG_store.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

DAG_analyze : DAG_analyze.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC) $(LIB_MIP)

DAG_analyze.o: DAG_analyze.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $(FLAG_MIP) $(INC_MIP) $< -o $@

%.o: %.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $< -o $@

#####

graph : gendot gensvg

gendot :
	dot -Tpng -O $(shell ls *.dot)

gensvg :
	@for FILE in $(shell ls *.png); do \
	  echo $$FILE; \
      convert -resize 70% $$FILE $$FILE; \
	  (display $$FILE &); \
	done

#####

clean :
	rm -f DAG_compose DAG_doxygen DAG_differentiate DAG_evaluate DAG_store DAG_external DAG_analyze
	rm -f *.o *.out *.eps *.dot *.png *.lp *~

