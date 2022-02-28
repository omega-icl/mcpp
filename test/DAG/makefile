PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

all: DAG_doxygen DAG_compose DAG_differentiate DAG_evaluate DAG_analyze

DAG_doxygen : DAG_doxygen.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

DAG_compose : DAG_compose.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

DAG_differentiate : DAG_differentiate.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

DAG_evaluate : DAG_evaluate.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

DAG_analyze : DAG_analyze.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

%.o: %.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $< -o $@

#####

graph : gendot gensvg

gendot :
	dot -Tpng -O $(shell ls *.dot)

gensvg :
	@for FILE in $(shell ls *.png); do \
	  echo $$FILE; \
	  (display $$FILE &); \
	done

#####

clean :
	rm -f DAG_compose DAG_doxygen DAG_differentiate DAG_evaluate DAG_analyze *.o *.out *.eps *.dot *.png *.lp *~

