PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

LIB_MC = -L$(PATH_MC)/lib $(LIB_FADBAD) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_LAPACK) $(LIB_HSL)
INC_MC = -I$(PATH_MC)/src/mc $(INC_FADBAD) $(INC_LAPACK) $(INC_PROFIL) $(INC_FILIB) $(INC_HSL)
FLAGS_MC = $(FLAGS_FILIB) $(FLAGS_HSL)

#FLAGS_OPT = -Ofast

PATH_CPLEX   = /opt/ibm/ILOG/CPLEX_Studio128/cplex
PATH_CONCERT = /opt/ibm/ILOG/CPLEX_Studio128/concert
LIB_CPLEX    = -L$(PATH_CPLEX)/lib/x86-64_linux/static_pic -lilocplex -lcplex \
               -L$(PATH_CONCERT)/lib/x86-64_linux/static_pic -lconcert \
               -lm -pthread -ldl
INC_CPLEX    = -I$(PATH_CPLEX)/include -I$(PATH_CONCERT)/include
FLAGS_CPLEX  = -m64 -fPIC -fexceptions -DIL_STD -Wno-ignored-attributes
#FLAGS_CPLEX  = -m64 -fPIC -fexceptions -DNOEBUG -DIL_STD

FLAGS_ALL = $(FLAGS_CPP) $(FLAGS_OPT) $(FLAGS_MC) -DMC__USE_CPLEX $(FLAGS_CPLEX)
LIB_ALL   = $(LIB_MC) $(LIB_CPLEX)
INC_ALL   = $(INC_MC) $(INC_CPLEX)

#PATH_GUROBI = $(GUROBI_HOME)
#LIB_GUROBI  = -L$(PATH_GUROBI)/lib -lgurobi_g++5.2 -lgurobi75 -pthread
#INC_GUROBI  = -I$(PATH_GUROBI)/include

#FLAGS_ALL = $(FLAGS_CPP) $(FLAGS_OPT) $(FLAGS_MC) -DMC__USE_GUROBI
#LIB_ALL   = $(LIB_MC) $(LIB_GUROBI)
#INC_ALL   = $(INC_MC) $(INC_GUROBI)

#####

all: DAG_compose DAG_differentiate DAG_evaluate DAG_analyze DAG_polyexpr

DAG_compose : DAG_compose.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

DAG_doxygen : DAG_doxygen.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

DAG_differentiate : DAG_differentiate.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

DAG_evaluate : DAG_evaluate.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

DAG_analyze : DAG_analyze.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

DAG_polyexpr : DAG_polyexpr.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

test_profil : test_profil.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

%.o: %.cpp 
	$(CPP) -c $(PROF) $(FLAGS_ALL) $(INC_ALL) -o $@ $<

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

