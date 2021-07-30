PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

LIB_MC = -L$(PATH_MC)/lib $(LIB_FADBAD) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_LAPACK) $(LIB_HSL)
INC_MC = -I$(PATH_MC)/src/mc $(INC_FADBAD) $(INC_LAPACK) $(INC_BOOST) $(INC_PROFIL) $(INC_FILIB) $(INC_HSL)
FLAGS_MC = $(FLAGS_FILIB) $(FLAGS_HSL)

#FLAGS_OPT = -Ofast

#PATH_CPLEX   = /opt/ibm/ILOG/CPLEX_Studio1210/cplex
#PATH_CONCERT = /opt/ibm/ILOG/CPLEX_Studio1210/concert
#LIB_CPLEX    = -L$(PATH_CPLEX)/lib/x86-64_linux/static_pic -lilocplex -lcplex \
#               -L$(PATH_CONCERT)/lib/x86-64_linux/static_pic -lconcert \
#               -lm -pthread -ldl
#INC_CPLEX    = -I$(PATH_CPLEX)/include -I$(PATH_CONCERT)/include
#FLAGS_CPLEX  = -m64 -fPIC -fexceptions -DIL_STD -Wno-ignored-attributes
##FLAGS_CPLEX  = -m64 -fPIC -fexceptions -DNOEBUG -DIL_STD

#FLAGS_ALL = $(FLAGS_CPP) $(FLAGS_OPT) $(FLAGS_MC) -DMC__USE_CPLEX $(FLAGS_CPLEX)
#LIB_ALL   = $(LIB_MC) $(LIB_CPLEX)
#INC_ALL   = $(INC_MC) $(INC_CPLEX)

PATH_GUROBI = $(GUROBI_HOME)
LIB_GUROBI  = -L$(PATH_GUROBI)/lib -lgurobi_g++5.2 -lgurobi91 -pthread
INC_GUROBI  = -I$(PATH_GUROBI)/include

FLAGS_ALL = $(FLAGS_CPP) $(FLAGS_OPT) $(FLAGS_MC) -DMC__USE_GUROBI
LIB_ALL   = $(LIB_MC) $(LIB_GUROBI)
INC_ALL   = $(INC_MC) $(INC_GUROBI)

#####

all: TM-1D TM-2D CM-1D CM-2D CM-SCM-2D SCM-QUAD-2D

TM-1D : TM-1D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

TM-2D : TM-2D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

CM-1D : CM-1D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

CM-2D : CM-2D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

SIRM-1D : SIRM-1D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

SIRM-2D : SIRM-2D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

CM-SCM-2D : CM-SCM-2D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

CM-SCM-3D : CM-SCM-3D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

SCM-SICM-2D : SCM-SICM-2D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

SCM-QUAD-2D : SCM-QUAD-2D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

SCM-QUAD-6D : SCM-QUAD-6D.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_ALL)

%.o: %.cpp 
	$(CPP) -c $(PROF) $(FLAGS_ALL) $(INC_ALL) -o $@ $<
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
	rm -f TM-1D TM-2D CM-1D CM-2D SIRM-1D SIRM-2D CM-SCM-2D SCM-SICM-2D SCM-QUAD-2D *.o *.out *.eps *.dot *.svg *~

