PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

PATH_GUROBI = $(GUROBI_HOME)
LIB_MIP     = -L$(PATH_GUROBI)/lib -lgurobi_g++5.2 -lgurobi100 -pthread
INC_MIP     = -I$(PATH_GUROBI)/include
FLAG_MIP    = -DMC__USE_GUROBI

#####

all: SQUAD-TEST SRED-TEST

SQUAD-TEST : SQUAD-TEST.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC) $(LIB_MIP)

SQUAD-TEST.o: SQUAD-TEST.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $(FLAG_MIP) $(INC_MIP) $< -o $@

SRED-TEST : SRED-TEST.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC) $(LIB_MIP)

SRED-TEST.o: SRED-TEST.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $(FLAG_MIP) $(INC_MIP) $< -o $@

%.o: %.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $< -o $@

#####

clean :
	rm -f SQUAD-TEST SRED-TEST
	rm -f *.o *.out *.eps *.dot *.svg *.lp *~

