# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
PATH_MC = $(shell cd ../../ ; pwd)
LIB_MC = -llapack -lblas
INC_MC = -I$(PATH_MC)/include
include $(PATH_MC)/src/makeoptions.mk

FLAG_ALL = $(FLAGS_CPP)
LIB_ALL = $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_SDPA) $(LIB_LAPACK)
INC_ALL = $(INC_MC) $(INC_FADBAD) $(INC_SDPA) $(INC_LAPACK) $(INC_PROFIL) $(INC_FILIB)

#####

all: test1 test2

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

