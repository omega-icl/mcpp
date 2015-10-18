PATH_MC = $(shell cd ../../ ; pwd)
LIB_MC = -L$(PATH_MC)/lib
INC_MC = -I$(PATH_MC)/include

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

test: test1 test2

test1: test1.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o test1 test1.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_LAPACK)
       
test1.o: test1.cpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_MC) $(INC_FADBAD) $(INC_LAPACK) $(INC_PROFIL) $(INC_FILIB) -o test1.o test1.cpp

test2: test2.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o test2 test2.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_LAPACK)
       
test2.o: test2.cpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_MC) $(INC_FADBAD) $(INC_LAPACK) $(INC_PROFIL) $(INC_FILIB) -o test2.o test2.cpp

#####

MC:
	(cd $(PATH_MC)/src; make; make install)

#####

clean:
	rm -f test? test?.o *.out *.eps
