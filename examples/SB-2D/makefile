PATH_MC = $(shell cd ../../ ; pwd)
LIB_MC = -L$(PATH_MC)/lib $(LIB_PROFIL) $(LIB_FILIB) $(LIB_LAPACK)
INC_MC = -I$(PATH_MC)/include $(INC_FADBAD) $(INC_LAPACK) $(INC_PROFIL) $(INC_FILIB)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

test: test1 test2

test% : test%.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_MC)

test%.o: test%.cpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_MC) -o $@ $<

#####

clean:
	rm -f test? *.o *.out *.eps *~
