# Executable name
EXE = TEST

PATH_MC = $(shell cd ../../ ; pwd)
INC_MC = -I$(PATH_MC)/include
OBJS = main.o

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

test : $(EXE) dispTest
	./$(EXE)

dispTest:
	@echo
	@(echo '***Testing MC++ library (ver.' $(version)')***')
	@echo

$(EXE) : $(OBJS)
	$(LINK) $(PROF) $(FLAGS_LINK) -o $(EXE) $(OBJS) $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_LAPACK) $(LIB_CPPUNIT)
      
main.o: main.cpp interval_test.hpp mccormick_test.hpp tmodel_test.hpp specbnd_test.hpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_FADBAD) -o main.o main.cpp

#####

clean :
	rm -f $(EXE) $(OBJS)
