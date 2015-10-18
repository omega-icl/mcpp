# Executable name
EXE = EIMGTEST

PATH_MC = $(shell cd ../../ ; pwd)
LIB_MC = -L$(PATH_MC)/lib
INC_MC = -I$(PATH_MC)/include
OBJS = main.o

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

$(EXE) : $(OBJS)
	$(LINK) $(PROF) $(FLAGS_LINK) -o $(EXE) $(OBJS) $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_SDPA) $(LIB_LAPACK)
       
main.o: main.cpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_MC) $(INC_FADBAD) $(INC_SDPA) $(INC_LAPACK) $(INC_PROFIL) $(INC_FILIB) -o main.o main.cpp

#####

clean :
	rm -f $(EXE) $(OBJS) *.eps
