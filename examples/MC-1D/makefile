# Executable name
EXE = MC-1D

PATH_MC = $(shell cd ../../ ; pwd)
LIB_MC = -L$(PATH_MC)/lib $(LIB_PROFIL) $(LIB_FILIB) $(LIB_LAPACK) $(LIB_FADBAD)
INC_MC = -I$(PATH_MC)/include $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_FADBAD)
OBJS = main.o

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

$(EXE) : $(OBJS)
	$(LINK) $(PROF) $(FLAGS_LINK) -o $(EXE) $(OBJS) $(LIB_MC)
       
main.o: main.cpp
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_MC) -o main.o main.cpp

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
	rm -f $(EXE) $(OBJS) *.out *.eps *.dot *.svg *~
