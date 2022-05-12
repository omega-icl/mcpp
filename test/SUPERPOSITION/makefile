PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

LIB_MC = -L$(PATH_MC)/lib $(LIB_PROFIL) $(LIB_FILIB) $(LIB_LAPACK) $(LIB_HSL)
#INC_MC = -I$(PATH_MC)/include $(INC_FADBAD) $(INC_LAPACK) $(INC_PROFIL) $(INC_FILIB)
INC_MC = -I$(PATH_MC)/src/mc $(INC_FADBAD) $(INC_LAPACK) $(INC_PROFIL) $(INC_FILIB) $(INC_HSL)

#####

all: test1 test2

test1 : test1.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_MC)

test2 : test2.o
	$(LINK) $(PROF) $(FLAGS_LINK) -o $@ $< $(LIB_MC)

%.o: %.cpp 
	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_MC) $< -o $@

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
	rm -f test1 *.o *.out *.eps *.dot *.svg *~



## Executable name
#EXE = test

#PATH_MC = $(shell cd ../../ ; pwd)
## Compilation options
#include $(PATH_MC)/src/makeoptions.mk

#LIB_MC = -L$(PATH_MC)/lib $(LIB_PROFIL) $(LIB_BOOST) $(LIB_LAPACK) $(LIB_HSL)
#INC_MC = -I$(PATH_MC)/src/mc $(INC_FADBAD) $(INC_LAPACK) $(INC_PROFIL) $(INC_BOOST) $(INC_HSL)
##LIB_MC = -L$(PATH_MC)/lib $(LIB_PROFIL) $(LIB_FILIB) $(LIB_LAPACK) $(LIB_FADBAD)
##INC_MC = -I$(PATH_MC)/include $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_FADBAD)
#OBJS = test1.o


######

#$(EXE) : $(OBJS)
#	echo $(PATH_MC)
#	$(LINK) $(PROF) $(FLAGS_LINK) -o $(EXE) $(OBJS) $(LIB_MC)
#       
#test1.o: test1.cpp
#	echo $(PATH_MC)
#	$(CPP) -c $(PROF) $(FLAGS_CPP) $(INC_MC) -o test1.o test1.cpp

######

#graph : gendot gensvg

#gendot :
#	dot -Tsvg -O $(shell ls *.dot)

#gensvg :
#	@for FILE in $(shell ls *.svg); do \
#	  echo $$FILE; \
#	  (display $$FILE &); \
#	done

######

#clean :
#	rm -f $(EXE) $(OBJS) *.out *.eps *.dot *.svg *~
#	

