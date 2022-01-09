PATH_MC = $(shell cd ../../ ; pwd)

# Compilation options
include $(PATH_MC)/src/makeoptions.mk

#####

all: QUAD-TEST SPARSEEXPR-TEST

QUAD-TEST : QUAD-TEST.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

SPARSEEXPR-TEST : SPARSEEXPR-TEST.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_MC)

%.o: %.cpp 
	$(CPP) -c $(FLAG_CPP) $(FLAG_MC) $(INC_MC) $< -o $@

#####

clean :
	rm -f QUAD-TEST SPARSEEXPR-TEST *.o *.out *.eps *.dot *.svg *~

